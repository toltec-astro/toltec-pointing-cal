#!/usr/bin/env python
"""Compute calibration ratios from an ECSV pointing file and flux estimates.

Expected SQLite schema (default table name: flux_estimates):

    CREATE TABLE flux_estimates (
        source TEXT NOT NULL,
        obs_date TEXT NOT NULL,     -- ISO-ish datetime, e.g. 2026-02-04 08:25:06
        freq_ghz REAL NOT NULL,
        flux_mjy REAL NOT NULL,
        PRIMARY KEY (source, obs_date, freq_ghz)
    );

The script reads source/date from ECSV metadata and measured flux from the `amp` column,
then queries flux estimates for each target frequency.
Calibration ratio is defined as: measured_flux_mjy / estimated_flux_mjy.
"""

from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sqlite3
import sys
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Protocol
from urllib.parse import urlencode
from urllib.request import urlopen

from astropy.table import Table

TARGET_FREQUENCIES_GHZ = (150.0, 220.0, 280.0)
DEFAULT_ARRAY_MAP = {0: 280.0, 1: 220.0, 2: 150.0}


@dataclass
class CalibrationResult:
    source: str
    obs_date_utc: str
    frequency_ghz: float
    array_index: int
    measured_flux_mjy: float
    measured_flux_err_mjy: Optional[float]
    estimated_flux_mjy: float
    estimate_obs_date_utc: str
    calibration_ratio: float
    backend_used: str
    sma_1mm_date_delta_days: Optional[float]
    alma_nearest_obs_days: Optional[float]


@dataclass
class FluxEstimate:
    flux_mjy: float
    estimate_obs_date_utc: str
    backend: str
    sma_1mm_date_delta_days: Optional[float] = None
    alma_nearest_obs_days: Optional[float] = None


class FluxEstimator(Protocol):
    def get_estimate(self, source: str, obs_date_utc: datetime, freq_ghz: float) -> FluxEstimate:
        ...


class SQLiteFluxDatabase:
    def __init__(self, db_path: Path, table_name: str = "flux_estimates") -> None:
        self._db_path = Path(db_path)
        self._table_name = table_name
        self._conn = sqlite3.connect(str(self._db_path))

    def close(self) -> None:
        self._conn.close()

    def get_estimate(self, source: str, obs_date_utc: datetime, freq_ghz: float) -> FluxEstimate:
        query = f"""
            SELECT flux_mjy, obs_date
            FROM {self._table_name}
            WHERE source = ?
              AND CAST(freq_ghz AS REAL) = ?
            ORDER BY ABS(julianday(obs_date) - julianday(?))
            LIMIT 1
        """

        row = self._conn.execute(query, (source, float(freq_ghz), obs_date_utc.isoformat(sep=" "))).fetchone()
        if row is None:
            raise LookupError(
                f"No flux estimate found for source={source!r} at frequency={freq_ghz:g} GHz in table {self._table_name!r}."
            )
        return FluxEstimate(
            flux_mjy=float(row[0]),
            estimate_obs_date_utc=str(row[1]),
            backend="sqlite",
        )


class AlmaFluxService:
    def __init__(self, base_url: str = "https://almascience.org/sc/flux", timeout_sec: float = 30.0) -> None:
        self._base_url = base_url
        self._timeout_sec = timeout_sec

    def get_estimate(self, source: str, obs_date_utc: datetime, freq_ghz: float) -> FluxEstimate:
        query_date = obs_date_utc.strftime("%d-%b-%Y")
        names_to_try = [source]
        if not source.upper().startswith("J"):
            names_to_try = [f"J{source}", source]

        failures: List[str] = []
        for name in names_to_try:
            print(
                f"Trying ALMA VO service with source name {name} at {freq_ghz:g} GHz for date {query_date}",
                file=sys.stderr,
            )
            status_code, date_used, flux_jy, nearest_obs_days = self._query(
                name=name, query_date=query_date, freq_ghz=freq_ghz
            )
            if status_code in (0, 1) and flux_jy is not None:
                print(
                    f"ALMA VO lookup succeeded for source name {name} (status={status_code})",
                    file=sys.stderr,
                )
                return FluxEstimate(
                    flux_mjy=flux_jy * 1000.0,
                    estimate_obs_date_utc=date_used,
                    backend="alma_vo",
                    alma_nearest_obs_days=nearest_obs_days,
                )
            failures.append(f"name={name}: status={status_code}, flux_jy={flux_jy}")

        details = "; ".join(failures) if failures else "no response rows"
        raise LookupError(f"ALMA VO returned no valid flux estimate for source={source!r}, freq={freq_ghz:g} GHz ({details}).")

    def _query(self, name: str, query_date: str, freq_ghz: float) -> tuple[int, str, Optional[float], Optional[float]]:
        query = urlencode(
            {
                "DATE": query_date,
                "FREQUENCY": f"{freq_ghz * 1.0e9:.1E}",
                "NAME": name,
            }
        )
        url = f"{self._base_url}?{query}"
        with urlopen(url, timeout=self._timeout_sec) as response:
            root = ET.fromstring(response.read())

        cells = [cell.text.strip() if cell.text else "" for cell in root.findall(".//TD")]
        if len(cells) < 6:
            raise LookupError(f"Unexpected ALMA VO response format for NAME={name!r}, FREQ={freq_ghz:g} GHz.")

        status_text = cells[0]
        date_text = cells[3]
        flux_text = cells[4]
        nearest_obs_text = cells[9] if len(cells) > 9 else ""
        try:
            status_code = int(status_text)
        except ValueError as exc:
            raise LookupError(
                f"ALMA VO returned non-integer status {status_text!r} for NAME={name!r}, FREQ={freq_ghz:g} GHz."
            ) from exc

        nearest_obs_days: Optional[float] = None
        if nearest_obs_text and nearest_obs_text.lower() != "null":
            try:
                nearest_obs_days = float(nearest_obs_text)
            except ValueError:
                nearest_obs_days = None

        if not flux_text or flux_text.lower() == "null":
            return status_code, date_text, None, nearest_obs_days
        return status_code, date_text, float(flux_text), nearest_obs_days


class FallbackFluxEstimator:
    def __init__(self, primary: FluxEstimator, fallback: FluxEstimator) -> None:
        self._primary = primary
        self._fallback = fallback

    def get_estimate(self, source: str, obs_date_utc: datetime, freq_ghz: float) -> FluxEstimate:
        try:
            return self._primary.get_estimate(source=source, obs_date_utc=obs_date_utc, freq_ghz=freq_ghz)
        except Exception:
            return self._fallback.get_estimate(source=source, obs_date_utc=obs_date_utc, freq_ghz=freq_ghz)


class SmaCatalogFluxService:
    CATALOG_URL = "http://sma1.sma.hawaii.edu/callist/callist.html"
    NU_1MM_GHZ = 299.792458

    def __init__(self, timeout_sec: float = 30.0, spectral_index: float = -0.7) -> None:
        self._timeout_sec = timeout_sec
        self._spectral_index = spectral_index
        self._cache: Optional[Dict[str, Dict[str, tuple[datetime, float]]]] = None

    def get_estimate(self, source: str, obs_date_utc: datetime, freq_ghz: float) -> FluxEstimate:
        catalog = self._load_catalog()
        candidates = self._source_candidates(source)

        matched_source: Optional[str] = None
        entry: Optional[Dict[str, tuple[datetime, float]]] = None
        for candidate in candidates:
            if candidate in catalog:
                matched_source = candidate
                entry = catalog[candidate]
                break

        if entry is None:
            raise LookupError(f"SMA catalog has no entry matching source={source!r} (candidates={candidates}).")

        if "1mm" not in entry:
            raise LookupError(f"SMA catalog source={matched_source!r} does not have a 1mm flux entry.")

        one_mm_date, one_mm_flux_jy = entry["1mm"]
        if one_mm_flux_jy <= 0.0:
            raise LookupError(f"SMA catalog has non-positive flux values for source={matched_source!r}.")

        estimate_jy = one_mm_flux_jy * (freq_ghz / self.NU_1MM_GHZ) ** self._spectral_index
        delta_days = (obs_date_utc - one_mm_date).total_seconds() / 86400.0

        print(
            f"SMA fallback used for source {matched_source} at {freq_ghz:g} GHz "
            f"(1mm_date={one_mm_date.date()}, alpha={self._spectral_index:.3f}, delta_days={delta_days:.2f})",
            file=sys.stderr,
        )
        return FluxEstimate(
            flux_mjy=estimate_jy * 1000.0,
            estimate_obs_date_utc=one_mm_date.strftime("%Y-%m-%d"),
            backend="sma_catalog",
            sma_1mm_date_delta_days=delta_days,
        )

    def _load_catalog(self) -> Dict[str, Dict[str, tuple[datetime, float]]]:
        if self._cache is not None:
            return self._cache

        with urlopen(self.CATALOG_URL, timeout=self._timeout_sec) as response:
            html_text = response.read().decode("iso-8859-1", errors="replace")

        table_match = re.search(r'<table align="center" border="1".*?</table>', html_text, flags=re.S | re.I)
        if not table_match:
            raise LookupError("Could not find SMA calibrator table in catalog page.")

        table_html = table_match.group(0)
        row_html_list = re.findall(r"<tr[^>]*>(.*?)</tr>", table_html, flags=re.S | re.I)
        catalog: Dict[str, Dict[str, tuple[datetime, float]]] = {}

        current_source: Optional[str] = None
        for row_html in row_html_list:
            raw_cells = re.findall(r"<td[^>]*>(.*?)</td>", row_html, flags=re.S | re.I)
            if not raw_cells:
                continue
            cells = [self._clean_cell(cell) for cell in raw_cells]

            band: Optional[str] = None
            date_text: Optional[str] = None
            flux_text: Optional[str] = None

            if len(cells) >= 8:
                # Full row with source columns and a band entry.
                current_source = cells[1]
                band = cells[4]
                date_text = cells[5]
                flux_text = cells[7]
            elif len(cells) >= 4 and current_source:
                # Continuation row (typically the second band for same source).
                band = cells[0]
                date_text = cells[1]
                flux_text = cells[3]

            if not current_source or not band or not date_text or not flux_text:
                continue

            band_key = self._normalize_band(band)
            if band_key != "1mm":
                continue

            date_val = self._parse_sma_date(date_text)
            flux_val = self._parse_flux_jy(flux_text)
            key = self._normalize_source(current_source)
            catalog.setdefault(key, {})[band_key] = (date_val, flux_val)

        self._cache = catalog
        return catalog

    @staticmethod
    def _clean_cell(text: str) -> str:
        no_tags = re.sub(r"<[^>]*>", "", text)
        return html.unescape(no_tags).replace("\xa0", " ").strip()

    @staticmethod
    def _normalize_source(source: str) -> str:
        return source.strip().upper().replace(" ", "")

    def _source_candidates(self, source: str) -> List[str]:
        s = self._normalize_source(source)
        candidates = [s]
        if s.startswith("J"):
            candidates.append(s[1:])
        else:
            candidates.append(f"J{s}")
        # Preserve order while deduping.
        seen = set()
        deduped: List[str] = []
        for cand in candidates:
            if cand not in seen:
                deduped.append(cand)
                seen.add(cand)
        return deduped

    @staticmethod
    def _normalize_band(text: str) -> str:
        t = text.strip().lower()
        if t in ("1mm",):
            return "1mm"
        if "850" in t:
            return "850um"
        return t

    @staticmethod
    def _parse_sma_date(text: str) -> datetime:
        return datetime.strptime(text.strip(), "%d %b %Y")

    @staticmethod
    def _parse_flux_jy(text: str) -> float:
        # e.g., "0.87 ± 0.03"
        match = re.search(r"[-+]?[0-9]*\.?[0-9]+", text)
        if not match:
            raise ValueError(f"Could not parse SMA flux from {text!r}")
        return float(match.group(0))


def parse_obs_datetime(raw: object) -> datetime:
    if raw is None:
        raise ValueError("ECSV metadata is missing observation date.")

    text = str(raw).strip()
    fmts = [
        "%Y-%m-%d.%H:%M:%S",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%dT%H:%M:%S",
        "%Y-%m-%d",
    ]
    for fmt in fmts:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            pass

    # Fall back to fromisoformat for fractional seconds/timezone offsets.
    try:
        return datetime.fromisoformat(text)
    except ValueError as exc:
        raise ValueError(f"Unsupported observation date format: {text!r}") from exc


def parse_array_map(spec: str) -> Dict[int, float]:
    mapping: Dict[int, float] = {}
    for item in spec.split(","):
        item = item.strip()
        if not item:
            continue
        if ":" not in item:
            raise ValueError(f"Invalid array map entry {item!r}. Expected format 'array:frequency'.")
        array_text, freq_text = item.split(":", maxsplit=1)
        mapping[int(array_text.strip())] = float(freq_text.strip())
    if not mapping:
        raise ValueError("Array map cannot be empty.")
    return mapping


def format_array_map(mapping: Dict[int, float]) -> str:
    items = sorted(mapping.items(), key=lambda kv: kv[0])
    return ",".join(f"{array}:{freq:g}" for array, freq in items)


def load_measurements(
    ecsv_path: Path, array_to_freq: Dict[int, float]
) -> tuple[str, datetime, Optional[str], Dict[float, Dict[str, float]]]:
    table = Table.read(str(ecsv_path), format="ascii.ecsv")

    source = table.meta.get("source")
    if not source:
        raise ValueError(f"ECSV metadata in {ecsv_path} is missing 'source'.")
    obsnum = table.meta.get("obsnum")

    obs_date_raw = table.meta.get("date") or table.meta.get("creation_date")
    obs_date = parse_obs_datetime(obs_date_raw)

    for col in ("array", "amp"):
        if col not in table.colnames:
            raise ValueError(f"Required column {col!r} not found in {ecsv_path}.")

    has_amp_err = "amp_err" in table.colnames

    measured_by_freq: Dict[float, Dict[str, float]] = {}
    for row in table:
        array_index = int(row["array"])
        if array_index not in array_to_freq:
            continue
        freq_ghz = float(array_to_freq[array_index])
        measured_by_freq[freq_ghz] = {
            "array_index": array_index,
            "amp": float(row["amp"]),
            "amp_err": float(row["amp_err"]) if has_amp_err else None,
        }

    return str(source), obs_date, (str(obsnum) if obsnum is not None else None), measured_by_freq


def resolve_output_path_with_obsnum(output_path: Path, obsnum: Optional[str]) -> Path:
    output_path = Path(output_path)
    if output_path.name != "ratios.json":
        return output_path
    if not obsnum:
        return output_path
    return output_path.with_name(f"ratios_obsnum_{obsnum}.json")


def compute_calibration(
    source: str,
    obs_date: datetime,
    measured_by_freq: Dict[float, Dict[str, float]],
    estimator: FluxEstimator,
    frequencies_ghz: Iterable[float],
) -> List[CalibrationResult]:
    results: List[CalibrationResult] = []

    for freq in frequencies_ghz:
        if freq not in measured_by_freq:
            raise ValueError(f"No measured flux in ECSV for {freq:g} GHz.")

        measured = measured_by_freq[freq]
        measured_flux_mjy = float(measured["amp"])
        estimate = estimator.get_estimate(source=source, obs_date_utc=obs_date, freq_ghz=freq)
        estimated_flux_mjy = estimate.flux_mjy
        estimate_obs_date = estimate.estimate_obs_date_utc
        if estimated_flux_mjy == 0.0:
            raise ZeroDivisionError(f"Estimated flux is zero at {freq:g} GHz; cannot compute ratio.")

        results.append(
            CalibrationResult(
                source=source,
                obs_date_utc=obs_date.isoformat(sep=" "),
                frequency_ghz=float(freq),
                array_index=int(measured["array_index"]),
                measured_flux_mjy=measured_flux_mjy,
                measured_flux_err_mjy=measured["amp_err"],
                estimated_flux_mjy=estimated_flux_mjy,
                estimate_obs_date_utc=estimate_obs_date,
                calibration_ratio=measured_flux_mjy / estimated_flux_mjy,
                backend_used=estimate.backend,
                sma_1mm_date_delta_days=estimate.sma_1mm_date_delta_days,
                alma_nearest_obs_days=estimate.alma_nearest_obs_days,
            )
        )

    return results


def print_results(results: List[CalibrationResult]) -> None:
    header = (
        "frequency_ghz",
        "array",
        "measured_mJy",
        "reference_mJy",
        "ratio(meas/ref)",
        "backend_used",
        "alma_nearest_obs_days",
        "sma_1mm_date_delta_days",
        "estimate_obs_date",
    )
    print("\t".join(header))
    for r in results:
        print(
            f"{r.frequency_ghz:.1f}\t{r.array_index}\t{r.measured_flux_mjy:.6g}\t"
            f"{r.estimated_flux_mjy:.6g}\t{r.calibration_ratio:.6g}\t{r.backend_used}\t"
            f"{'' if r.alma_nearest_obs_days is None else f'{r.alma_nearest_obs_days:.3f}'}\t"
            f"{'' if r.sma_1mm_date_delta_days is None else f'{r.sma_1mm_date_delta_days:.2f}'}\t{r.estimate_obs_date_utc}"
        )


def write_results(results: List[CalibrationResult], output_path: Path) -> None:
    output_path = Path(output_path)
    rows = [asdict(r) for r in results]

    if output_path.suffix.lower() == ".json":
        output_path.write_text(json.dumps(rows, indent=2), encoding="utf-8")
        return

    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare measured ECSV fluxes to ALMA VO estimates with SMA(1mm+alpha)/SQLite fallback at 150/220/280 GHz."
    )
    parser.add_argument("ecsv", type=Path, help="Path to pointing ECSV file")
    parser.add_argument("db", type=Path, nargs="?", help="Optional path to SQLite flux-estimate database")
    parser.add_argument(
        "--backend",
        choices=("alma", "sqlite"),
        default="alma",
        help="Primary flux estimate backend. Default: alma (ALMA VO service).",
    )
    parser.add_argument(
        "--table",
        default="flux_estimates",
        help="Database table name containing source/date/frequency/flux columns (default: flux_estimates)",
    )
    parser.add_argument(
        "--array-map",
        default=format_array_map(DEFAULT_ARRAY_MAP),
        help="Mapping from ECSV array index to frequency GHz, e.g. '0:150,1:220,2:280'",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output file (.csv or .json). If omitted, prints table to stdout only.",
    )
    parser.add_argument(
        "--sma-spectral-index",
        type=float,
        default=-0.7,
        help="Spectral index alpha for SMA 1mm extrapolation (S_nu ~ nu^alpha). Default: -0.7.",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    array_map = parse_array_map(args.array_map)
    source, obs_date, obsnum, measured_by_freq = load_measurements(args.ecsv, array_map)

    sqlite_db = SQLiteFluxDatabase(args.db, table_name=args.table) if args.db else None
    if args.backend == "sqlite" and sqlite_db is None:
        raise ValueError("SQLite backend selected but no database path was provided.")

    if args.backend == "sqlite":
        estimator: FluxEstimator = sqlite_db  # type: ignore[assignment]
    else:
        alma_estimator = AlmaFluxService()
        sma_estimator = SmaCatalogFluxService(spectral_index=args.sma_spectral_index)
        estimator = FallbackFluxEstimator(primary=alma_estimator, fallback=sma_estimator)
        if sqlite_db is not None:
            estimator = FallbackFluxEstimator(primary=estimator, fallback=sqlite_db)

    try:
        results = compute_calibration(
            source=source,
            obs_date=obs_date,
            measured_by_freq=measured_by_freq,
            estimator=estimator,
            frequencies_ghz=TARGET_FREQUENCIES_GHZ,
        )
    finally:
        if sqlite_db is not None:
            sqlite_db.close()

    print_results(results)

    if args.output:
        final_output_path = resolve_output_path_with_obsnum(args.output, obsnum)
        if final_output_path != args.output:
            print(f"Writing output to {final_output_path} (derived from obsnum={obsnum})", file=sys.stderr)
        write_results(results, final_output_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
