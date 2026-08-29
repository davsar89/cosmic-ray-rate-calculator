from __future__ import annotations

import argparse
import bisect
import csv
import hashlib
import math
import re
import site
import statistics
import sys
import sysconfig
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Literal, Sequence


CSV_HEADER = ("energy_eV", "flux_per_m2_sr_s_GeV")
MAX_SAMPLE_POINTS = 1_000_000

# The CSV schema accepts plain ASCII decimal or scientific notation only.
# Bare float() is laxer (underscores, Unicode digits, inf/nan spellings), so
# cells are checked against this pattern first.
_NUMERIC_CELL_RE = re.compile(
    r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\Z", re.ASCII
)
TailModel = Literal["truncate", "last-segment"]
Geometry = Literal["flat", "hemisphere", "full_sky"]

# Break points for the documented piecewise approximation. The rates and
# exponents are derived from the selected spectrum rather than duplicated as
# hard-coded constants.
APPROX_BREAKS_EV = (
    1.0e10,
    1.0e12,
    2.0e15,
    5.0e16,
    3.0e18,
    3.0e19,
    1.0e20,
)


def _find_default_csv() -> Path:
    """Find data.csv in a source checkout or in known install data locations.

    Wheel data-files land in the active install scheme's data root, which is
    sys.prefix only for venv-style installs; user and system schemes differ,
    so several candidates are checked.
    """
    source_tree_path = Path(__file__).with_name("data.csv")
    share_suffix = Path("share") / "cosmic-ray-rate-calculator" / "data.csv"

    candidates = [source_tree_path]
    data_root = sysconfig.get_path("data")
    if data_root:
        candidates.append(Path(data_root) / share_suffix)
    try:
        user_base = site.getuserbase()
    except (AttributeError, OSError):
        user_base = None
    if user_base:
        candidates.append(Path(user_base) / share_suffix)
    candidates.append(Path(sys.prefix) / share_suffix)

    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return source_tree_path


DEFAULT_CSV = _find_default_csv()


@dataclass(frozen=True)
class Spectrum:
    """Validated log-log spectrum plus precomputed segment integrals."""

    energies_eV: tuple[float, ...]
    fluxes_per_m2_sr_s_GeV: tuple[float, ...]
    segment_gammas: tuple[float, ...]
    segment_integrals_eV: tuple[float, ...]
    cumulative_integrals_eV: tuple[float, ...]
    source_path: Path
    source_sha256: str
    dataset_sha256: str

    @property
    def min_energy_eV(self) -> float:
        return self.energies_eV[0]

    @property
    def max_energy_eV(self) -> float:
        return self.energies_eV[-1]


@dataclass(frozen=True)
class Approximation:
    """Piecewise power-law approximation calibrated to one spectrum/tail model."""

    breaks_eV: tuple[float, ...]
    rates_per_s_per_km2: tuple[float, ...]
    exponents: tuple[float, ...]
    dataset_sha256: str
    tail_model: TailModel


# ---------------------------------------------------------------------------
# Validation and data loading
# ---------------------------------------------------------------------------


def _require_positive_finite(value: float, name: str) -> float:
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a finite number greater than zero")
    return value


def _normalise_tail_model(tail_model: str) -> TailModel:
    if tail_model not in {"truncate", "last-segment"}:
        raise ValueError("tail_model must be one of: 'truncate', 'last-segment'")
    return tail_model  # type: ignore[return-value]


def _normalise_geometry(geometry: str) -> Geometry:
    if geometry not in {"flat", "hemisphere", "full_sky"}:
        raise ValueError("geometry must be one of: 'flat', 'hemisphere', 'full_sky'")
    return geometry  # type: ignore[return-value]


def _integrate_power_law_from_reference(
    e_lo_eV: float,
    e_hi_eV: float,
    f_lo: float,
    gamma: float,
) -> float:
    """Integrate f_lo * (E / e_lo_eV)**gamma over [e_lo_eV, e_hi_eV]."""
    if not (0.0 < e_lo_eV < e_hi_eV):
        raise ValueError("power-law integration bounds must satisfy 0 < lower < upper")
    _require_positive_finite(f_lo, "f_lo")
    if not math.isfinite(gamma):
        raise ValueError("power-law exponent must be finite")

    log_ratio = math.log(e_hi_eV / e_lo_eV)
    exponent_plus_one = gamma + 1.0

    try:
        if math.isclose(exponent_plus_one, 0.0, rel_tol=0.0, abs_tol=1.0e-12):
            integral = f_lo * e_lo_eV * log_ratio
        else:
            # expm1 is stable when gamma is close to -1 and avoids constructing a
            # coefficient involving very large powers of absolute energy.
            integral = (
                f_lo
                * e_lo_eV
                * math.expm1(exponent_plus_one * log_ratio)
                / exponent_plus_one
            )
    except OverflowError as exc:
        raise ValueError("power-law segment integral overflowed") from exc

    if not math.isfinite(integral) or integral <= 0.0:
        raise ValueError("power-law segment integral is not finite and positive")
    return integral


def _integrate_power_law_segment(
    e_lo_eV: float,
    e_hi_eV: float,
    f_lo: float,
    f_hi: float,
) -> float:
    """Integrate a log-log segment defined by two positive endpoints."""
    _require_positive_finite(e_lo_eV, "e_lo_eV")
    _require_positive_finite(e_hi_eV, "e_hi_eV")
    if e_lo_eV >= e_hi_eV:
        raise ValueError("segment energies must satisfy e_lo_eV < e_hi_eV")
    _require_positive_finite(f_lo, "f_lo")
    _require_positive_finite(f_hi, "f_hi")
    gamma = math.log(f_hi / f_lo) / math.log(e_hi_eV / e_lo_eV)
    return _integrate_power_law_from_reference(e_lo_eV, e_hi_eV, f_lo, gamma)


def _parse_spectrum(path: Path, text: str, digest: str) -> Spectrum:
    meaningful_lines = [
        (line_number, raw_line)
        for line_number, raw_line in enumerate(text.splitlines(), start=1)
        if raw_line.strip() and not raw_line.lstrip().startswith("#")
    ]
    if not meaningful_lines:
        raise ValueError(f"{path}: no CSV header or data rows found")

    header_line_number, header_line = meaningful_lines[0]
    try:
        header = tuple(cell.strip() for cell in next(csv.reader([header_line])))
    except csv.Error as exc:
        raise ValueError(f"{path}:{header_line_number}: invalid CSV header: {exc}") from exc

    if header != CSV_HEADER:
        expected = ",".join(CSV_HEADER)
        actual = ",".join(header)
        raise ValueError(
            f"{path}:{header_line_number}: expected CSV header '{expected}', got '{actual}'"
        )

    records: list[tuple[float, float, int]] = []
    for line_number, raw_line in meaningful_lines[1:]:
        try:
            row = next(csv.reader([raw_line]))
        except csv.Error as exc:
            raise ValueError(f"{path}:{line_number}: invalid CSV row: {exc}") from exc

        if len(row) != 2:
            raise ValueError(
                f"{path}:{line_number}: expected 2 columns, found {len(row)}"
            )

        energy_text = row[0].strip()
        flux_text = row[1].strip()
        if not (_NUMERIC_CELL_RE.match(energy_text) and _NUMERIC_CELL_RE.match(flux_text)):
            raise ValueError(
                f"{path}:{line_number}: energy and flux must be numeric "
                "(plain decimal or scientific notation)"
            )
        energy_eV = float(energy_text)
        flux = float(flux_text)

        _require_positive_finite(energy_eV, f"{path}:{line_number}: energy_eV")
        _require_positive_finite(
            flux, f"{path}:{line_number}: flux_per_m2_sr_s_GeV"
        )
        records.append((energy_eV, flux, line_number))

    if len(records) < 2:
        raise ValueError(f"{path}: at least two numeric spectrum rows are required")

    # Row order is not part of the schema. Sorting avoids the old failure mode
    # where reversed data silently caused energy and flux columns to be swapped.
    records.sort(key=lambda record: record[0])
    for previous, current in zip(records, records[1:]):
        if previous[0] == current[0]:
            raise ValueError(
                f"{path}: duplicate energy_eV value {current[0]:.12g} "
                f"on input lines {previous[2]} and {current[2]}"
            )

    energies = tuple(record[0] for record in records)
    fluxes = tuple(record[1] for record in records)

    gammas: list[float] = []
    segment_integrals: list[float] = []
    for e_lo, e_hi, f_lo, f_hi in zip(energies, energies[1:], fluxes, fluxes[1:]):
        gamma = math.log(f_hi / f_lo) / math.log(e_hi / e_lo)
        if not math.isfinite(gamma):
            raise ValueError(f"{path}: non-finite log-log slope in spectrum")
        gammas.append(gamma)
        segment_integrals.append(
            _integrate_power_law_from_reference(e_lo, e_hi, f_lo, gamma)
        )

    cumulative = [0.0] * len(energies)
    for index in range(len(segment_integrals) - 1, -1, -1):
        cumulative[index] = cumulative[index + 1] + segment_integrals[index]

    semantic_hasher = hashlib.sha256()
    for energy_eV, flux in zip(energies, fluxes):
        semantic_hasher.update(energy_eV.hex().encode("ascii"))
        semantic_hasher.update(b",")
        semantic_hasher.update(flux.hex().encode("ascii"))
        semantic_hasher.update(b"\n")

    return Spectrum(
        energies_eV=energies,
        fluxes_per_m2_sr_s_GeV=fluxes,
        segment_gammas=tuple(gammas),
        segment_integrals_eV=tuple(segment_integrals),
        cumulative_integrals_eV=tuple(cumulative),
        source_path=path,
        source_sha256=digest,
        dataset_sha256=semantic_hasher.hexdigest(),
    )


# Parsed spectra keyed on (resolved path, content SHA-256). Keying on the
# content digest instead of stat metadata means a rewritten file is always
# re-parsed, even when its size and mtime survive the rewrite.
_SPECTRUM_CACHE_MAX_ENTRIES = 32
_spectrum_cache: dict[tuple[str, str], Spectrum] = {}


def load_spectrum(csv_path: Path | str = DEFAULT_CSV) -> Spectrum:
    """Load, validate, sort, and precompute a spectrum from the defined CSV schema."""
    path = Path(csv_path).expanduser().resolve()
    if path.is_dir():
        raise ValueError(f"Spectrum path is not a regular file: {path}")
    try:
        raw_bytes = path.read_bytes()
    except OSError as exc:
        raise ValueError(f"Could not read spectrum file {path}: {exc.strerror or exc}") from exc
    digest = hashlib.sha256(raw_bytes).hexdigest()

    key = (str(path), digest)
    cached = _spectrum_cache.get(key)
    if cached is not None:
        return cached

    try:
        # utf-8-sig drops a leading BOM before comment/header filtering sees it.
        text = raw_bytes.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise ValueError(f"{path}: CSV file must be UTF-8 text") from exc
    spectrum = _parse_spectrum(path, text, digest)
    if len(_spectrum_cache) >= _SPECTRUM_CACHE_MAX_ENTRIES:
        _spectrum_cache.clear()
    _spectrum_cache[key] = spectrum
    return spectrum


# ---------------------------------------------------------------------------
# Spectrum integration and ideal acceptance models
# ---------------------------------------------------------------------------


def _finite_integral_to_csv_max_eV(threshold_eV: float, spectrum: Spectrum) -> float:
    if threshold_eV == spectrum.max_energy_eV:
        return 0.0

    segment_index = bisect.bisect_right(spectrum.energies_eV, threshold_eV) - 1
    e_lo = spectrum.energies_eV[segment_index]
    e_hi = spectrum.energies_eV[segment_index + 1]
    f_lo = spectrum.fluxes_per_m2_sr_s_GeV[segment_index]
    gamma = spectrum.segment_gammas[segment_index]

    try:
        f_at_threshold = f_lo * (threshold_eV / e_lo) ** gamma
    except OverflowError as exc:
        raise ValueError("interpolated flux overflowed") from exc
    if not math.isfinite(f_at_threshold) or f_at_threshold <= 0.0:
        raise ValueError("interpolated flux is not finite and positive")
    partial_segment = _integrate_power_law_from_reference(
        threshold_eV, e_hi, f_at_threshold, gamma
    )
    return partial_segment + spectrum.cumulative_integrals_eV[segment_index + 1]


def _last_segment_tail_integral_eV(threshold_eV: float, spectrum: Spectrum) -> float:
    gamma = spectrum.segment_gammas[-1]
    if gamma >= -1.0:
        raise ValueError(
            "The final log-log segment has slope gamma >= -1, so its integral "
            "to infinity does not converge"
        )

    e_max = spectrum.max_energy_eV
    f_max = spectrum.fluxes_per_m2_sr_s_GeV[-1]
    if threshold_eV <= e_max:
        e_start = e_max
        f_start = f_max
    else:
        e_start = threshold_eV
        try:
            f_start = f_max * (threshold_eV / e_max) ** gamma
        except OverflowError as exc:
            raise ValueError("extrapolated tail flux overflowed") from exc

    tail = f_start * e_start / (-(gamma + 1.0))
    if not math.isfinite(tail) or tail <= 0.0:
        raise ValueError("Last-segment tail integral is not finite and positive")
    return tail


def _integrated_intensity_above_from_spectrum(
    threshold_eV: float,
    spectrum: Spectrum,
    tail_model: TailModel,
) -> float:
    _require_positive_finite(threshold_eV, "threshold_eV")

    if threshold_eV < spectrum.min_energy_eV:
        raise ValueError(
            f"threshold_eV must be at least {spectrum.min_energy_eV:.6g} eV; "
            "low-energy extrapolation is not implemented"
        )

    if tail_model == "truncate":
        if threshold_eV > spectrum.max_energy_eV:
            raise ValueError(
                f"threshold_eV exceeds the CSV maximum of "
                f"{spectrum.max_energy_eV:.6g} eV; choose "
                "tail_model='last-segment' only if that extrapolation is intended"
            )
        integral_eV = _finite_integral_to_csv_max_eV(threshold_eV, spectrum)
    else:
        if threshold_eV < spectrum.max_energy_eV:
            integral_eV = _finite_integral_to_csv_max_eV(threshold_eV, spectrum)
            integral_eV += _last_segment_tail_integral_eV(
                spectrum.max_energy_eV, spectrum
            )
        else:
            integral_eV = _last_segment_tail_integral_eV(threshold_eV, spectrum)

    # The CSV flux is differential per GeV while the energy coordinates are eV:
    # dE_GeV = dE_eV / 1e9.
    return integral_eV * 1.0e-9


def integrated_intensity_above(
    threshold_eV: float,
    csv_path: Path | str = DEFAULT_CSV,
    *,
    tail_model: TailModel = "truncate",
) -> float:
    """Return the threshold-integrated intensity in m^-2 sr^-1 s^-1.

    ``tail_model='truncate'`` integrates only through the largest tabulated
    energy. It does not claim that the physical flux is zero beyond that point.

    ``tail_model='last-segment'`` explicitly extends the final CSV log-log slope
    to infinity. This is a model-dependent extrapolation and can dominate the
    result near the upper end of the data.
    """
    normalised_tail_model = _normalise_tail_model(tail_model)
    spectrum = load_spectrum(csv_path)
    return _integrated_intensity_above_from_spectrum(
        threshold_eV, spectrum, normalised_tail_model
    )


def _geometry_factor_sr(geometry: Geometry) -> float:
    if geometry == "flat":
        # Integral of cos(theta) dOmega over a downward hemisphere.
        return math.pi
    if geometry == "hemisphere":
        # Direction-independent effective area over 2*pi sr. This is not the
        # projected acceptance of an ordinary flat detector.
        return 2.0 * math.pi
    return 4.0 * math.pi


def _rate_above_from_spectrum(
    threshold_eV: float,
    spectrum: Spectrum,
    *,
    area_km2: float,
    geometry: Geometry,
    tail_model: TailModel,
) -> float:
    _require_positive_finite(area_km2, "area_km2")
    intensity = _integrated_intensity_above_from_spectrum(
        threshold_eV, spectrum, tail_model
    )
    area_m2 = area_km2 * 1.0e6
    rate = intensity * _geometry_factor_sr(geometry) * area_m2
    if not math.isfinite(rate):
        raise ValueError("calculated rate is not finite; reduce the supplied area")
    return rate


def rate_above(
    threshold_eV: float,
    area_km2: float = 1.0,
    geometry: Geometry = "flat",
    csv_path: Path | str = DEFAULT_CSV,
    *,
    tail_model: TailModel = "truncate",
) -> float:
    """Return an ideal crossing-rate estimate in events/s.

    This function applies only a geometric/angular factor to an isotropic
    intensity. It does not include detector efficiency, dead time, trigger or
    reconstruction response, atmospheric shower response, or an energy- and
    direction-dependent effective area.

    ``geometry='flat'`` is a one-sided horizontal surface exposed to a downward
    isotropic intensity and uses the projected factor pi*A.

    ``geometry='hemisphere'`` and ``geometry='full_sky'`` treat ``area_km2`` as a
    direction-independent effective area integrated over 2*pi or 4*pi sr. They
    are aperture conventions, not ordinary flat-detector geometries.
    """
    normalised_geometry = _normalise_geometry(geometry)
    normalised_tail_model = _normalise_tail_model(tail_model)
    spectrum = load_spectrum(csv_path)
    return _rate_above_from_spectrum(
        threshold_eV,
        spectrum,
        area_km2=area_km2,
        geometry=normalised_geometry,
        tail_model=normalised_tail_model,
    )


def rate_above_from_aperture(
    threshold_eV: float,
    aperture_m2_sr: float,
    csv_path: Path | str = DEFAULT_CSV,
    *,
    tail_model: TailModel = "truncate",
) -> float:
    """Return events/s for a caller-supplied, energy-independent aperture.

    ``aperture_m2_sr`` is the integrated effective area in m^2 sr. For realistic
    instruments it generally depends on energy and direction; callers should
    integrate that response separately rather than using this scalar shortcut.
    """
    _require_positive_finite(aperture_m2_sr, "aperture_m2_sr")
    intensity = integrated_intensity_above(
        threshold_eV, csv_path=csv_path, tail_model=tail_model
    )
    rate = intensity * aperture_m2_sr
    if not math.isfinite(rate):
        raise ValueError("calculated rate is not finite; reduce the supplied aperture")
    return rate


# ---------------------------------------------------------------------------
# Piecewise approximation
# ---------------------------------------------------------------------------


def _validate_breaks(
    breaks_eV: Sequence[float], spectrum: Spectrum, tail_model: TailModel
) -> tuple[float, ...]:
    breaks = tuple(float(value) for value in breaks_eV)
    if len(breaks) < 2:
        raise ValueError("At least two approximation break energies are required")
    for index, value in enumerate(breaks):
        _require_positive_finite(value, f"breaks_eV[{index}]")
    if any(left >= right for left, right in zip(breaks, breaks[1:])):
        raise ValueError("Approximation break energies must be strictly increasing")
    if breaks[0] < spectrum.min_energy_eV:
        raise ValueError(
            f"First approximation break must be at least "
            f"{spectrum.min_energy_eV:.6g} eV"
        )
    if tail_model == "truncate" and breaks[-1] >= spectrum.max_energy_eV:
        raise ValueError(
            "For a truncated spectrum, the final approximation break must be "
            "strictly below the CSV maximum so the reference rate remains positive"
        )
    return breaks


def _derive_default_breaks(spectrum: Spectrum, tail_model: TailModel) -> tuple[float, ...]:
    """Restrict APPROX_BREAKS_EV to break energies the dataset can support."""
    breaks = tuple(
        value
        for value in APPROX_BREAKS_EV
        if value >= spectrum.min_energy_eV
        and (tail_model != "truncate" or value < spectrum.max_energy_eV)
    )
    if len(breaks) < 2:
        raise ValueError(
            "The dataset energy range is too narrow for the built-in "
            "approximation break points; supply explicit break energies "
            "with --approx-breaks-ev"
        )
    return breaks


@lru_cache(maxsize=32)
def _calibrate_approximation_cached(
    spectrum: Spectrum,
    breaks_eV: tuple[float, ...],
    tail_model: TailModel,
) -> Approximation:
    rates = tuple(
        _rate_above_from_spectrum(
            threshold,
            spectrum,
            area_km2=1.0,
            geometry="flat",
            tail_model=tail_model,
        )
        for threshold in breaks_eV
    )
    if any(not math.isfinite(rate) or rate <= 0.0 for rate in rates):
        raise ValueError("Approximation calibration produced a non-positive rate")

    exponents = tuple(
        -math.log(rate_hi / rate_lo) / math.log(e_hi / e_lo)
        for e_lo, e_hi, rate_lo, rate_hi in zip(
            breaks_eV, breaks_eV[1:], rates, rates[1:]
        )
    )
    return Approximation(
        breaks_eV=breaks_eV,
        rates_per_s_per_km2=rates,
        exponents=exponents,
        dataset_sha256=spectrum.dataset_sha256,
        tail_model=tail_model,
    )


def calibrate_piecewise_approximation(
    csv_path: Path | str = DEFAULT_CSV,
    *,
    breaks_eV: Sequence[float] = APPROX_BREAKS_EV,
    tail_model: TailModel = "truncate",
) -> Approximation:
    """Explicitly calibrate the piecewise approximation to a selected dataset."""
    normalised_tail_model = _normalise_tail_model(tail_model)
    spectrum = load_spectrum(csv_path)
    validated_breaks = _validate_breaks(breaks_eV, spectrum, normalised_tail_model)
    return _calibrate_approximation_cached(
        spectrum, validated_breaks, normalised_tail_model
    )


def _approximation_matches(
    approximation: Approximation, spectrum: Spectrum, tail_model: TailModel
) -> bool:
    return (
        approximation.dataset_sha256 == spectrum.dataset_sha256
        and approximation.tail_model == tail_model
    )


@lru_cache(maxsize=1)
def _default_approximation() -> Approximation:
    # The bundled dataset is treated as immutable within a process; paths that
    # depend on the dataset re-check _approximation_matches against a freshly
    # loaded spectrum, so a mid-process edit raises instead of mismatching.
    return calibrate_piecewise_approximation(DEFAULT_CSV)


def approx_rate_above_flat_km2(
    threshold_eV: float,
    area_km2: float = 1.0,
    *,
    approximation: Approximation | None = None,
) -> float:
    """Evaluate a calibrated piecewise power-law rate approximation.

    With no ``approximation`` argument, the function uses the bundled data.csv,
    the explicit CSV-truncation convention, and ``APPROX_BREAKS_EV``. A custom
    dataset must first be passed to :func:`calibrate_piecewise_approximation`.
    """
    _require_positive_finite(threshold_eV, "threshold_eV")
    _require_positive_finite(area_km2, "area_km2")
    selected = approximation if approximation is not None else _default_approximation()

    if (
        len(selected.breaks_eV) < 2
        or len(selected.rates_per_s_per_km2) != len(selected.breaks_eV)
        or len(selected.exponents) != len(selected.breaks_eV) - 1
    ):
        raise ValueError(
            "approximation is structurally invalid: it needs at least two breaks, "
            "one rate per break, and one exponent per segment"
        )

    if not (selected.breaks_eV[0] <= threshold_eV <= selected.breaks_eV[-1]):
        raise ValueError(
            f"Approximation is calibrated for {selected.breaks_eV[0]:.6g} to "
            f"{selected.breaks_eV[-1]:.6g} eV"
        )

    if threshold_eV == selected.breaks_eV[-1]:
        rate = area_km2 * selected.rates_per_s_per_km2[-1]
    else:
        segment_index = bisect.bisect_right(selected.breaks_eV, threshold_eV) - 1
        e0 = selected.breaks_eV[segment_index]
        r0 = selected.rates_per_s_per_km2[segment_index]
        exponent = selected.exponents[segment_index]
        try:
            rate = area_km2 * r0 * (threshold_eV / e0) ** (-exponent)
        except OverflowError as exc:
            raise ValueError("calculated approximation overflowed") from exc

    if not math.isfinite(rate):
        raise ValueError("calculated approximation is not finite")
    return rate


def _validate_sampling_range(
    e_min_eV: float,
    e_max_eV: float,
    n_points: int,
    approximation: Approximation,
) -> None:
    _require_positive_finite(e_min_eV, "e_min_eV")
    _require_positive_finite(e_max_eV, "e_max_eV")
    if e_min_eV >= e_max_eV:
        raise ValueError("e_min_eV must be smaller than e_max_eV")
    if (
        isinstance(n_points, bool)
        or not isinstance(n_points, int)
        or n_points < 2
        or n_points > MAX_SAMPLE_POINTS
    ):
        raise ValueError(
            f"n_points must be an integer between 2 and {MAX_SAMPLE_POINTS}"
        )
    if e_min_eV < approximation.breaks_eV[0] or e_max_eV > approximation.breaks_eV[-1]:
        raise ValueError(
            "Comparison range must lie inside the approximation calibration range "
            f"[{approximation.breaks_eV[0]:.6g}, "
            f"{approximation.breaks_eV[-1]:.6g}] eV"
        )


def _log_spaced_values(e_min_eV: float, e_max_eV: float, n_points: int) -> list[float]:
    log_min = math.log10(e_min_eV)
    log_span = math.log10(e_max_eV) - log_min
    values = [
        10.0 ** (log_min + index * log_span / (n_points - 1))
        for index in range(n_points)
    ]
    # 10**log10(x) can land one ulp outside [e_min_eV, e_max_eV], which would
    # fail range checks that the caller's inputs satisfy. Pin the endpoints.
    values[0] = e_min_eV
    values[-1] = e_max_eV
    return values


def _linear_quantile(sorted_values: Sequence[float], probability: float) -> float:
    if not sorted_values:
        raise ValueError("Cannot compute a quantile of an empty sequence")
    position = probability * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def _select_approximation_for_spectrum(
    spectrum: Spectrum,
    tail_model: TailModel,
    approximation: Approximation | None,
) -> Approximation:
    selected = approximation
    if selected is None:
        default = _default_approximation()
        if not _approximation_matches(default, spectrum, tail_model):
            raise ValueError(
                "The built-in approximation is tied to the bundled data.csv and "
                "tail_model='truncate'. Explicitly recalibrate an approximation "
                "for this dataset/tail model first."
            )
        selected = default

    if not _approximation_matches(selected, spectrum, tail_model):
        raise ValueError(
            "The supplied approximation was calibrated for a different dataset "
            "or tail model"
        )
    return selected


def _sampled_rates(
    csv_path: Path | str,
    e_min_eV: float,
    e_max_eV: float,
    n_points: int,
    tail_model: TailModel,
    approximation: Approximation | None,
) -> tuple[TailModel, list[tuple[float, float, float]]]:
    """Sample (energy, reference rate, approximation rate) triples.

    Shared by compare_formula and make_comparison_plot so both always use the
    same validation and the same flat 1 km^2 reference convention.
    """
    normalised_tail_model = _normalise_tail_model(tail_model)
    spectrum = load_spectrum(csv_path)
    selected = _select_approximation_for_spectrum(
        spectrum, normalised_tail_model, approximation
    )
    _validate_sampling_range(e_min_eV, e_max_eV, n_points, selected)

    samples: list[tuple[float, float, float]] = []
    for energy in _log_spaced_values(e_min_eV, e_max_eV, n_points):
        reference = _rate_above_from_spectrum(
            energy,
            spectrum,
            area_km2=1.0,
            geometry="flat",
            tail_model=normalised_tail_model,
        )
        approx = approx_rate_above_flat_km2(energy, approximation=selected)
        samples.append((energy, reference, approx))
    return normalised_tail_model, samples


def compare_formula(
    csv_path: Path | str = DEFAULT_CSV,
    e_min_eV: float = 1.0e10,
    e_max_eV: float = 1.0e20,
    n_points: int = 10_001,
    *,
    tail_model: TailModel = "truncate",
    approximation: Approximation | None = None,
) -> tuple[float, float, float]:
    """Return median, linear p95, and worst sampled relative error.

    Thresholds are uniformly spaced in log10(E). These are sampled statistics,
    not mathematical global bounds.
    """
    _, samples = _sampled_rates(
        csv_path, e_min_eV, e_max_eV, n_points, tail_model, approximation
    )
    rel_errors = sorted(
        abs(approx / reference - 1.0) for _, reference, approx in samples
    )
    return (
        statistics.median(rel_errors),
        _linear_quantile(rel_errors, 0.95),
        rel_errors[-1],
    )


def make_comparison_plot(
    output_path: Path | str,
    csv_path: Path | str = DEFAULT_CSV,
    e_min_eV: float = 1.0e10,
    e_max_eV: float = 1.0e20,
    n_points: int = 600,
    *,
    tail_model: TailModel = "truncate",
    approximation: Approximation | None = None,
) -> None:
    """Write a reference-versus-approximation plot."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError(
            "Plotting requires Matplotlib; install the project with the 'plot' extra"
        ) from exc

    normalised_tail_model, samples = _sampled_rates(
        csv_path, e_min_eV, e_max_eV, n_points, tail_model, approximation
    )
    energies = [energy for energy, _, _ in samples]
    reference_rates = [reference for _, reference, _ in samples]
    approx_rates = [approx for _, _, approx in samples]
    rel_errors_pct = [
        100.0 * (approx / reference - 1.0) for _, reference, approx in samples
    ]

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(8.0, 8.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    reference_label = (
        "CSV-truncated reference (log-log segments)"
        if normalised_tail_model == "truncate"
        else "CSV + last-segment tail reference"
    )
    ax_top.loglog(energies, reference_rates, label=reference_label, linewidth=2.2)
    ax_top.loglog(
        energies,
        approx_rates,
        "--",
        label="Calibrated piecewise power law",
        linewidth=2.0,
    )
    ax_top.set_ylabel("Ideal rate above threshold [events / s / km^2]")
    ax_top.grid(True, which="both", alpha=0.3)
    ax_top.legend()
    ax_top.set_title("Cosmic-ray threshold rate: flat 1 km^2 acceptance")

    ax_bottom.semilogx(energies, rel_errors_pct, linewidth=1.8)
    ax_bottom.axhline(0.0, linewidth=1.0, alpha=0.6)
    ax_bottom.set_xlabel("Threshold energy [eV]")
    ax_bottom.set_ylabel("Approx - reference [%]")
    ax_bottom.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(Path(output_path), dpi=180, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Command-line interface
# ---------------------------------------------------------------------------


def _positive_float_argument(value: str) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be a number") from exc
    if not math.isfinite(parsed) or parsed <= 0.0:
        raise argparse.ArgumentTypeError("must be finite and greater than zero")
    return parsed


def _breaks_argument(value: str) -> tuple[float, ...]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) < 2:
        raise argparse.ArgumentTypeError(
            "must list at least two comma-separated energies"
        )
    breaks: list[float] = []
    for part in parts:
        try:
            parsed = float(part)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(
                "must be comma-separated numbers"
            ) from exc
        if not math.isfinite(parsed) or parsed <= 0.0:
            raise argparse.ArgumentTypeError(
                "energies must be finite and greater than zero"
            )
        breaks.append(parsed)
    if any(left >= right for left, right in zip(breaks, breaks[1:])):
        raise argparse.ArgumentTypeError("energies must be strictly increasing")
    return tuple(breaks)


def _check_points_argument(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("must be an integer") from exc
    if not 2 <= parsed <= MAX_SAMPLE_POINTS:
        raise argparse.ArgumentTypeError(
            f"must be between 2 and {MAX_SAMPLE_POINTS}"
        )
    return parsed


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Ideal cosmic-ray crossing-rate estimate from a tabulated spectrum. "
            "The default reference integral is explicitly truncated at the CSV maximum."
        )
    )
    parser.add_argument(
        "threshold_eV",
        type=_positive_float_argument,
        help="Energy threshold in eV (finite and > 0)",
    )
    parser.add_argument(
        "--area-km2",
        type=_positive_float_argument,
        default=None,
        help="Geometric/effective area in km^2 (default: 1)",
    )
    parser.add_argument(
        "--geometry",
        choices=["flat", "hemisphere", "full_sky"],
        default=None,
        help=(
            "Angular convention: projected one-sided flat surface, or a "
            "direction-independent effective area over 2*pi/4*pi sr (default: flat)"
        ),
    )
    parser.add_argument(
        "--aperture-m2-sr",
        type=_positive_float_argument,
        help=(
            "Use an energy-independent aperture in m^2 sr instead of --area-km2/--geometry"
        ),
    )
    parser.add_argument(
        "--tail-model",
        choices=["truncate", "last-segment"],
        default="truncate",
        help=(
            "Above the CSV maximum: truncate explicitly (default), or extrapolate "
            "the final log-log segment to infinity"
        ),
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="Spectrum CSV using the documented header schema (default: bundled data.csv)",
    )
    parser.add_argument(
        "--recalibrate-approximation",
        action="store_true",
        help=(
            "Explicitly derive approximation coefficients from --csv and --tail-model; "
            "required for a modified dataset or non-default tail model"
        ),
    )
    parser.add_argument(
        "--approx-breaks-ev",
        type=_breaks_argument,
        default=None,
        metavar="E1,E2,...",
        help=(
            "Comma-separated, strictly increasing break energies in eV for "
            "--recalibrate-approximation (default: built-in breaks, restricted "
            "to the dataset range when necessary)"
        ),
    )
    parser.add_argument(
        "--check-formula",
        action="store_true",
        help=(
            "Print sampled median / 95th-percentile / worst relative approximation error"
        ),
    )
    parser.add_argument(
        "--check-points",
        type=_check_points_argument,
        default=10_001,
        help="Number of log-spaced thresholds used by --check-formula (default: 10001)",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        help="Write a reference-versus-approximation plot",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.aperture_m2_sr is not None and (
        args.area_km2 is not None or args.geometry is not None
    ):
        parser.error("--aperture-m2-sr cannot be combined with --area-km2 or --geometry")

    if args.approx_breaks_ev is not None and not args.recalibrate_approximation:
        parser.error("--approx-breaks-ev requires --recalibrate-approximation")

    area_km2 = args.area_km2 if args.area_km2 is not None else 1.0
    geometry: Geometry = args.geometry if args.geometry is not None else "flat"
    tail_model: TailModel = args.tail_model

    try:
        spectrum = load_spectrum(args.csv)

        approximation: Approximation | None = None
        if args.recalibrate_approximation:
            if args.approx_breaks_ev is not None:
                breaks_eV = args.approx_breaks_ev
            else:
                breaks_eV = _derive_default_breaks(spectrum, tail_model)
                if breaks_eV != APPROX_BREAKS_EV:
                    print(
                        "Note: approximation break points restricted to the "
                        "dataset range: "
                        + ", ".join(f"{value:.6g}" for value in breaks_eV)
                        + " eV"
                    )
            approximation = calibrate_piecewise_approximation(
                args.csv, breaks_eV=breaks_eV, tail_model=tail_model
            )
        else:
            # The bundled approximation applies only when the selected dataset
            # and tail model match it; a missing/broken bundled data.csv just
            # means no shortcut is available, not a hard failure.
            try:
                candidate = _default_approximation()
            except (OSError, ValueError):
                candidate = None
            if candidate is not None and _approximation_matches(
                candidate, spectrum, tail_model
            ):
                approximation = candidate

        if args.aperture_m2_sr is not None:
            reference_rate = rate_above_from_aperture(
                args.threshold_eV,
                args.aperture_m2_sr,
                csv_path=args.csv,
                tail_model=tail_model,
            )
        else:
            reference_rate = _rate_above_from_spectrum(
                args.threshold_eV,
                spectrum,
                area_km2=area_km2,
                geometry=geometry,
                tail_model=tail_model,
            )

        reference_name = (
            "CSV-truncated reference rate"
            if tail_model == "truncate"
            else "CSV + last-segment-tail rate"
        )
        print(f"{reference_name}: {reference_rate:.6g} events / s")

        can_scale_flat_approximation = (
            args.aperture_m2_sr is None and geometry == "flat"
        )
        if approximation is None:
            print(
                "Approximation             : n/a (use --recalibrate-approximation "
                "for this dataset/tail model)"
            )
        elif not can_scale_flat_approximation:
            print(
                "Approximation             : n/a (documented fit uses flat geometry)"
            )
        elif (
            approximation.breaks_eV[0]
            <= args.threshold_eV
            <= approximation.breaks_eV[-1]
        ):
            approx_rate = approx_rate_above_flat_km2(
                args.threshold_eV,
                area_km2=area_km2,
                approximation=approximation,
            )
            rel_error = abs(approx_rate / reference_rate - 1.0)
            print(f"Piecewise approximation   : {approx_rate:.6g} events / s")
            print(f"Sample-point relative error: {100.0 * rel_error:.3f} %")
        else:
            print(
                "Approximation             : n/a (threshold outside calibration range)"
            )

        if args.check_formula:
            if approximation is None:
                raise ValueError(
                    "--check-formula needs an approximation matched to the selected "
                    "dataset/tail model; add --recalibrate-approximation"
                )
            median, p95, worst = compare_formula(
                csv_path=args.csv,
                n_points=args.check_points,
                tail_model=tail_model,
                approximation=approximation,
            )
            print(
                "Sampled formula error      : "
                f"median={100.0 * median:.2f} %, "
                f"p95={100.0 * p95:.2f} %, "
                f"worst={100.0 * worst:.2f} % "
                f"({args.check_points} thresholds uniform in log10(E))"
            )

        if args.plot is not None:
            if approximation is None:
                raise ValueError(
                    "--plot needs an approximation matched to the selected dataset/"
                    "tail model; add --recalibrate-approximation"
                )
            make_comparison_plot(
                output_path=args.plot,
                csv_path=args.csv,
                tail_model=tail_model,
                approximation=approximation,
            )
            print(f"Plot written               : {args.plot}")

    except (OSError, RuntimeError, ValueError) as exc:
        parser.error(str(exc))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
