from __future__ import annotations

import math
import os
import subprocess
import sys
from decimal import Decimal, localcontext
from pathlib import Path

import pytest

import cosmic_ray_rate as crr


SCRIPT = Path(crr.__file__).resolve()


def write_spectrum(
    path: Path,
    rows: list[tuple[float, float]],
    *,
    comments: tuple[str, ...] = (),
) -> Path:
    lines = [*(f"# {comment}" for comment in comments), ",".join(crr.CSV_HEADER)]
    lines.extend(f"{energy:.17g},{flux:.17g}" for energy, flux in rows)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_default_spectrum_loads() -> None:
    spectrum = crr.load_spectrum()

    assert len(spectrum.energies_eV) == 74
    assert spectrum.min_energy_eV == pytest.approx(255_087_529.14182016)
    assert spectrum.max_energy_eV == pytest.approx(3.3552442243395184e20)
    assert all(
        left < right
        for left, right in zip(spectrum.energies_eV, spectrum.energies_eV[1:])
    )


def test_reference_rate_regression_and_approximation() -> None:
    reference = crr.rate_above(1.0e15)
    approximation = crr.approx_rate_above_flat_km2(1.0e15)

    assert reference == pytest.approx(3.4742062203714448, rel=1.0e-11)
    assert abs(approximation / reference - 1.0) < 0.18


def test_calibration_recovers_a_power_law(tmp_path: Path) -> None:
    # F(E) = (E/GeV)^-2; its infinite integral is (E/GeV)^-1.
    path = write_spectrum(tmp_path / "power.csv", [(1e9, 1.0), (1e11, 1e-4)])
    approximation = crr.calibrate_piecewise_approximation(
        path, breaks_eV=(1e9, 3e9, 1e10, 1e11), tail_model="last-segment"
    )
    assert approximation.exponents == pytest.approx([1.0] * 3, rel=1e-12)
    for threshold in (1e9, 2e9, 3e9, 5e10, 1e11):
        expected = math.pi * 1e6 * (1e9 / threshold)
        actual = crr.approx_rate_above_flat_km2(threshold, approximation=approximation)
        assert actual == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_default_fit_is_continuous_and_decreasing() -> None:
    approximation = crr.calibrate_piecewise_approximation()
    assert len(approximation.exponents) == 5
    assert all(exponent > 0.0 for exponent in approximation.exponents)
    for energy in approximation.breaks_eV[1:-1]:
        left = crr.approx_rate_above_flat_km2(math.nextafter(energy, 0.0))
        right = crr.approx_rate_above_flat_km2(energy)
        assert left == pytest.approx(right, rel=1e-12, abs=0.0)


def test_power_law_integral_and_explicit_tail(tmp_path: Path) -> None:
    # f(E) = E^-2 on [1e9, 1e10] in the CSV units.
    csv_path = write_spectrum(
        tmp_path / "power_law.csv",
        [(1.0e9, 1.0e-18), (1.0e10, 1.0e-20)],
    )

    truncated = crr.integrated_intensity_above(2.0e9, csv_path)
    extrapolated = crr.integrated_intensity_above(
        2.0e9, csv_path, tail_model="last-segment"
    )

    # Integral in eV is 1/E_th - 1/E_max; multiply by 1e-9 for dE_GeV.
    assert truncated == pytest.approx(4.0e-19, rel=1.0e-11, abs=0.0)
    assert extrapolated == pytest.approx(5.0e-19, rel=1.0e-11, abs=0.0)


@pytest.mark.parametrize("gamma", [-3.7, -2.0, -1.0000000001, -1.0, -0.9999999999, 0.0, 2.0])
def test_power_laws_against_high_precision_integrals(tmp_path: Path, gamma: float) -> None:
    # F(E) = (E / GeV)**gamma per m^2 sr s GeV on [1, 100] GeV.
    path = write_spectrum(tmp_path / "power.csv", [(1e9, 1.0), (1e11, 100.0**gamma)])
    with localcontext() as ctx:
        ctx.prec = 60
        q = Decimal(str(gamma)) + 1
        for threshold_GeV in (1.0, 2.5, 10.0, 99.0):
            low, high = Decimal(str(threshold_GeV)), Decimal(100)
            expected = (high / low).ln() if q == 0 else (high**q - low**q) / q
            actual = crr.integrated_intensity_above(threshold_GeV * 1e9, path)
            assert actual == pytest.approx(float(expected), rel=1e-12, abs=0.0)
        if gamma < -1.0:
            # Above the table, the infinite-tail integral has a known slope.
            expected_tail = float(-(Decimal(200)**q) / q)
            actual = crr.integrated_intensity_above(2e11, path, tail_model="last-segment")
            # Near -1 the tail is ill-conditioned against rounding the input flux.
            tolerance = 1e-5 if abs(gamma + 1) < 1e-8 else 1e-12
            assert actual == pytest.approx(expected_tail, rel=tolerance, abs=0.0)
        else:
            with pytest.raises(ValueError, match="does not converge"):
                crr.integrated_intensity_above(2e11, path, tail_model="last-segment")


def test_integral_resolves_last_float_below_endpoint() -> None:
    spectrum = crr.load_spectrum()
    upper = spectrum.max_energy_eV
    lower = math.nextafter(upper, 0.0)
    # Over one ulp the flux is constant to floating-point precision.
    expected = spectrum.fluxes_per_m2_sr_s_GeV[-1] * (upper - lower) * 1e-9
    assert crr.integrated_intensity_above(lower) == pytest.approx(expected, rel=1e-12, abs=0.0)


def test_rate_is_continuous_decreasing_and_has_the_correct_derivative() -> None:
    spectrum = crr.load_spectrum()
    thresholds = crr._log_spaced_values(spectrum.min_energy_eV, spectrum.max_energy_eV, 101)
    rates = [crr.rate_above(e) for e in thresholds]
    assert all(left > right for left, right in zip(rates, rates[1:]))
    for i, energy in enumerate(spectrum.energies_eV[1:-1], start=1):
        delta = energy * 1e-6
        left, right = crr.rate_above(energy - delta), crr.rate_above(energy + delta)
        flux = spectrum.fluxes_per_m2_sr_s_GeV[i]
        # Fundamental theorem: dR/dE_eV = -pi * 10^-3 * F(E) for 1 km^2.
        assert (left - right) / (2 * delta) == pytest.approx(math.pi * 1e-3 * flux, rel=1e-5, abs=0.0)


def test_csv_endpoint_is_zero_only_for_the_truncated_definition() -> None:
    spectrum = crr.load_spectrum()

    assert crr.rate_above(spectrum.max_energy_eV, tail_model="truncate") == 0.0
    with pytest.raises(ValueError, match="exceeds the CSV maximum"):
        crr.rate_above(1.01 * spectrum.max_energy_eV, tail_model="truncate")

    endpoint_tail = crr.rate_above(
        spectrum.max_energy_eV, tail_model="last-segment"
    )
    assert endpoint_tail == pytest.approx(5.193586955416152e-12, rel=1.0e-11, abs=0.0)
    assert crr.rate_above(1.0e20, tail_model="last-segment") > crr.rate_above(1.0e20)


def test_nonconvergent_last_segment_is_rejected(tmp_path: Path) -> None:
    csv_path = write_spectrum(
        tmp_path / "nonconvergent.csv",
        [(1.0, 1.0), (10.0, 10.0**-0.5)],
    )

    with pytest.raises(ValueError, match="does not converge"):
        crr.integrated_intensity_above(2.0, csv_path, tail_model="last-segment")


def test_loader_uses_schema_and_sorts_rows(tmp_path: Path) -> None:
    rows = [(1.0, 10.0), (2.0, 4.0), (4.0, 1.0)]
    first = write_spectrum(tmp_path / "first.csv", rows, comments=("first",))
    second = write_spectrum(
        tmp_path / "second.csv", list(reversed(rows)), comments=("different",)
    )

    spectrum_a = crr.load_spectrum(first)
    spectrum_b = crr.load_spectrum(second)

    assert spectrum_b.energies_eV == (1.0, 2.0, 4.0)
    assert spectrum_b.fluxes_per_m2_sr_s_GeV == (10.0, 4.0, 1.0)
    assert spectrum_a == spectrum_b


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("wrong,header\n1,2\n2,1\n", "expected CSV header"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1,0\n2,1\n", "greater than zero"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1,nan\n2,1\n", "must be numeric"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1,2\n1,1\n", "duplicate energy"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1,2,3\n2,1\n", "expected 2 columns"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1_0,2.5\n2,1\n", "must be numeric"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n١٢e12,3.0\n2,1\n", "must be numeric"),
        ("energy_eV,flux_per_m2_sr_s_GeV\n1,infinity\n2,1\n", "must be numeric"),
        ('energy_eV,flux_per_m2_sr_s_GeV\n1,"2\n2,1\n', "invalid CSV row"),
    ],
)
def test_loader_rejects_invalid_csv(tmp_path: Path, contents: str, message: str) -> None:
    path = tmp_path / "bad.csv"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        crr.load_spectrum(path)


def test_bom_with_leading_comment_loads(tmp_path: Path) -> None:
    content = (
        "# provenance comment\n"
        + ",".join(crr.CSV_HEADER)
        + "\n1.0,2.0\n2.0,1.0\n"
    )
    path = tmp_path / "bom.csv"
    path.write_bytes(b"\xef\xbb\xbf" + content.encode("utf-8"))

    spectrum = crr.load_spectrum(path)

    assert spectrum.energies_eV == (1.0, 2.0)


def test_cache_serves_new_content_after_same_size_same_mtime_rewrite(
    tmp_path: Path,
) -> None:
    header = ",".join(crr.CSV_HEADER)
    path = tmp_path / "spectrum.csv"
    path.write_text(f"{header}\n1.0,8.0\n2.0,4.0\n", encoding="utf-8")
    first = crr.load_spectrum(path)
    stat = path.stat()

    # Filesystem metadata must not determine which physical spectrum is used.
    path.write_text(f"{header}\n1.0,6.0\n2.0,3.0\n", encoding="utf-8")
    os.utime(path, ns=(stat.st_atime_ns, stat.st_mtime_ns))
    second = crr.load_spectrum(path)

    assert first.fluxes_per_m2_sr_s_GeV == (8.0, 4.0)
    assert second.fluxes_per_m2_sr_s_GeV == (6.0, 3.0)


def test_compare_formula_reaches_custom_range_endpoints() -> None:
    # 2195841772600371.0 does not survive a 10**log10(x) round trip; before the
    # endpoint clamp this raised "Approximation is calibrated for ..." even
    # though the requested range equalled the calibration range.
    upper = 2195841772600371.0
    approximation = crr.calibrate_piecewise_approximation(
        breaks_eV=(1.0e10, 1.0e14, upper)
    )

    median, p95, worst = crr.compare_formula(
        e_min_eV=1.0e10,
        e_max_eV=upper,
        n_points=101,
        approximation=approximation,
    )

    assert 0.0 <= median <= p95 <= worst


def test_structurally_invalid_approximation_raises_value_error() -> None:
    empty = crr.Approximation(
        breaks_eV=(),
        rates_per_s_per_km2=(),
        exponents=(),
        spectrum=crr.load_spectrum(),
        tail_model="truncate",
    )

    with pytest.raises(ValueError, match="structurally invalid"):
        crr.approx_rate_above_flat_km2(1.0e15, approximation=empty)


@pytest.mark.parametrize("value", [0.0, -1.0, math.inf, -math.inf, math.nan])
def test_rate_rejects_invalid_area(value: float) -> None:
    with pytest.raises(ValueError, match="area_km2"):
        crr.rate_above(1.0e15, area_km2=value)


@pytest.mark.parametrize("value", [0.0, -1.0, math.inf, -math.inf, math.nan])
def test_integrator_rejects_invalid_threshold(value: float) -> None:
    with pytest.raises(ValueError, match="threshold_eV"):
        crr.integrated_intensity_above(value)


def test_scalar_aperture_and_angular_conventions() -> None:
    flat = crr.rate_above(1.0e15, geometry="flat")
    hemisphere = crr.rate_above(1.0e15, geometry="hemisphere")
    full_sky = crr.rate_above(1.0e15, geometry="full_sky")
    aperture = crr.rate_above_from_aperture(1.0e15, math.pi * 1.0e6)

    assert hemisphere == pytest.approx(2.0 * flat)
    assert full_sky == pytest.approx(4.0 * flat)
    assert aperture == pytest.approx(flat)

    with pytest.raises(ValueError, match="aperture_m2_sr"):
        crr.rate_above_from_aperture(1.0e15, 0.0)


def test_custom_data_requires_and_supports_explicit_recalibration(tmp_path: Path) -> None:
    default = crr.load_spectrum()
    doubled_rows = [
        (energy, 2.0 * flux)
        for energy, flux in zip(
            default.energies_eV, default.fluxes_per_m2_sr_s_GeV
        )
    ]
    custom = write_spectrum(tmp_path / "doubled.csv", doubled_rows)
    calibrated = crr.calibrate_piecewise_approximation(custom)

    assert crr.rate_above(1.0e15, csv_path=custom) == pytest.approx(
        2.0 * crr.rate_above(1.0e15)
    )
    assert crr.approx_rate_above_flat_km2(
        1.0e15, approximation=calibrated
    ) == pytest.approx(2.0 * crr.approx_rate_above_flat_km2(1.0e15))

    with pytest.raises(ValueError, match="Explicitly recalibrate"):
        crr.compare_formula(csv_path=custom, n_points=11)

    median, p95, worst = crr.compare_formula(
        csv_path=custom, n_points=101, approximation=calibrated
    )
    assert 0.0 <= median <= p95 <= worst < 0.20


@pytest.mark.parametrize("n_points", [0, 1, 1.5, True, 10**9])
def test_comparison_requires_a_bounded_integer_sample_count(n_points: object) -> None:
    with pytest.raises(ValueError, match="integer between 2 and"):
        crr.compare_formula(n_points=n_points)  # type: ignore[arg-type]


def test_sampled_error_statistics_converge() -> None:
    coarse = crr.compare_formula(n_points=1_001)
    fine = crr.compare_formula(n_points=10_001)
    # At least the previous six-piece fit's accuracy, across the same range.
    assert 0.0 <= fine[0] <= fine[1] <= fine[2]
    assert fine[0] < 0.0544
    assert fine[1] < 0.15
    assert fine[2] < 0.18
    assert coarse == pytest.approx(fine, abs=0.002)


def test_cli_reports_invalid_input_without_a_traceback() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "1e15", "--area-km2", "-1"],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert "greater than zero" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_cli_rejects_underflow_without_a_traceback(capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as excinfo:
        crr.main(["1e20", "--area-km2", "1e-320"])
    assert excinfo.value.code == 2
    assert "floating-point range" in capsys.readouterr().err


def test_cli_does_not_silently_use_default_fit_for_custom_data(tmp_path: Path) -> None:
    default = crr.load_spectrum()
    custom = write_spectrum(
        tmp_path / "scaled.csv",
        [
            (energy, 3.0 * flux)
            for energy, flux in zip(
                default.energies_eV, default.fluxes_per_m2_sr_s_GeV
            )
        ],
    )

    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "1e15", "--csv", str(custom)],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0
    assert "Approximation             : n/a" in completed.stdout
    assert "recalibrate-approximation" in completed.stdout


def test_cli_survives_missing_default_csv_with_custom_data(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    default = crr.load_spectrum()
    custom = write_spectrum(
        tmp_path / "scaled.csv",
        [
            (energy, 2.0 * flux)
            for energy, flux in zip(
                default.energies_eV, default.fluxes_per_m2_sr_s_GeV
            )
        ],
    )
    monkeypatch.setattr(crr, "DEFAULT_CSV", tmp_path / "missing.csv")

    assert crr.main(["1e15", "--csv", str(custom)]) == 0
    assert "Approximation             : n/a" in capsys.readouterr().out

    assert (
        crr.main(["1e15", "--csv", str(custom), "--recalibrate-approximation"])
        == 0
    )
    assert "Piecewise approximation" in capsys.readouterr().out


def test_cli_recalibrates_narrow_dataset_with_derived_breaks(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    custom = write_spectrum(
        tmp_path / "narrow.csv",
        [(1.0e11, 1.0e-2), (1.0e15, 1.0e-10), (1.0e19, 1.0e-20)],
    )

    assert (
        crr.main(["1e15", "--csv", str(custom), "--recalibrate-approximation", "--check-formula"])
        == 0
    )
    output = capsys.readouterr().out
    assert "Note: approximation break points restricted" in output
    assert "Piecewise approximation" in output
    assert "Sampled formula error" in output

    assert (
        crr.main(
            [
                "1e13",
                "--csv",
                str(custom),
                "--recalibrate-approximation",
                "--approx-breaks-ev",
                "1e12,1e14,1e16",
                "--check-formula",
            ]
        )
        == 0
    )
    assert "Piecewise approximation" in capsys.readouterr().out


def test_default_approximation_tracks_edits(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spectrum = crr.load_spectrum()
    rows = list(zip(spectrum.energies_eV, spectrum.fluxes_per_m2_sr_s_GeV))
    path = write_spectrum(tmp_path / "default.csv", rows)
    monkeypatch.setattr(crr, "DEFAULT_CSV", path)
    before = crr.approx_rate_above_flat_km2(1e15)
    write_spectrum(path, [(e, 2 * f) for e, f in rows])
    assert crr.approx_rate_above_flat_km2(1e15) == pytest.approx(2 * before)


def test_narrow_calibration_plot(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    approximation = crr.calibrate_piecewise_approximation(breaks_eV=(1e12, 1e14, 1e16))
    output = tmp_path / "comparison.png"
    crr.make_comparison_plot(output, approximation=approximation, n_points=51)
    assert output.is_file()


def test_cli_breaks_flag_requires_recalibration() -> None:
    with pytest.raises(SystemExit) as excinfo:
        crr.main(["1e15", "--approx-breaks-ev", "1e12,1e14"])
    assert excinfo.value.code == 2


def test_cli_formula_check_needs_a_matching_approximation(tmp_path: Path) -> None:
    custom = write_spectrum(
        tmp_path / "small.csv",
        [(1.0e10, 1.0), (1.0e12, 1.0e-4), (1.0e20, 1.0e-20), (2.0e20, 1.0e-22)],
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "1e12",
            "--csv",
            str(custom),
            "--check-formula",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert "add --recalibrate-approximation" in completed.stderr
    assert "Traceback" not in completed.stderr
