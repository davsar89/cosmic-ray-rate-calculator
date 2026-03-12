from __future__ import annotations

import argparse
import math
from pathlib import Path


DEFAULT_CSV = Path(__file__).with_name("data.csv")

# Approximate rate for a flat horizontal 1 km^2 detector exposed to an isotropic
# downward flux. Units: events / s / km^2 above E_th.
APPROX_BREAKS_EV = [
    1.0e10,
    1.0e12,
    2.0e15,
    5.0e16,
    3.0e18,
    3.0e19,
    1.0e20,
]
APPROX_RATES_PER_S_PER_KM2 = [
    4.061354085353714e8,
    3.2976226042660396e5,
    1.029360329532377,
    1.8695676052064019e-3,
    2.1958700384821697e-7,
    2.6415245285191667e-9,
    1.0436638137472337e-10,
]
APPROX_EXPONENTS = [
    1.5452349516439396,
    1.6678532251023768,
    1.9606179347727197,
    2.210236356647733,
    1.9197519852386478,
    2.6837847870606364,
]


def load_spectrum(csv_path: Path) -> tuple[list[float], list[float]]:
    col1: list[float] = []
    col2: list[float] = []

    for raw_line in csv_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        left, right = [part.strip() for part in line.split(",", 1)]
        col1.append(float(left))
        col2.append(float(right))

    if len(col1) < 2:
        raise ValueError(f"No numeric data found in {csv_path}")

    # The header in data.csv says "F, E", but the numeric content behaves like
    # "E (eV), F ((m^2 sr s GeV)^-1)". Use the monotonic trend in the data.
    col1_up = all(a < b for a, b in zip(col1, col1[1:]))
    col2_up = all(a < b for a, b in zip(col2, col2[1:]))
    col1_down = all(a > b for a, b in zip(col1, col1[1:]))
    col2_down = all(a > b for a, b in zip(col2, col2[1:]))

    if col1_up and col2_down:
        energies_eV, fluxes = col1, col2
    elif col2_up and col1_down:
        energies_eV, fluxes = col2, col1
    else:
        raise ValueError(
            "Could not infer which CSV column is energy. Expected one column to "
            "increase with row number and the other to decrease."
        )

    return energies_eV, fluxes


def _integrate_power_law_segment(
    e_lo_eV: float,
    e_hi_eV: float,
    f_lo: float,
    f_hi: float,
) -> float:
    gamma = math.log(f_hi / f_lo) / math.log(e_hi_eV / e_lo_eV)
    coeff = f_lo / (e_lo_eV ** gamma)

    if math.isclose(gamma, -1.0, rel_tol=0.0, abs_tol=1e-12):
        return coeff * math.log(e_hi_eV / e_lo_eV)

    return coeff * (e_hi_eV ** (gamma + 1.0) - e_lo_eV ** (gamma + 1.0)) / (gamma + 1.0)


def integrated_intensity_above(
    threshold_eV: float,
    csv_path: Path = DEFAULT_CSV,
) -> float:
    """
    Returns the exact threshold-integrated intensity in units:
        m^-2 sr^-1 s^-1

    The CSV flux is differential in GeV, while the energy axis is in eV, so the
    final integral carries an extra factor dE_GeV = dE_eV / 1e9.
    """
    energies_eV, fluxes = load_spectrum(csv_path)

    if not (energies_eV[0] <= threshold_eV <= energies_eV[-1]):
        raise ValueError(
            f"threshold_eV must lie inside [{energies_eV[0]:.6g}, {energies_eV[-1]:.6g}] eV"
        )

    integral_eV = 0.0
    for e1, e2, f1, f2 in zip(energies_eV, energies_eV[1:], fluxes, fluxes[1:]):
        if threshold_eV >= e2:
            continue

        seg_lo = max(threshold_eV, e1)
        seg_hi = e2
        if seg_hi <= seg_lo:
            continue

        if seg_lo == e1:
            f_lo = f1
        else:
            gamma = math.log(f2 / f1) / math.log(e2 / e1)
            f_lo = f1 * (seg_lo / e1) ** gamma

        integral_eV += _integrate_power_law_segment(seg_lo, seg_hi, f_lo, f2)

    return integral_eV * 1.0e-9


def rate_above(
    threshold_eV: float,
    area_km2: float = 1.0,
    geometry: str = "flat",
    csv_path: Path = DEFAULT_CSV,
) -> float:
    """
    Returns the event rate above threshold in events / s.

    geometry='flat':
        flat horizontal detector, downward isotropic flux -> multiply by pi
    geometry='hemisphere':
        integrated over 2 pi sr with no cos(theta) projection
    geometry='full_sky':
        integrated over 4 pi sr with no projection
    """
    intensity = integrated_intensity_above(threshold_eV, csv_path=csv_path)

    if geometry == "flat":
        solid_angle_factor = math.pi
    elif geometry == "hemisphere":
        solid_angle_factor = 2.0 * math.pi
    elif geometry == "full_sky":
        solid_angle_factor = 4.0 * math.pi
    else:
        raise ValueError("geometry must be one of: 'flat', 'hemisphere', 'full_sky'")

    area_m2 = area_km2 * 1.0e6
    return intensity * solid_angle_factor * area_m2


def approx_rate_above_flat_km2(threshold_eV: float, area_km2: float = 1.0) -> float:
    """
    Piecewise power-law approximation to rate_above(..., geometry='flat').

    Validity range:
        1e10 eV <= threshold_eV <= 1e20 eV

    Typical error is a few to ~15 percent versus the exact interpolation-based
    integral from data.csv, with the worst case staying below ~20 percent over
    this range.
    """
    if not (APPROX_BREAKS_EV[0] <= threshold_eV <= APPROX_BREAKS_EV[-1]):
        raise ValueError(
            f"Approximation is calibrated for {APPROX_BREAKS_EV[0]:.0e} to "
            f"{APPROX_BREAKS_EV[-1]:.0e} eV"
        )

    for i, exponent in enumerate(APPROX_EXPONENTS):
        e0 = APPROX_BREAKS_EV[i]
        e1 = APPROX_BREAKS_EV[i + 1]
        if e0 <= threshold_eV <= e1:
            r0 = APPROX_RATES_PER_S_PER_KM2[i]
            return area_km2 * r0 * (threshold_eV / e0) ** (-exponent)

    return area_km2 * APPROX_RATES_PER_S_PER_KM2[-1]


def compare_formula(
    csv_path: Path = DEFAULT_CSV,
    e_min_eV: float = 1.0e10,
    e_max_eV: float = 1.0e20,
    n_points: int = 200,
) -> tuple[float, float, float]:
    rel_errors: list[float] = []
    for i in range(n_points):
        t = i / (n_points - 1)
        energy = 10.0 ** (
            math.log10(e_min_eV) + t * (math.log10(e_max_eV) - math.log10(e_min_eV))
        )
        exact = rate_above(energy, geometry="flat", csv_path=csv_path)
        approx = approx_rate_above_flat_km2(energy)
        rel_errors.append(abs(approx / exact - 1.0))

    rel_errors.sort()
    median = rel_errors[len(rel_errors) // 2]
    p95 = rel_errors[int(0.95 * (len(rel_errors) - 1))]
    worst = rel_errors[-1]
    return median, p95, worst


def make_comparison_plot(
    output_path: Path,
    csv_path: Path = DEFAULT_CSV,
    e_min_eV: float = 1.0e10,
    e_max_eV: float = 1.0e20,
    n_points: int = 300,
) -> None:
    import matplotlib.pyplot as plt

    energies: list[float] = []
    exact_rates: list[float] = []
    approx_rates: list[float] = []
    rel_errors_pct: list[float] = []

    for i in range(n_points):
        t = i / (n_points - 1)
        energy = 10.0 ** (
            math.log10(e_min_eV) + t * (math.log10(e_max_eV) - math.log10(e_min_eV))
        )
        exact = rate_above(energy, geometry="flat", csv_path=csv_path)
        approx = approx_rate_above_flat_km2(energy)

        energies.append(energy)
        exact_rates.append(exact)
        approx_rates.append(approx)
        rel_errors_pct.append(100.0 * (approx / exact - 1.0))

    fig, (ax_top, ax_bottom) = plt.subplots(
        2,
        1,
        figsize=(8.0, 8.0),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]},
    )

    ax_top.loglog(energies, exact_rates, label="Exact: CSV + log-log integration", linewidth=2.2)
    ax_top.loglog(
        energies,
        approx_rates,
        "--",
        label="Analytical: piecewise power law",
        linewidth=2.0,
    )
    ax_top.set_ylabel("Rate above threshold [events / s / km^2]")
    ax_top.grid(True, which="both", alpha=0.3)
    ax_top.legend()
    ax_top.set_title("Cosmic-ray rate above threshold")

    ax_bottom.semilogx(energies, rel_errors_pct, linewidth=1.8)
    ax_bottom.axhline(0.0, color="black", linewidth=1.0, alpha=0.6)
    ax_bottom.set_xlabel("Threshold energy [eV]")
    ax_bottom.set_ylabel("Approx - Exact [%]")
    ax_bottom.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cosmic-ray rate above a threshold, based on data.csv"
    )
    parser.add_argument("threshold_eV", type=float, help="Energy threshold in eV")
    parser.add_argument(
        "--area-km2",
        type=float,
        default=1.0,
        help="Detector area in km^2 (default: 1)",
    )
    parser.add_argument(
        "--geometry",
        choices=["flat", "hemisphere", "full_sky"],
        default="flat",
        help="Angular integration convention for the exact result (default: flat)",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help="Path to the CSV file (default: ./data.csv)",
    )
    parser.add_argument(
        "--check-formula",
        action="store_true",
        help="Print median / 95th-percentile / worst relative error of the approximation",
    )
    parser.add_argument(
        "--plot",
        type=Path,
        help="Write a log-log comparison plot (flat geometry, 1e10 to 1e20 eV)",
    )
    args = parser.parse_args()

    exact_rate = rate_above(
        args.threshold_eV,
        area_km2=args.area_km2,
        geometry=args.geometry,
        csv_path=args.csv,
    )
    print(f"Exact rate   : {exact_rate:.12g} events / s")

    if args.geometry == "flat" and APPROX_BREAKS_EV[0] <= args.threshold_eV <= APPROX_BREAKS_EV[-1]:
        approx_rate = approx_rate_above_flat_km2(args.threshold_eV, area_km2=args.area_km2)
        rel_error = abs(approx_rate / exact_rate - 1.0)
        print(f"Approx rate  : {approx_rate:.12g} events / s")
        print(f"Rel. error   : {100.0 * rel_error:.3f} %")
    elif args.geometry == "flat":
        print("Approx rate  : n/a (formula calibrated for 1e10 to 1e20 eV)")

    if args.check_formula:
        median, p95, worst = compare_formula(csv_path=args.csv)
        print(f"Formula error: median={100.0 * median:.2f} %, p95={100.0 * p95:.2f} %, worst={100.0 * worst:.2f} %")

    if args.plot is not None:
        make_comparison_plot(output_path=args.plot, csv_path=args.csv)
        print(f"Plot written  : {args.plot}")


if __name__ == "__main__":
    main()
