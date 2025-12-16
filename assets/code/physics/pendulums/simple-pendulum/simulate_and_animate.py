from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np

from simple_pendulum_animator import SimplePendulumAnimator
from simple_pendulum import SimplePendulum


out_dir = Path("./assets/media/physics/pendulums/simple-pendulum/videos").resolve()


def make_saa_pendulums():
    ap = SimplePendulumAnimator()

    ls = [9.81, 2, 9.81]
    thetas = [10, 10, 20]
    omegas = [0, 0, 20]
    t_f = 30
    dt = 0.02

    for l, theta0, omega0 in zip(ls, thetas, omegas):
        sp = SimplePendulum(l=l)
        results = sp.simulate(theta0, omega0, t_f, dt)

        out_path = (
            out_dir / f"simple_pend-{str(l).replace('.', '_')}-{theta0}-{omega0}.mp4"
        )
        ap.animate_single_simulation(
            results,
            "deg",
            skip_frames=3,
            out_path=out_path,
        )


def make_comparisons():
    ap = SimplePendulumAnimator()

    sp_saa = SimplePendulum(9.81)
    sp_rk4 = SimplePendulum(9.81)

    thetas = [10.0, 20.0, 30.0, 45.0, 60.0, 90.0, 120.0, 150.0, 179.0]
    omegas = [0 for _ in thetas]
    t_f = 40
    dt = 0.02

    for theta0, omega0 in zip(thetas, omegas):
        results_saa = sp_saa.simulate(theta0, omega0, t_f, dt, "saa")
        results_rk4 = sp_rk4.simulate(theta0, omega0, t_f, dt, "rk4")

        out_path = out_dir / f"compare_saa_rk4-{theta0}-{omega0}.mp4"

        ap.animate_comparison_simulations(
            results_saa,
            results_rk4,
            "deg",
            skip_frames=3,
            fps=50,
            out_path=out_path
        )

def make_separatrix_pendulum():
    ap = SimplePendulumAnimator()

    l = 9.81
    theta0 = 0.0
    omega0 = 2*180/np.pi
    t_f = 30
    dt = 0.0005

    sp = SimplePendulum(l=l)
    results = sp.simulate(theta0, omega0, t_f, dt, "rk4")

    out_path = (
        out_dir / f"separatrix_pend-{theta0:.0f}-{omega0:.0f}.mp4"
    )
    ap.animate_single_simulation(
        results,
        "deg",
        skip_frames=135,
        out_path=out_path,
    )

def make_open_pendulum():
    ap = SimplePendulumAnimator()

    l = 9.81
    theta0 = 0.0
    omega0 = 2*180/np.pi + 0.4
    t_f = 30
    dt = 0.01

    sp = SimplePendulum(l=l)
    results = sp.simulate(theta0, omega0, t_f, dt, "rk4")

    out_path = (
        out_dir / f"open_pend-{theta0:.0f}-{omega0:.0f}.mp4"
    )
    ap.animate_single_simulation(
        results,
        "rad",
        skip_frames=6,
        out_path=out_path,
    )

def _set_dark_theme_rcparams() -> None:
    rc_dict: dict[str, Any] = {
        "figure.facecolor": "black",
        "figure.titlesize": 24,
        "axes.edgecolor": "white",
        "axes.labelcolor": "white",
        "axes.titlecolor": "white",
        "axes.facecolor": "black",
        "axes.linewidth": 1.5,
        "xtick.color": "white",
        "ytick.color": "white",
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "mathtext.rm": "STIXGeneral",
        "mathtext.it": "STIXGeneral:italic",
        "mathtext.bf": "STIXGeneral:bold",
        "font.size": 18,
        "legend.frameon": True,
    }
    plt.rcParams.update(rc_dict)


def elliptic_K_agm(k: float, tol: float = 1e-12) -> float:
    if not (0.0 <= k < 1.0):
        raise ValueError("AGM method requires 0 <= k < 1")

    a = 1.0
    b = float(np.sqrt(1.0 - k * k))

    while abs(a - b) > tol:
        a, b = (a + b) / 2.0, float(np.sqrt(a * b))

    return float(np.pi / (2.0 * a))


def plot_period_and_percent_error_vs_angle(
    lengths: Iterable[float] = (5.0, 9.81, 15.0),
    *,
    g: float = 9.81,
    max_angle_deg: float = 179.0,
    n_points: int = 600,
    tol: float = 1e-12,
    out_path: str | Path | None = None,
) -> None:
    """
    Plot (1) period vs initial angle and (2) percent error of small-angle period vs exact.

    Percent error is defined as:
        100 * |T_small - T_exact| / T_exact
    """
    if g <= 0:
        raise ValueError("g must be > 0")
    if not (0.0 < max_angle_deg < 180.0):
        raise ValueError("max_angle_deg must be in (0, 180)")
    if n_points < 3:
        raise ValueError("n_points must be >= 3")

    lengths = tuple(float(l) for l in lengths)
    for l in lengths:
        if l <= 0:
            raise ValueError("all lengths must be > 0")

    if out_path is None:
        out_path = out_dir.parent / "figures/period_vs_amplitude.png"
    else:
        out_path = Path(str(out_path)).resolve()
        if out_path.suffix != ".png":
            raise ValueError(
                f"out_path must have suffix '.png', got {out_path.suffix!r}"
            )

    _set_dark_theme_rcparams()

    theta0_deg = np.linspace(0.0, max_angle_deg, n_points)
    theta0_rad = np.deg2rad(theta0_deg)

    fig, (ax_top, ax_bot) = plt.subplots(
        2,
        1,
        figsize=(11, 9),
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.6]},
    )

    ax_top.set_ylabel(r"Period $T$ (s)")
    ax_top.set_title(r"Exact nonlinear period (AGM) vs small-angle approximation")

    ax_bot.set_xlabel(r"Initial angle $\theta_0$ (deg)")
    ax_bot.set_ylabel(r"Error (%)")
    ax_bot.set_title(r"Percent error of small-angle period")

    for ax in (ax_top, ax_bot):
        ax.axhline(color="w", lw=1.0, alpha=0.6)
        ax.axvline(color="w", lw=1.0, alpha=0.6)

    # Ticks and gridlines for bottom plot
    ax_bot.minorticks_on()
    ax_bot.xaxis.set_major_locator(plt.MultipleLocator(15))
    ax_bot.xaxis.set_minor_locator(plt.MultipleLocator(5))
    ax_bot.yaxis.set_major_locator(plt.MultipleLocator(5))
    ax_bot.yaxis.set_minor_locator(plt.MultipleLocator(1))

    ax_bot.grid(
        which="major",
        color="white",
        alpha=0.15,
        linewidth=0.8,
    )
    ax_bot.grid(
        which="minor",
        color="white",
        alpha=0.06,
        linewidth=0.5,
    )

    colors = plt.cm.cool(np.linspace(0.15, 0.85, len(lengths)))

    l_ref = lengths[0]
    T_exact_ref = np.empty_like(theta0_rad)
    for i, th in enumerate(theta0_rad):
        k = float(np.sin(th / 2.0))
        K = elliptic_K_agm(k, tol=tol)
        T_exact_ref[i] = 4.0 * np.sqrt(l_ref / g) * K

    T_small_ref = float(2.0 * np.pi * np.sqrt(l_ref / g))
    err_pct = 100.0 * np.abs(T_small_ref - T_exact_ref) / T_exact_ref

    for l, color in zip(lengths, colors):
        T_exact = np.empty_like(theta0_rad)
        for i, th in enumerate(theta0_rad):
            k = float(np.sin(th / 2.0))
            K = elliptic_K_agm(k, tol=tol)
            T_exact[i] = 4.0 * np.sqrt(l / g) * K

        T_small = float(2.0 * np.pi * np.sqrt(l / g))
        T_small_arr = np.full_like(theta0_deg, T_small, dtype=float)

        ax_top.plot(
            theta0_deg,
            T_exact,
            color=color,
            lw=2.4,
            label=rf"Exact, $\ell = {l:g}\,$m",
        )
        ax_top.plot(
            theta0_deg,
            T_small_arr,
            color=color,
            lw=1.8,
            linestyle="--",
            alpha=0.75,
            label=rf"Small-angle, $\ell = {l:g}\,$m",
        )

    ax_bot.plot(
        theta0_deg,
        err_pct,
        color="#f15555",
        lw=2.4,
        label=r"Small-angle error (universal)",
    )

    leg1 = ax_top.legend(
        loc="upper left", fontsize=12, facecolor="black", edgecolor="white"
    )
    for t in leg1.get_texts():
        t.set_color("white")

    leg2 = ax_bot.legend(
        loc="upper left", fontsize=12, facecolor="black", edgecolor="white"
    )
    for t in leg2.get_texts():
        t.set_color("white")

    ax_bot.set_yticks(list(range(0, 80, 10)))


    fig.tight_layout()
    fig.savefig(out_path)




if __name__ == "__main__":
    # make_saa_pendulums()
    # make_comparisons()
    # plot_period_and_percent_error_vs_angle()
    # make_separatrix_pendulum()
    make_open_pendulum()
