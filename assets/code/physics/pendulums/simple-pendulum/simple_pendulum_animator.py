from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal, TypeAlias

import matplotlib

matplotlib.use("Agg")

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.text import Text
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from simple_pendulum import SimplePendulumResult


AngleUnit: TypeAlias = Literal["deg", "rad"]


@dataclass(frozen=True)
class SimplePendulumFigAxes:
    """
    Bundle of Matplotlib objects used by the simple pendulum animations.

    Layout
    ------
    Single simulation:
        ax1: pendulum geometry (standalone)
        ax2: phase portrait (omega vs theta)
        ax3: angle vs time
        time_text: time annotation drawn on ax3 (top-left)

    Comparison:
        ax1: phase portrait for simulation 1 with inset pendulum overlay
        ax2: phase portrait for simulation 2 with inset pendulum overlay
        ax3: angle vs time (both curves)
        time_text: time annotation drawn on ax3 (top-left)

    Notes
    -----
    `ax1` meaning depends on whether comparison mode is enabled.
    """

    fig: Figure
    ax1: Axes
    ax2: Axes
    ax3: Axes
    time_text: Text


class PendulumAnimator:
    """
    Base class for pendulum animation helpers.

    This class centralizes shared styling and numeric utilities so that
    specialized pendulum animators (simple pendulum, driven pendulum, etc.)
    can inherit a consistent look and behavior.
    """

    def __init__(self) -> None:
        self._set_dark_theme()

    def _set_dark_theme(self) -> None:
        """
        Apply a dark theme and typography settings via Matplotlib rcParams.
        """
        rc_dict = {
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
        }
        plt.rcParams.update(rc_dict)

    def _padded_limits(
        self, vmin: float, vmax: float, pad_frac: float = 0.08
    ) -> tuple[float, float]:
        """
        Compute [min, max] limits padded by a fraction of the data span.
        """
        if vmin == vmax:
            pad = 1.0 if vmin == 0 else abs(vmin) * pad_frac
            return vmin - pad, vmax + pad
        span = vmax - vmin
        pad = span * pad_frac
        return vmin - pad, vmax + pad

    def _make_pendulum_inset(
        self,
        ax_phase: Axes,
        color: str,
        *,
        width: str = "100%",
        height: str = "100%",
        loc: str = "center",
        borderpad: float = 0.9,
        bob_size: float = 22,
        rod_lw: float = 3.5,
    ) -> tuple[Axes, Line2D, Line2D]:
        ax_inset = inset_axes(
            ax_phase,
            width=width,
            height=height,
            loc=loc,
            borderpad=borderpad,
        )
        ax_inset.set_facecolor("none")
        ax_inset.set_xticks([])
        ax_inset.set_yticks([])
        for spine in ax_inset.spines.values():
            spine.set_visible(False)

        ax_inset.set_xlim(-1.2, 1.2)
        ax_inset.set_ylim(-1.2, 1.2)
        ax_inset.set_aspect("equal", adjustable="box")

        (rod,) = ax_inset.plot([], [], color="w", lw=rod_lw, zorder=10)
        (bob,) = ax_inset.plot(
            [],
            [],
            marker="o",
            linestyle="None",
            color=color,
            markersize=bob_size,
            markeredgecolor="#bbbbbb",
            markeredgewidth=1.5,
            zorder=11,
        )
        return ax_inset, rod, bob

    def _make_pendulum_axes(
        self,
        ax: Axes,
        color: str,
        *,
        bob_size: float = 22,
        rod_lw: float = 3.2,
    ) -> tuple[Line2D, Line2D]:
        ax.set_title("Pendulum")
        ax.set_xticks([])
        ax.set_yticks([])

        ax.set_xlim(-1.2, 1.2)
        ax.set_ylim(-1.2, 1.2)
        ax.set_aspect("equal", adjustable="box")

        (rod,) = ax.plot([], [], color="w", lw=rod_lw, zorder=10)
        (bob,) = ax.plot(
            [],
            [],
            marker="o",
            linestyle="None",
            color=color,
            markersize=bob_size,
            markeredgecolor="#bbbbbb",
            markeredgewidth=1.5,
            zorder=11,
        )
        return rod, bob

    def _make_time_text(self, ax: Axes) -> Text:
        return ax.text(
            -0.115,
            1.224,
            "",
            transform=ax.transAxes,
            ha="left",
            va="top",
            color="white",
            fontsize=16,
            bbox=dict(
                boxstyle="round,pad=0.6",
                facecolor="black",
                edgecolor="white",
                linewidth=1.2,
            ),
        )

    def _format_phase_axes(self, ax: Axes, *, angle_unit: AngleUnit, title: str) -> None:
        ax.set_title(title)
        ax.set_xlabel(f"Angle ({angle_unit})")
        ax.set_ylabel(f"Angular Speed ({angle_unit}/s)")
        ax.axhline(color="w", lw=0.7)
        ax.axvline(color="w", lw=0.7)

    def _format_angle_time_axes(self, ax: Axes, *, angle_unit: AngleUnit, title: str) -> None:
        ax.set_title(title)
        ax.set_xlabel("Time (s)")
        ax.set_ylabel(f"Angle ({angle_unit})")
        ax.axhline(color="w", lw=0.7)

    def _suptitle_single(
        self,
        *,
        l: float,
        theta0: float,
        omega0: float,
        period: float | None,
        angle_unit: AngleUnit,
    ) -> str:
        sep = r"$\qquad$"
        pieces = [
            rf"$l = {l:g}m$",
            rf"$\theta_0 = {theta0:.1f}\ \mathrm{{{angle_unit}}}$",
            rf"$\omega_0 = {omega0:.1f}\ \mathrm{{{angle_unit}}}/\mathrm{{s}}$",
        ]
        if period is not None:
            pieces.append(f"$T=${period:.2f}s")
        return sep.join(pieces)

    def _suptitle_compare(
        self,
        *,
        l1: float,
        theta10: float,
        omega10: float,
        period1: float | None,
        l2: float,
        theta20: float,
        omega20: float,
        period2: float | None,
        angle_unit: AngleUnit,
    ) -> str:
        sep = r"$\qquad$"
        suptitle1_pieces = [
                rf"$l_1 = {l1:g}m$",
                rf"$\theta_{{1,0}} = {theta10:.1f}\ \mathrm{{{angle_unit}}}$",
                rf"$\omega_{{1,0}} = {omega10:.1f}\ \mathrm{{{angle_unit}}}/\mathrm{{s}}$",
            ]
        if period1 is not None:
            suptitle1_pieces.append(f"$T_1=${period1:.2f}s")
        line1 = sep.join(suptitle1_pieces)

        suptitle2_pieces = [
                rf"$l_2 = {l2:g}m$",
                rf"$\theta_{{2,0}} = {theta20:.1f}\ \mathrm{{{angle_unit}}}$",
                rf"$\omega_{{2,0}} = {omega20:.1f}\ \mathrm{{{angle_unit}}}/\mathrm{{s}}$",
            ]
        if period2 is not None:
            suptitle2_pieces.append(f"$T_2=${period2:.2f}s")
        line2 = sep.join(suptitle2_pieces)
        return line1 + "\n" + line2


class SimplePendulumAnimator(PendulumAnimator):
    """
    Animator for `SimplePendulumResult` using Matplotlib FuncAnimation.

    Two modes:
    - Single: pendulum + phase portrait + angle vs time (all separate axes)
    - Compare: phase1 with inset overlay + phase2 with inset overlay + angle vs time
    """

    def _get_figs_and_axes(
        self,
        simulation_results: SimplePendulumResult,
        simulation_results_compare: SimplePendulumResult | None = None,
        angle_unit: AngleUnit = "deg",
    ) -> SimplePendulumFigAxes:
        l = simulation_results.l
        algorithm = simulation_results.algorithm
        time = simulation_results.time
        period = simulation_results.period

        if angle_unit == "deg":
            theta = simulation_results.theta_deg
            omega = simulation_results.omega_deg
        else:
            theta = simulation_results.theta_rad
            omega = simulation_results.omega_rad

        fig = plt.figure(figsize=(10, 9))

        compare_mode = simulation_results_compare is not None

        gs = fig.add_gridspec(
            2,
            2,
            height_ratios=[1.0, 1.0],
            width_ratios=[1.0, 1.0],
            hspace=0.35,
            wspace=0.30,
        )

        ax_top_left = fig.add_subplot(gs[0, 0])
        ax_top_right = fig.add_subplot(gs[0, 1])
        ax_bottom = fig.add_subplot(gs[1, :])

        ax_top_left.set_box_aspect(1.0)
        ax_top_right.set_box_aspect(1.0)

        if not compare_mode:
            suptitle = self._suptitle_single(
                l=float(l),
                theta0=float(theta[0]),
                omega0=float(omega[0]),
                period=period,
                angle_unit=angle_unit,
            )
            fig.suptitle(suptitle, color="w")

            self._format_phase_axes(
                ax_top_right,
                angle_unit=angle_unit,
                title="Angular Speed vs. Angle",
            )
            self._format_angle_time_axes(
                ax_bottom,
                angle_unit=angle_unit,
                title="Angle vs. Time",
            )

            time_text = self._make_time_text(ax_bottom)

            x0, x1 = self._padded_limits(float(np.min(theta)), float(np.max(theta)), pad_frac=0.08)
            y0, y1 = self._padded_limits(float(np.min(omega)), float(np.max(omega)), pad_frac=0.08)
            mn, mx = min(x0, y0), max(x1, y1)
            ax_top_right.set_xlim(mn, mx)
            ax_top_right.set_ylim(mn, mx)
            ax_top_right.set_aspect("equal", adjustable="box")

            x0, x1 = self._padded_limits(float(np.min(time)), float(np.max(time)), pad_frac=0.03)
            y0, y1 = self._padded_limits(float(np.min(theta)), float(np.max(theta)), pad_frac=0.03)
            ax_bottom.set_xlim(x0, x1)
            ax_bottom.set_ylim(y0, y1)



            fig.tight_layout()
            plt.subplots_adjust(top=0.84)

            return SimplePendulumFigAxes(
                fig=fig,
                ax1=ax_top_left,   # pendulum axis in single mode
                ax2=ax_top_right,  # phase axis in single mode
                ax3=ax_bottom,     # angle-vs-time axis in single mode
                time_text=time_text,
            )

        l2 = simulation_results_compare.l
        algorithm2 = simulation_results_compare.algorithm
        period2 = simulation_results_compare.period

        if angle_unit == "deg":
            theta2, omega2 = simulation_results_compare.theta_deg, simulation_results_compare.omega_deg
        else:
            theta2, omega2 = simulation_results_compare.theta_rad, simulation_results_compare.omega_rad

        suptitle = self._suptitle_compare(
            l1=float(l),
            theta10=float(theta[0]),
            omega10=float(omega[0]),
            period1=period,
            l2=float(l2),
            theta20=float(theta2[0]),
            omega20=float(omega2[0]),
            period2=period2,
            angle_unit=angle_unit,
        )
        fig.suptitle(suptitle, color="w")

        self._format_phase_axes(
            ax_top_left,
            angle_unit=angle_unit,
            title=f"Pendulum 1 ({algorithm})",
        )
        self._format_phase_axes(
            ax_top_right,
            angle_unit=angle_unit,
            title=f"Pendulum 2 ({algorithm2})",
        )
        self._format_angle_time_axes(
            ax_bottom,
            angle_unit=angle_unit,
            title="Angle vs. Time",
        )

        time_text = self._make_time_text(ax_bottom)

        x0, x1 = self._padded_limits(float(np.min(theta)), float(np.max(theta)), pad_frac=0.08)
        y0, y1 = self._padded_limits(float(np.min(omega)), float(np.max(omega)), pad_frac=0.08)
        mn, mx = min(x0, y0), max(x1, y1)
        ax_top_left.set_xlim(mn, mx)
        ax_top_left.set_ylim(mn, mx)
        ax_top_left.set_aspect("equal", adjustable="box")

        x0, x1 = self._padded_limits(float(np.min(theta2)), float(np.max(theta2)), pad_frac=0.08)
        y0, y1 = self._padded_limits(float(np.min(omega2)), float(np.max(omega2)), pad_frac=0.08)
        mn, mx = min(x0, y0), max(x1, y1)
        ax_top_right.set_xlim(mn, mx)
        ax_top_right.set_ylim(mn, mx)
        ax_top_right.set_aspect("equal", adjustable="box")

        x0, x1 = self._padded_limits(float(np.min(time)), float(np.max(time)), pad_frac=0.03)
        y0, y1 = self._padded_limits(
            float(min(np.min(theta), np.min(theta2))),
            float(max(np.max(theta), np.max(theta2))),
            pad_frac=0.03,
        )
        ax_bottom.set_xlim(x0, x1)
        ax_bottom.set_ylim(y0, y1)

        fig.tight_layout()
        plt.subplots_adjust(top=0.84)

        return SimplePendulumFigAxes(
            fig=fig,
            ax1=ax_top_left,
            ax2=ax_top_right,
            ax3=ax_bottom,
            time_text=time_text,
        )

    def animate_single_simulation(
        self,
        simulation_results: SimplePendulumResult,
        angle_unit: AngleUnit = "deg",
        color: str = "#65C6E3",
        skip_frames: int = 5,
        fps: int = 40,
        out_path: str | Path | None = None,
        blit: bool = False,
    ) -> None:
        l = simulation_results.l
        if l <= 0:
            raise ValueError("l must be > 0")
        if skip_frames <= 0:
            raise ValueError("skip_frames must be > 0")
        if fps <= 0:
            raise ValueError("fps must be > 0")

        if out_path is None:
            out_path = "./simple_pendulum.mp4"
        out_path = Path(out_path).expanduser().resolve()
        if out_path.suffix.lower() != ".mp4":
            raise ValueError(f"out_path must have '.mp4' suffix, got {out_path.suffix!r}")

        time = simulation_results.time
        theta_rad = simulation_results.theta_rad

        if angle_unit == "deg":
            theta_plot = simulation_results.theta_deg
            omega_plot = simulation_results.omega_deg
        elif angle_unit == "rad":
            theta_plot = simulation_results.theta_rad
            omega_plot = simulation_results.omega_rad
        else:
            raise ValueError(f"angle_unit must be one of ['deg', 'rad'], got {angle_unit!r}")

        fig_axes = self._get_figs_and_axes(simulation_results, None, angle_unit)

        ax_pendulum = fig_axes.ax1
        ax_phase = fig_axes.ax2
        ax_angle = fig_axes.ax3
        time_text = fig_axes.time_text

        frame_idx = np.arange(0, len(time), skip_frames, dtype=int)

        phase_curve, = ax_phase.plot([], [], color="coral", lw=2.0)
        phase_dot, = ax_phase.plot(
            [], [], 
            linestyle="None", 
            color=color,
            markersize=12,
            marker="o", 
            markeredgecolor="#bbbbbb",
            markeredgewidth=1.3, 
        )

        angle_curve, = ax_angle.plot([], [], color=color, lw=2.3)

        rod, bob = self._make_pendulum_axes(ax_pendulum, color)

        def init() -> tuple[Line2D, Line2D, Line2D, Line2D, Line2D, Text]:
            phase_curve.set_data([], [])
            phase_dot.set_data([], [])
            angle_curve.set_data([], [])
            rod.set_data([], [])
            bob.set_data([], [])
            time_text.set_text("")
            return phase_curve, phase_dot, angle_curve, rod, bob, time_text

        def update(frame_number: int) -> tuple[Line2D, Line2D, Line2D, Line2D, Line2D, Text]:
            i = int(frame_idx[frame_number])

            phase_curve.set_data(theta_plot[: i + 1], omega_plot[: i + 1])
            phase_dot.set_data([theta_plot[i]], [omega_plot[i]])

            angle_curve.set_data(time[: i + 1], theta_plot[: i + 1])

            th = float(theta_rad[i])
            x = float(np.sin(th))
            y = -float(np.cos(th))
            rod.set_data([0.0, x], [0.0, y])
            bob.set_data([x], [y])

            time_text.set_text(f"t = {time[i]:04.1f}s")
            return phase_curve, phase_dot, angle_curve, rod, bob, time_text

        ani = animation.FuncAnimation(
            fig_axes.fig,
            update,
            frames=len(frame_idx),
            init_func=init,
            interval=int(1000 / fps),
            blit=blit,
        )

        writer = animation.FFMpegWriter(fps=fps)
        ani.save(out_path, writer=writer)
        plt.close(fig_axes.fig)

    def animate_comparison_simulations(
        self,
        simulation_results1: SimplePendulumResult,
        simulation_results2: SimplePendulumResult,
        angle_unit: AngleUnit = "deg",
        color1: str = "#65C6E3",
        color2: str = "#E225DC",
        skip_frames: int = 5,
        fps: int = 40,
        out_path: str | Path | None = None,
        blit: bool = False,
    ) -> None:
        l1, l2 = simulation_results1.l, simulation_results2.l
        if l1 <= 0 or l2 <= 0:
            raise ValueError("l1 and l2 must both be > 0")

        if skip_frames <= 0:
            raise ValueError("skip_frames must be > 0")
        if fps <= 0:
            raise ValueError("fps must be > 0")

        if out_path is None:
            out_path = "./simple_pendulum.mp4"
        out_path = Path(out_path).expanduser().resolve()
        if out_path.suffix.lower() != ".mp4":
            raise ValueError(f"out_path must have '.mp4' suffix, got {out_path.suffix!r}")

        time = simulation_results1.time
        theta1_rad, theta2_rad = simulation_results1.theta_rad, simulation_results2.theta_rad
        alg1, alg2 = simulation_results1.algorithm, simulation_results2.algorithm

        if angle_unit == "deg":
            theta1_plot, theta2_plot = simulation_results1.theta_deg, simulation_results2.theta_deg
            omega1_plot, omega2_plot = simulation_results1.omega_deg, simulation_results2.omega_deg
        elif angle_unit == "rad":
            theta1_plot, theta2_plot = simulation_results1.theta_rad, simulation_results2.theta_rad
            omega1_plot, omega2_plot = simulation_results1.omega_rad, simulation_results2.omega_rad
        else:
            raise ValueError(f"angle_unit must be one of ['deg', 'rad'], got {angle_unit!r}")

        fig_axes = self._get_figs_and_axes(simulation_results1, simulation_results2, angle_unit)

        ax1 = fig_axes.ax1
        ax2 = fig_axes.ax2
        ax3 = fig_axes.ax3
        time_text = fig_axes.time_text

        frame_idx = np.arange(0, len(time), skip_frames, dtype=int)

        phase_curve1, = ax1.plot([], [], color="coral", lw=2.0)
        phase_dot1, = ax1.plot(
            [], [], 
            linestyle="None", 
            color=color1,
            markersize=12,
            marker="o", 
            markeredgecolor="#bbbbbb",
            markeredgewidth=1.3, 
        )

        phase_curve2, = ax2.plot([], [], color="coral", lw=2.0)
        phase_dot2, = ax2.plot(
            [], [], 
            linestyle="None", 
            color=color2,
            markersize=12,
            marker="o", 
            markeredgecolor="#bbbbbb",
            markeredgewidth=1.3, 
        )

        angle_curve1, = ax3.plot([], [], color=color1, lw=2.3, label=alg1)
        angle_curve2, = ax3.plot([], [], color=color2, lw=2.3, label=alg2)

        legend = ax3.legend(
            loc="upper right",
            frameon=True,
            facecolor="black",
            edgecolor="white",
            fontsize=14,
        )
        for text in legend.get_texts():
            text.set_color("white")

        _, rod1, bob1 = self._make_pendulum_inset(ax1, color1)
        _, rod2, bob2 = self._make_pendulum_inset(ax2, color2)

        def init() -> tuple[
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Text,
        ]:
            phase_curve1.set_data([], [])
            phase_dot1.set_data([], [])
            phase_curve2.set_data([], [])
            phase_dot2.set_data([], [])
            angle_curve1.set_data([], [])
            angle_curve2.set_data([], [])
            rod1.set_data([], [])
            bob1.set_data([], [])
            rod2.set_data([], [])
            bob2.set_data([], [])
            time_text.set_text("")
            return (
                phase_curve1,
                phase_dot1,
                phase_curve2,
                phase_dot2,
                angle_curve1,
                angle_curve2,
                rod1,
                bob1,
                rod2,
                bob2,
                time_text,
            )

        def update(frame_number: int) -> tuple[
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Line2D,
            Text,
        ]:
            i = int(frame_idx[frame_number])

            phase_curve1.set_data(theta1_plot[: i + 1], omega1_plot[: i + 1])
            phase_dot1.set_data([theta1_plot[i]], [omega1_plot[i]])

            phase_curve2.set_data(theta2_plot[: i + 1], omega2_plot[: i + 1])
            phase_dot2.set_data([theta2_plot[i]], [omega2_plot[i]])

            angle_curve1.set_data(time[: i + 1], theta1_plot[: i + 1])
            angle_curve2.set_data(time[: i + 1], theta2_plot[: i + 1])

            th1 = float(theta1_rad[i])
            x1 = float(np.sin(th1))
            y1 = -float(np.cos(th1))
            rod1.set_data([0.0, x1], [0.0, y1])
            bob1.set_data([x1], [y1])

            th2 = float(theta2_rad[i])
            x2 = float(np.sin(th2))
            y2 = -float(np.cos(th2))
            rod2.set_data([0.0, x2], [0.0, y2])
            bob2.set_data([x2], [y2])

            time_text.set_text(f"t = {time[i]:04.1f}s")

            return (
                phase_curve1,
                phase_dot1,
                phase_curve2,
                phase_dot2,
                angle_curve1,
                angle_curve2,
                rod1,
                bob1,
                rod2,
                bob2,
                time_text,
            )

        ani = animation.FuncAnimation(
            fig_axes.fig,
            update,
            frames=len(frame_idx),
            init_func=init,
            interval=int(1000 / fps),
            blit=blit,
        )

        writer = animation.FFMpegWriter(fps=fps)
        ani.save(out_path, writer=writer)
        plt.close(fig_axes.fig)
