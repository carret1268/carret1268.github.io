from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from typing import Any, Literal, TypeAlias

import numpy as np
import numpy.typing as npt


FloatArray: TypeAlias = npt.NDArray[np.floating]
AlgorithmFunc: TypeAlias = Callable[
    [float, float, FloatArray],
    tuple[FloatArray, FloatArray],
]

Algorithm: TypeAlias = Literal["saa", "rk4"]


@dataclass(frozen=True)
class SimplePendulumResult:
    """
    Immutable container for simple pendulum simulation results.

    All state variables are stored in SI units (radians and seconds) to
    avoid ambiguity during numerical computation. Convenience properties
    are provided for accessing angular quantities in degrees.
    """

    l: float
    algorithm: Algorithm
    time: FloatArray
    theta_rad: FloatArray
    omega_rad: FloatArray
    period: float | None

    @property
    def theta_deg(self) -> FloatArray:
        return self.theta_rad * (180.0 / np.pi)

    @property
    def omega_deg(self) -> FloatArray:
        return self.omega_rad * (180.0 / np.pi)


class SimplePendulum:
    """
    Simulator for an ideal simple pendulum.

    This class models a point mass suspended by a massless, rigid rod in a
    uniform gravitational field. Multiple solution strategies can be
    selected, including analytical approximations and numerical
    integrators.
    """

    def __init__(
        self,
        l: float = 1.0,
        g: float = 9.81,
        array_dtype: np.dtype[Any] = np.dtype(np.float64),
    ):
        """
        Initialize a simple pendulum model.

        Parameters
        ----------
        l : float, optional
            Length of the pendulum (meters). Must be positive.
        g : float, optional
            Gravitational acceleration (m/s^2). Must be positive.
        array_dtype : numpy.dtype, optional
            NumPy dtype used for internally generated arrays.

        Raises
        ------
        ValueError
            If `l` or `g` is non-positive.
        """
        if l <= 0:
            raise ValueError("l must be > 0")
        if g <= 0:
            raise ValueError("g must be > 0")

        self.l = float(l)
        self.g = float(g)
        self.omega_n = float(np.sqrt(self.g / self.l))
        self.array_dtype: np.dtype[Any] = array_dtype

        self.algorithm_map: dict[str, AlgorithmFunc] = {
            "saa": self._sol_small_angle_approx,
            "rk4": self._sol_rk4,  # stub for now
        }

    def _calculate_period(
        self,
        omega: FloatArray,
        time: FloatArray,
    ) -> float | None:
        if len(omega) != len(time):
            raise ValueError("omega and time must have the same length")
        if len(omega) < 3:
            return None

        eps = 1e-12

        def _sign(x: float) -> int:
            if x > eps:
                return 1
            if x < -eps:
                return -1
            return 0

        crossings: list[tuple[float, int]] = []  # (t_cross, direction), direction: -1 (+ to -), +1 (- to +)

        s_prev = _sign(float(omega[0]))
        for i in range(1, len(omega)):
            s = _sign(float(omega[i]))

            if s_prev == 0:
                s_prev = s
                continue

            if s == 0:
                continue

            if s != s_prev:
                # linear interpolate crossing time between i-1 and i
                w0 = float(omega[i - 1])
                w1 = float(omega[i])
                t0 = float(time[i - 1])
                t1 = float(time[i])
                if w1 == w0:
                    t_cross = t0
                else:
                    frac = -w0 / (w1 - w0)
                    t_cross = t0 + frac * (t1 - t0)

                direction = 1 if (s_prev < 0 and s > 0) else -1
                crossings.append((t_cross, direction))
                s_prev = s

                # stop once we have two crossings of the same direction
                if len(crossings) >= 2 and crossings[-1][1] == crossings[0][1]:
                    return crossings[-1][0] - crossings[0][0]
            else:
                s_prev = s

        return None

    def _sol_small_angle_approx(
        self,
        theta0_rad: float,
        omega0_rad_s: float,
        t: FloatArray,
    ) -> tuple[FloatArray, FloatArray]:
        """
        Analytical solution using the small-angle approximation.

        Assumes sin(theta) ≈ theta, reducing the equation of motion to that
        of a simple harmonic oscillator with natural frequency
        sqrt(g / l).

        Parameters
        ----------
        theta0_rad : float
            Initial angular displacement in radians.
        omega0_rad_s : float
            Initial angular velocity in radians per second.
        t : FloatArray
            Time array at which to evaluate the solution.

        Returns
        -------
        tuple[FloatArray, FloatArray]
            Tuple of angular displacement and angular velocity arrays,
            both evaluated at times `t`.
        """
        w = self.omega_n

        theta0 = self.array_dtype.type(theta0_rad)
        omega0 = self.array_dtype.type(omega0_rad_s)

        theta = theta0 * np.cos(w * t) + (omega0 / w) * np.sin(w * t)

        omega = -theta0 * w * np.sin(w * t) + omega0 * np.cos(w * t)

        return theta, omega

    def _sol_rk4(
        self,
        theta0_rad: float,
        omega0_rad_s: float,
        t: FloatArray,
    ) -> tuple[FloatArray, FloatArray]:
        """
        Numerical solution using a fourth-order Runge-Kutta integrator.

        This integrates the full nonlinear simple pendulum equation of motion:

            dtheta/dt = omega
            domega/dt = -(g/l) * sin(theta)

        Parameters
        ----------
        theta0_rad : float
            Initial angular displacement in radians.
        omega0_rad_s : float
            Initial angular velocity in radians per second.
        t : FloatArray
            Monotonically increasing time array.

        Returns
        -------
        tuple[FloatArray, FloatArray]
            Tuple of (theta, omega) arrays evaluated at times `t`.

        Raises
        ------
        ValueError
            If `t` has fewer than 2 points or is not evenly spaced.
        """
        if t.size < 2:
            raise ValueError(
                "t must contain at least 2 time points for RK4 integration"
            )

        # infer dt from the time array
        dt = float(t[1] - t[0])
        if dt <= 0:
            raise ValueError("t must be strictly increasing")

        dtype = self.array_dtype
        theta = np.empty(t.shape, dtype=dtype)
        omega = np.empty(t.shape, dtype=dtype)

        theta[0] = dtype.type(theta0_rad)
        omega[0] = dtype.type(omega0_rad_s)

        c = dtype.type(self.g / self.l)  # constant factor (g/l)

        def f_theta(_th: float, _om: float) -> float:
            return _om

        def f_omega(_th: float, _om: float) -> float:
            # _om is unused here, but kept for symmetry/readability
            return -c * np.sin(_th)

        dt_t = dtype.type(dt)
        HALF = dtype.type(0.5)
        SIXTH = dtype.type(1.0 / 6.0)

        for i in range(t.size - 1):
            th = theta[i]
            om = omega[i]

            k1_th = f_theta(th, om)
            k1_om = f_omega(th, om)

            th2 = th + HALF * dt_t * k1_th
            om2 = om + HALF * dt_t * k1_om
            k2_th = f_theta(th2, om2)
            k2_om = f_omega(th2, om2)

            th3 = th + HALF * dt_t * k2_th
            om3 = om + HALF * dt_t * k2_om
            k3_th = f_theta(th3, om3)
            k3_om = f_omega(th3, om3)

            th4 = th + dt_t * k3_th
            om4 = om + dt_t * k3_om
            k4_th = f_theta(th4, om4)
            k4_om = f_omega(th4, om4)

            theta[i + 1] = th + dt_t * SIXTH * (
                k1_th + 2.0 * k2_th + 2.0 * k3_th + k4_th
            )
            omega[i + 1] = om + dt_t * SIXTH * (
                k1_om + 2.0 * k2_om + 2.0 * k3_om + k4_om
            )

        return theta, omega

    def simulate(
        self,
        theta0_deg: float,
        omega0_deg_s: float,
        t_f: float,
        dt: float,
        algorithm: Algorithm = "saa",
    ) -> SimplePendulumResult:
        """
        Simulate the motion of a simple pendulum.

        Parameters
        ----------
        theta0_deg : float
            Initial angular displacement in degrees.
        omega0_deg_s : float
            Initial angular velocity in degrees per second.
        t_f : float
            Final simulation time in seconds.
        dt : float
            Time step size in seconds.
        algorithm : Algorithm, optional
            Algorithm to use (e.g., ``"saa"`` or ``"rk4"``).
            "saa" = "small angle approximation".

        Returns
        -------
        SimplePendulumResult
            Immutable container holding the simulation time series.

        Raises
        ------
        ValueError
            If `t_f` or `dt` is non-positive, or if an unknown solution
            identifier is provided.
        """
        if t_f <= 0:
            raise ValueError("t_f must be > 0")
        if dt <= 0:
            raise ValueError("dt must be > 0")

        try:
            solver = self.algorithm_map[algorithm]
        except KeyError as e:
            raise ValueError(f"Unknown algorithm: {algorithm!r}") from e

        time = np.arange(0.0, t_f, dt, dtype=self.array_dtype)

        theta0_rad = float(np.deg2rad(theta0_deg))
        omega0_rad_s = float(np.deg2rad(omega0_deg_s))

        theta, omega = solver(theta0_rad, omega0_rad_s, time)

        return SimplePendulumResult(
            l=self.l,
            algorithm=algorithm,
            time=time,
            theta_rad=theta.astype(self.array_dtype, copy=False),
            omega_rad=omega.astype(self.array_dtype, copy=False),
            period=self._calculate_period(omega, time)
        )

    
