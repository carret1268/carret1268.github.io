---
layout: post
title: "Simple Pendulum"
categories: [physics, classical-mechanics, numerical-methods]
tags: [pendulum, simple-harmonic-motion, differential-equations, rk4, phase-space, energy-conservation, elliptic-integrals, simulation, python, matplotlib]
---

## Introduction

In this post I’m going to build up the simple pendulum in a few layers:

1. Derive the exact (nonlinear) equation of motion from Newton’s 2nd law.
2. Use the small angle approximation to get a closed-form solution and the familiar period formula.
3. Simulate the motion (first with the small angle solution, then with RK4 for the exact equation) and compare what changes as amplitude increases.
4. Explain why the small-angle phase diagram is an ellipse (two viewpoints: parametric geometry and conservation of energy).
5. Derive the exact nonlinear period in terms of an elliptic integral, show how to compute it efficiently (AGM), and connect that to the separatrix.

---

## $\oint$ Equations of Motion

### $\unicode{x222F}$ Deriving the Equations of Motion

Consider a simple pendulum: a "bob" of mass $m$ connected to a perfectly rigid and massless rod of length $l$ whose other end is attached to a pivot point. There is no friction or air resistance. We displace the bob so the rod makes an angle $\theta$ relative to the vertical, then let it go. We want to predict how $\theta(t)$ evolves.

We begin by drawing a free body diagram (using polar coordinates).

![Free body diagram of simple pendulum]({{ '/assets/media/physics/pendulums/simple-pendulum/figures/simplePendFBD.svg' | relative_url }})

To derive the equations of motion, we write Newton's second law in the radial and tangential directions. Because the rod is rigid, the radius is constant: $r=l$, so $\dot{r}=0$ and $\ddot{r}=0$. However, the bob *does* have radial acceleration due to its changing direction (centripetal acceleration), so the radial equation is

$$
  \sum F_r = m(\ddot{r} - r\dot{\theta}^2) = -ml\dot{\theta}^2 = -T_r + mg\cos\theta
$$

\begin{equation}\label{eq:tension}
  \Rightarrow T_r = m(g\cos\theta + l\dot{\theta}^2),
\end{equation}

which is the tension in our rod - we see that the tension due to gravity is a maximum when the pendulum hangs straight down ($\theta = 0$), and the quicker it rotates, the more tension there is.

Now we consider the tangential direction. The tangential acceleration is $a_\theta = l\ddot{\theta}$. The tangential component of gravity is $-mg\sin\theta$ (restoring toward $\theta=0$). So Newton’s second law gives

$$
  \sum F_\theta = ma_\theta = ml\ddot{\theta} = -mg\sin\theta
$$

$$
  \Rightarrow \frac{d^2\theta}{dt^2}=-\frac{g}{l}\sin\theta = -\omega_n^2\sin\theta
$$

\begin{equation}
  \label{eq:simple-pend}
  \Rightarrow \frac{d^2\theta}{dt^2} = -\omega_n^2\sin\theta
\end{equation}

Where <em>g=9.81m/s<sup>2</sup></em> is acceleration due to gravity, and we have defined the natural angular frequency

$$
  \omega_n = \sqrt{\frac{g}{l}}
$$

I pause here to note how fascinating it is that the angular frequency of such a pendulum is strictly a function of its length - of course this is a consequence of the equally surprising result that all object fall at the same rate regardless of their mass.

Eq. \eqref{eq:simple-pend} is the exact second-order nonlinear differential equation of motion. There *is* an analytical solution, but it involves elliptic integrals (and is not a simple elementary function). We’ll take two practical routes:

1. Use the small angle approximation to get an analytic solution.
2. Solve the exact equation numerically and compare.

For the uninitiated, the small angle approximation says that for small angles (in radians),

$$
  \sin\theta\approx\theta
$$

Perhaps the following table can convince you of its validity:

| **$\theta$ (degrees)** | **$\theta$ (radians)** | **$\sin\theta$** | **Percent Error (%)** |
|----|----|----|----|
| 0   | 0.0000 | 0.0000 | 0.00 |
| 1   | 0.0175 | 0.0175 | 0.00 |
| 5   | 0.0873 | 0.0872 | 0.13 |
| 10  | 0.1745 | 0.1736 | 0.52 |
| 20  | 0.3491 | 0.3420 | 2.03 |
| 30  | 0.5236 | 0.5000 | 4.51 |
| 40  | 0.6981 | 0.6428 | 7.92 |
| 50  | 0.8727 | 0.7660 | 12.23 |
| 90  | 1.5708 | 1.0000 | 36.34 |
| 120 | 2.0944 | 0.8660 | 58.66 |

It is clear that the approximation is quite good for small angles, and gets worse as the angle grows. This means that simplifying our differential equation to

\begin{equation}\label{eq:homogeneous}
  \frac{d^2\theta}{dt^2}=-\omega_n^2\theta
\end{equation}

will act as a good approximation as long as we keep the initial angular displacement fairly small. With this simplification, we can solve the motion analytically with standard methods.

---

## $\oint$ Small Angle Approximation Solution

### $\unicode{x222F}$ Solving the Equations of Motion

First, rewrite \eqref{eq:homogeneous} as

$$
  \frac{d^2\theta}{dt^2} + \omega_n^2\theta = 0
$$

Which is the famous simple harmonic oscillator, whose characteristic equation is

$$
  r^2 + \omega_n^2 = 0 \qquad \Rightarrow \qquad r = \pm i\omega_n
$$

So the general solution is

\begin{equation}
  \label{eq:general_soln}
  \theta(t) = A\cos(\omega_n t) + B\sin(\omega_n t)
\end{equation}

To solve for $A$ and $B$ we apply initial conditions. Let

$$
  \theta(t=0) = \theta_0 \qquad \Rightarrow \qquad \dot{\theta}(t=0) = \omega_0
$$

where $\theta_0$ is the initial angular displacement and $\omega_0$ is the initial angular velocity. Using this with Eq. \eqref{eq:general_soln} yields

$$
  \theta(t=0) = \theta_0 = A
$$

Now for $B$, take the derivative of Eq. \eqref{eq:general_soln} with respect to time:

$$
  \dot{\theta} = -A\omega_n\sin(\omega_n t) + B\omega_n\cos(\omega_n t)
$$

Using $\dot{\theta}(t=0)=\omega_0$ gives

$$
  \dot{\theta}(t=0) = \omega_0 = B\omega_n
$$

$$
  \Rightarrow B = \frac{\omega_0}{\omega_n}
$$

Hence, our solution is

\begin{equation}\label{eq:final_general_soln}
  \theta(t) = \theta_0\cos(\omega_n t) + \frac{\omega_0}{\omega_n}\sin(\omega_n t)
\end{equation}

Notice how the right-most term is the contribution of a non-zero initial angular velocity.

And for the sake of completeness, we will write its derivative,

\begin{equation}\label{eq:derivative_of_general_soln}
  \dot{\theta} = \theta_0 \omega_n \sin(\omega_n t) + \omega_0\cos(\omega_n t)
\end{equation}

---

## $\oint$ Determining the Small Angle Period $T_0$

From the small-angle model, we can immediately estimate the period. Recall the relationship between frequency $f$ and angular frequency:

$$
f = \frac{1}{2\pi}\omega_n = \frac{1}{2\pi}\sqrt{\frac{g}{l}}
$$

$$
T_0 = \frac{1}{f} = 2\pi\sqrt{\frac{l}{g}}
$$

So if we let the length of our pendulum be $l=9.81$m, we should expect $T_0 = 2\pi\approx 6.28$ seconds (since then $\omega_n=1$). This is not a function of the mass, but strictly the length and gravity. Again: not intuitive at first, but very real.

---

## $\oint$ Simulating with Small Angle Approximation

Now we’ll simulate a pendulum released from some initial conditions, and see how well the small-angle prediction holds.

We will start by importing the tools we need, building some custom `TypeAlias`, and constructing a custom return type for our simulations.

[View 'assets/code/physics/pendulums/simple-pendulum/simple_pendulum.py' on GitHub](https://github.com/carret1268/carret1268.github.io/blob/main/assets/code/physics/pendulums/simple-pendulum/simple_pendulum.py)

```python
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
```

Our script begins by enabling postponed evaluation of annotations (`from __future__ import annotations`), which allows type hints to reference classes and aliases without runtime ordering issues. Several type aliases are then defined to clarify intent:

- `FloatArray` standardizes the representation of NumPy floating-point arrays using `numpy.typing`.
- `AlgorithmFunc` describes the callable signature expected of solution backends: given initial conditions and a time array, the function returns angle and angular velocity arrays.
- `Algorithm` constrains solver selection to a small, explicit set of valid identifiers (e.g., "saa" = small angle approximation, "rk4" = 4th order Runge-Kutta), preventing accidental misuse.

The `SimplePendulumResult` dataclass provides an immutable container for simulation output: time samples, angular displacement $\theta(t)$, and angular velocity $\omega(t)$ stored in radians, plus convenience properties for degrees.

Now we define the `SimplePendulum` class within the same file:

```python
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
        of a simple harmonic oscillator with natural frequency sqrt(g / l).
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
        raise NotImplementedError("RK4 solver not implemented yet")

    def simulate(
        self,
        theta0_deg: float,
        omega0_deg_s: float,
        t_f: float,
        dt: float,
        algorithm: Algorithm = "saa",
    ) -> SimplePendulumResult:
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
```

At construction time, `SimplePendulum` validates and stores the pendulum length `l` and gravitational acceleration `g`, computes the natural frequency $\omega_n=\sqrt{g/l}$, and records an internal array dtype that controls the precision of all generated NumPy arrays. The small dispatch table (`algorithm_map`) makes it straightforward to add new solution strategies.

Beyond basic simulation, `_calculate_period` estimates the oscillation period from successive zero-crossings of angular velocity. The crossing time is refined using linear interpolation between samples. If the simulation window is too short to observe a full cycle, it returns `None` rather than guessing.

---

## $\oint$ Animations and What to Look For

### $\unicode{x222F}$ Animating

Creating animations with `matplotlib` can be quite a challenge for the uninitiated, and deserves more space and focus than I have to give it here. I’ve chosen to omit the animation code from this post for the sake of brevity, but you may view it here:

[View 'assets/code/physics/pendulums/simple-pendulum/simple_pendulum_animator.py' on GitHub](https://github.com/carret1268/carret1268.github.io/blob/main/assets/code/physics/pendulums/simple-pendulum/simple_pendulum_animator.py)

### $\unicode{x222F}$ Script for Simulating and Animating

This script generates all figures in this post:

[View 'assets/code/physics/pendulums/simple-pendulum/simulate_and_animate.py' on GitHub](https://github.com/carret1268/carret1268.github.io/blob/main/assets/code/physics/pendulums/simple-pendulum/simulate_and_animate.py)

### $\unicode{x222F}$ Animation Output for Small Angle Approximation

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/simple_pend-9_81-10-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

Firstly, note that the period estimate from `SimplePendulum._calculate_period` is $T\approx 6.28$s, matching the small-angle prediction for $l=9.81$m.

Secondly, look at the phase diagram in the top right: it’s a closed curve (and in this special case, a circle). For the small-angle solution with $\omega_0=0$ we have

$$
x(t) = \theta(t) = \theta_0\cos(\omega_n t) \qquad y(t) = \dot{\theta}(t) = -\theta_0\omega_n\sin(\omega_n t)
$$

This is exactly the parametric form of an ellipse. The semi-axes are

$$
a = \theta_0 \qquad b = \theta_0\omega_n,
$$

which means the eccentricity is defined as

\begin{equation}\label{eq:eccentricity}
  e = \sqrt{1 - \frac{a^2}{b^2}} = \sqrt{1 - \frac{1}{\omega_n^2}}
\end{equation}

So if we choose $l$ such that $\omega_n=1$, then $a=b$, $\,e=0$ and the ellipse becomes a circle. That’s what you’re seeing in the $l=9.81$m animation.

Now see what happens when we set $l = 2$ (so $\omega_n>1$ and $b>a$):

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/simple_pend-2-10-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

Now the phase diagram is an ellipse that’s “taller” in $\omega$ because the maximum angular speed scales with $\omega_n$. Also notice that the period decreases due to the shorter rod.

Before moving on to the exact equation, let’s see how changing the initial angular velocity changes the curves. We should see a phase shift in each plot, but the period should not change. With $\omega_0\neq 0$ the parametric equations become

$$
x(t) = \theta(t) = \theta_0\cos(\omega_n t) + \frac{\omega_0}{\omega_n}\sin(\omega_n t)
$$

$$
y(t) = \dot{\theta}(t) = -\theta_0\omega_n\sin(\omega_n t) + \omega_0\cos(\omega_n t)
$$

If you square both equations and add them appropriately, you again get an ellipse:

$$
  \frac{x^2}{\theta_0^2 + \left(\frac{\omega_0}{\omega_n}\right)^2} + \frac{y^2}{\omega_n^2\left(\theta_0^2 + \left(\frac{\omega_0}{\omega_n}\right)^2\right)} = 1
$$

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/simple_pend-9_81-20-20.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

---

## $\oint$ Why the phase diagram is an Ellipse (Small Angle)

### $\unicode{x222f}$ Why Ellipse?

At this point, if you asked why the small-angle phase diagram is elliptical, I could say “because sine and cosine,” and that would be true - it is, after all, just the parametric equation of an ellipse. I could also say something broader, like “because of Newton’s laws, the conservative nature of the forces acting on the pendulum, and the geometry of the problem,” which is also correct. To me, neither answer feels especially insightful. So instead, I’ll approach this from a different angle - through conservation of energy - which I hope you’ll find a bit more intellectually satisfying.

For an ideal pendulum, the total mechanical energy is kinetic plus potential. With a point mass at radius $l$, the moment of inertia about the pivot is

$$
  I = ml^2
$$

and the height relative to the lowest point is $h = l(1-\cos\theta)$. So

$$
  E = T + U = \frac{1}{2}I\omega^2 + mgh
$$

\begin{equation}\label{eq:consv-energy}
  E = \frac{1}{2}ml^2\omega^2 + mgl(1 - \cos\theta)
\end{equation}

Energy is conserved, so the phase-space trajectory is just the set of points $(\theta, \omega)$ that satisfy Eq. \eqref{eq:consv-energy} for some constant $E$.

For small angular displacements, when $\vert\theta\vert \ll 1$, we can approximate

$$
  \cos\theta \approx 1 - \frac{\theta^2}{2}
$$

which makes the potential energy approximately quadratic:

$$
  U = mgl(1-\cos\theta) \approx \frac{1}{2}mgl\theta^2
$$

Plugging this back into Eq. \eqref{eq:consv-energy} gives

$$
  E \approx \frac{1}{2}ml^2\omega^2 + \frac{1}{2}mgl\theta^2
$$

Divide by $\frac{1}{2}ml^2$:

\begin{equation}\label{eq:conv-energy-result}
  \varepsilon = \omega^2 + \frac{g}{l}\theta^2
\end{equation}

where $\varepsilon = \frac{2E}{ml^2}$ is a constant. Rearranging makes the ellipse obvious:

$$
  1 = \frac{1}{\varepsilon}\omega^2 + \frac{g}{l\varepsilon}\theta^2
$$

$$
  \Rightarrow \frac{\omega^2}{a^2} + \frac{\theta^2}{b^2} = 1
$$

Where

$$
  a = \sqrt{\varepsilon} \qquad b = \sqrt{\frac{\varepsilon l}{g}}
$$

So for $l > g$, $\,0 < \vert\omega_n\vert < 1$ the eccentricity is

$$
  e = \sqrt{1 - \frac{g}{l}} = \sqrt{1 - \omega_n^2}
$$

for $l \leq g$, $\,\vert\omega_n\vert \geq 1$ the eccentricity is exactly what we found in Eq. \eqref{eq:eccentricity}

$$
  e = \sqrt{1 - \frac{l}{g}} = \sqrt{1 - \frac{1}{\omega_n^2}}
$$

So the small-angle phase diagram is elliptical because (under the small-angle approximation) the energy contours are quadratic in $(\theta,\omega)$. To put it plainly, the ellipses fall are a result of the energy of the system being conserved.

Just for fun, we can retrieve our small-angle equation of motion by differentiating Eq. \eqref{eq:conv-energy-result} with respect to time. Since $\varepsilon$ is constant:

$$
  0 = \frac{d}{dt}\left(\omega^2 + \frac{g}{l}\theta^2\right)
$$

Applying the chain rule:

$$
  2\omega\frac{d\omega}{dt} + 2\frac{g}{l}\theta\frac{d\theta}{dt} = 0
$$

Using $\omega=\dot{\theta}$ and $\dot{\omega}=\ddot{\theta}$ yields

$$
  \Rightarrow \frac{d^2\theta}{dt^2} + \omega_n^2\theta = 0
$$

---

## $\oint$ Solving the Exact Equation Numerically

We’ve made great progress, but I want to see (visually) when the small angle approximation fails, and how it degrades at larger amplitudes. Start again with Eq. \eqref{eq:simple-pend}:

$$
  \frac{d^2\theta}{dt^2} = -\omega_n^2\sin\theta
$$

To solve this second-order ODE numerically, we convert it into two coupled first-order ODEs. Let

\begin{equation}\label{eq:dtheta}
  \frac{d\theta}{dt} = \omega
\end{equation}

and therefore

\begin{equation}\label{eq:domega}
  \frac{d\omega}{dt} = -\omega_n^2\sin\theta
\end{equation}

Now we can solve $(\theta,\omega)$ using many techniques. For this post, we’ll use 4th order Runge-Kutta (RK4).

### $\unicode{x222F}$ RK4 Implementation

We update `_sol_rk4` as follows:

```python
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
```

Given an initial angle and angular velocity, the routine infers a fixed time step from the supplied time array, allocates output arrays, and then advances the state one step at a time. At each step, four intermediate slope evaluations ($k_1$ through $k_4$) are computed for both $\theta$ and $\omega$, then combined in the standard RK4 weighted average.

---

## $\oint$ Comparing Approximate and Exact Solutions

Let’s look at the same example with an initial angular displacement of 10 degrees.

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-10.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

They are almost identical. First, notice that the RK4 period is very slightly longer, e.g. $T_2 \approx 6.30$s compared to $T_1 \approx 6.28$s for the small angle approximation (numbers will vary slightly with timestep and how the period is estimated). Visually, this is hard to see until late in the simulation, where the angle vs. time curves start to drift out of phase.

Now we increase the initial displacement. Watch (1) how the period changes, and (2) how the RK4 phase diagram deforms away from an ellipse.

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-20.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-30.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-45.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-60.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-90.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-120.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-150.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/compare_saa_rk4-179.0-0.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

As the initial angular displacement increases, the exact solution clearly shows the anharmonic nature of the pendulum: the period is not strictly a function of length alone, but also depends on amplitude.

So at what point does the small angle approximation stop being “good”? I’ll still leave that for you to decide, but now you’ve got a way to quantify it (and visualize it).

---

## $\oint$ The Exact Period $T(\theta_0)$ of the Nonlinear Simple Pendulum

Let’s revisit Eq. \eqref{eq:consv-energy}:

$$
  E = \frac{1}{2}ml^2\dot{\theta}^2 + mgl(1 - \cos\theta)
$$

We define a scaled energy variable

$$
  \varepsilon = \frac{E}{mgl},
$$

so dividing by $mgl$ gives

\begin{equation}\label{eq:naturalish-energy}
  \varepsilon = \frac{l}{2g}\dot{\theta}^2 + (1 - \cos\theta).
\end{equation}

Since energy is conserved, evaluate $\varepsilon$ at a turning point $\theta = \theta_0$, where $\dot{\theta}=0$:

$$
  \varepsilon = 1 - \cos\theta_0.
$$

Substitute back into Eq. \eqref{eq:naturalish-energy} and solve for $\dot{\theta}$:

$$
  1 - \cos\theta_0
  =
  \frac{l}{2g}\dot{\theta}^2 + (1 - \cos\theta),
$$

which simplifies to

$$
  \frac{l}{2g}\dot{\theta}^2 = \cos\theta - \cos\theta_0.
$$

Therefore,

\begin{equation}\label{eq:dottheta}
  \dot{\theta} = \pm \sqrt{ \frac{2g}{l}(\cos\theta - \cos\theta_0) }.
\end{equation}

Although the motion is not simple harmonic, the oscillation is symmetric. So we compute the time required to move from $\theta = 0$ to $\theta = \theta_0$ and multiply by four:

$$
  \frac{T}{4}
  =
  \int_0^{\theta_0}
  \frac{d\theta}{\sqrt{ \frac{2g}{l}(\cos\theta - \cos\theta_0) }}.
$$

Factor constants:

$$
  \frac{T}{4}
  =
  \sqrt{\frac{l}{2g}}
  \int_0^{\theta_0}
  \frac{d\theta}{\sqrt{ \cos\theta - \cos\theta_0 }}.
$$

Using

$$
  \cos\theta - \cos\theta_0
  =
  2\left(\sin^2\frac{\theta_0}{2} - \sin^2\frac{\theta}{2}\right),
$$

and defining

$$
  k = \sin\frac{\theta_0}{2}
$$

we obtain

\begin{equation}\label{eq:quarter-period}
  \frac{T}{4} = \sqrt{\frac{l}{g}} \int_0^{\theta_0} \frac{d\theta}{2\sqrt{k^2 - \sin^2\frac{\theta}{2}}}.
\end{equation}

To write this in standard form, make the substitution

\begin{equation}\label{eq:sintheta-sinphi}
  \sin\frac{\theta}{2} = k\sin\phi.
\end{equation}

As $\theta \to 0$, $\phi \to 0$, and as $\theta \to \theta_0$, $\phi \to \pi/2$. Differentiating:

$$
  \frac{1}{2}\cos\frac{\theta}{2}\,d\theta = k\cos\phi\,d\phi,
$$

so

$$
  d\theta = \frac{2k\cos\phi}{\cos\frac{\theta}{2}}\,d\phi.
$$

Substitute into Eq. \eqref{eq:quarter-period}:

\begin{equation}\label{eq:quarter-period2}
  \frac{T}{4} = \sqrt{\frac{l}{g}} \int_0^{\pi/2} \frac{k\cos\phi}{\cos\frac{\theta}{2}\sqrt{k^2 - \sin^2\frac{\theta}{2}}} \,d\phi.
\end{equation}

Simplify using

\begin{equation}\label{eq:simp1}
  \cos\frac{\theta}{2} = \sqrt{1 - \sin^2\frac{\theta}{2}} = \sqrt{1 - k^2\sin^2\phi},
\end{equation}

and

\begin{equation}\label{eq:simp2}
  \sqrt{k^2 - \sin^2\frac{\theta}{2}} = \sqrt{k^2 - k^2\sin^2\phi} = k\cos\phi.
\end{equation}

Then the integral reduces to

$$
  \frac{T}{4}
  =
  \sqrt{\frac{l}{g}}
  \int_0^{\pi/2}
  \frac{d\phi}{\sqrt{1 - k^2\sin^2\phi}}.
$$

Thus the exact period is

\begin{equation}\label{eq:elliptic-integral-form}
  T(\theta_0) = 4\sqrt{\frac{l}{g}}
  K\left(\sin\frac{\theta_0}{2}\right).
\end{equation}

where

$$
  K(k) = \int_0^{\pi/2} \frac{d\phi}{\sqrt{1 - k^2\sin^2\phi}}
$$

is there complete elliptic integral of the first kind.

In the small-angle limit $\theta_0 \to 0$, $k \to 0$ and $K(0) = \pi/2$, so Eq. \eqref{eq:elliptic-integral-form} reduces to

$$
  T_0 \approx 2\pi\sqrt{\frac{l}{g}} = \frac{2\pi}{\omega_n}
$$

which is exactly what we predicted earlier when we invoked the small angle approximation to determine our period.

Furthermore, $K(k)$ can be calculated efficiently using the arithmetic-geometric mean:

$$
  K(k) = \frac{\pi}{2\text{agm}(1, \sqrt{1-k^2})}
$$

where $\text{agm}(x, y)$ is defined as the mutual limit of the arithmetic and geometric mean sequences:

$$
  a_0 = x \qquad g_0 = y
$$

$$
  a_{n+1} = \frac{1}{2}(a_n + g_n) \qquad g_{n+1} = \sqrt{a_n g_n}
$$

These converge very quickly. In Python:

```python
def elliptic_K_agm(k: float, tol: float = 1e-12) -> float:
    if not (0.0 <= k < 1.0):
        raise ValueError("AGM method requires 0 <= k < 1")

    a = 1.0
    b = float(np.sqrt(1.0 - k * k))

    while abs(a - b) > tol:
        a, b = (a + b) / 2.0, float(np.sqrt(a * b))

    return float(np.pi / (2.0 * a))
```

This breaks for $k \ge 1$ (i.e., $\theta_0 \ge \pi$) because the elliptic integral diverges there, which is exactly what we should expect physically.

Now here is the period-vs-amplitude plot comparing the nonlinear period to the small-angle prediction (plus percent error):

[View 'assets/code/physics/pendulums/simple-pendulum/simulate_and_animate.py' on GitHub](https://github.com/carret1268/carret1268.github.io/blob/main/assets/code/physics/pendulums/simple-pendulum/simulate_and_animate.py)

![Period vs. amplitude of simple pendulum]({{ '/assets/media/physics/pendulums/simple-pendulum/figures/period_vs_amplitude.png' | relative_url }})

Increasing the pendulum length uniformly scales the period by $\sqrt{l}$, shifting all curves upward without changing their shape. The divergence between the small-angle approximation and the exact nonlinear period depends only on the initial angle, not on the length.

We see that there is < 5% error at an amplitude of 30 degrees, and just over 10% error at 45 degrees. The error grows faster as the amplitude approaches 180 degrees. The divergence near 180 degrees is a direct manifestation of the separatrix: the system spends an ever-increasing amount of time near the unstable upright equilibrium. In more direct terms, as $\theta_0 \rightarrow \pi$, $T(\theta_0) \rightarrow \infty$.

---

## $\oint$ The Separatrix

The separatrix is the phase-space trajectory corresponding to the critical energy at which the pendulum just reaches the upright position. It divides bounded oscillations from unbounded rotations and marks a qualitative change in the system’s dynamics. Motion along the separatrix takes infinite time, which explains the divergence of the pendulum period as the initial angle approaches 180 degrees. In phase space, the separatrix appears as a boundary curve separating closed and open energy contours.

Recall Eq. \eqref{eq:consv-energy}:

$$
  E = \frac{1}{2}ml^2\omega^2 + mgl(1 - \cos\theta)
$$

The critical energy is the energy needed to reach the top ($\theta=\pi$) with zero speed:

$$
  E_{\text{sep}} = 2mgl
$$

There are 3 cases to consider.

### Separatrix: $E < E_{\text{sep}}$

This regime, which we've already seen plenty of examples for, is characterized by:

- Closed-loop phase diagram
- Oscillatory motion
- Finite period

### Separatrix: $E = E_{\text{sep}}$

We have already seen the system’s energy approach the separatrix energy as the initial angle approaches $\theta_0 \approx 180^\circ$. But what happens if we give the system just enough energy such that when the pendulum reaches $\theta=\pi$ radians, it has $\omega = 0$?

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/separatrix_pend-0-115.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

At this critical energy, the phase diagram undergoes a qualitative change: instead of closed loops, the trajectory becomes a non-closed curve that asymptotically approaches the unstable fixed point $(\theta=\pi,\ \omega=0)$. This trajectory is the separatrix.

We now show explicitly that motion along the separatrix requires infinite time.

Recall Eq. \eqref{eq:dottheta}:

$$  
  \dot{\theta}
  =
  \pm \sqrt{\frac{2g}{\ell}\bigl(\cos\theta - \cos\theta_0\bigr)} .
$$

At the separatrix, the turning point is the upright position,
$\theta_0 = \pi$, for which $\cos\theta_0 = -1$. Equation
\eqref{eq:dottheta} therefore reduces to

$$
  \dot{\theta}
  =
  \pm \sqrt{\frac{2g}{\ell}\bigl(\cos\theta + 1\bigr)} .
$$

Using the trigonometric identity
$$
  \cos\theta + 1 = 2\cos^2\left(\frac{\theta}{2}\right) ,
$$
we obtain

$$
  \dot{\theta}
  =
  \pm 2\sqrt{\frac{g}{\ell}}\cos\frac{\theta}{2} .
$$

To analyze the motion near the top, write
$$
  \theta(t) = \pi - \delta(t) ,
$$
where $\delta(t) \ll 1$ measures the small angular deviation from the
upright position. Differentiating gives
$$
  \dot{\delta} = -\dot{\theta} .
$$

Next, expand the cosine term for small $\delta$:

$$
  \cos\frac{\theta}{2}
  =
  \cos\left(\frac{\pi}{2} - \frac{\delta}{2}\right)
  =
  \sin\frac{\delta}{2}
  \approx \frac{\delta}{2} ,
  \qquad (\delta \ll 1).
$$

Substituting into the expression for $\dot{\theta}$ yields

$$
  \dot{\theta}
  \approx
  \pm \sqrt{\frac{g}{\ell}}\,\delta .
$$

Using $\dot{\delta} = -\dot{\theta}$, we obtain

$$
  \dot{\delta}
  =
  \mp \sqrt{\frac{g}{\ell}}\,\delta .
$$

Since the pendulum approaches the top as $\delta \to 0^+$, we must have
$\dot{\delta} < 0$, which selects the physical branch

$$
  \dot{\delta}
  =
  -\sqrt{\frac{g}{\ell}}\,\delta .
$$

This is a first-order separable differential equation:

$$
  \frac{1}{\delta}\,d\delta
  =
  -\sqrt{\frac{g}{\ell}}\,dt .
$$

Integrating both sides gives

$$
  \ln\vert\delta\vert
  =
  -\sqrt{\frac{g}{\ell}}\,t + C ,
$$

or equivalently,

$$
  \delta(t)
  =
  \delta_0\,e^{-\sqrt{g/\ell}\,t} ,
$$

where $\delta_0 = \delta(t=0)$.

Thus, the angular distance from the upright position decays exponentially but never vanishes at any finite time. In the limit $\delta \to 0$, we require $t \to \infty$. That is, motion along the separatrix takes an infinite amount of time to reach $\theta = \pi$, explaining the divergence of the nonlinear pendulum period.

### Separatrix: $E > E_{\text{sep}}$

For $E > E_{\text{sep}}$, the motion is rotational. Starting from Eq. \eqref{eq:naturalish-energy} and solving for $\omega=\dot{\theta}$ gives

$$
  \omega(\theta) = \pm \sqrt{ \frac{2g}{l} (\varepsilon - 1 + \cos\theta) },
$$

with $\varepsilon = \frac{E}{mgl}$. Since $E > 2mgl$, we have $\varepsilon > 2$. Because $\cos\theta \in [-1,1]$:

$$
  \varepsilon - 1 + \cos\theta \in [\varepsilon - 2, \varepsilon],
$$

which is strictly positive. Therefore:

- $\omega(\theta)$ is defined for all $\theta$
- $\omega$ never hits 0 (no turning points)
- the curve has two smooth branches, + and - (counter-clockwise and clockwise)
- it repeats in $\theta$ with a period of $2\pi$

We can bound $\vert\omega\vert$ in terms of energy:

$$
  \sqrt{ \frac{2g}{l}(\varepsilon - 2) } \leq \vert \omega(\theta) \vert \leq \sqrt{\frac{2g}{l}\varepsilon}
$$

So as $\theta$ increases, the phase curve looks like a smooth, wavy strip oscillating between these bounds, repeating every $2\pi$ radians.

<video
  class="lazy-video"
  controls
  preload="auto"
  playsinline
  muted
  width="85%"
  data-src="/assets/media/physics/pendulums/simple-pendulum/videos/open_pend-0-115.mp4"
>
  Sorry, your browser does not support embedded videos.
</video>

---

## Summary

Here’s the quick recap:

- We derived the exact equation of motion from Newton’s 2nd law:

  $$
    \frac{d^2\theta}{dt^2} = -\omega_n^2\sin\theta
  $$

- Using $\sin\theta\approx\theta$, we got the small-angle harmonic oscillator and the closed-form solution:
  
  $$
    \theta(t) = \theta_0\cos(\omega_n t) + \frac{\omega_0}{\omega_n}\sin(\omega_n t)
  $$

- That immediately gives the classic small-angle period:
  $$
  T_0 = 2\pi\sqrt{\frac{l}{g}}
  $$
- We implemented the small-angle solution and an RK4 solver for the full nonlinear system, then compared them across amplitudes. The big takeaway: as amplitude increases, the exact period grows and the phase diagram deforms away from an ellipse.
- We explained the small-angle ellipse in phase space two ways: (1) parametric sine/cosine geometry and (2) conservation of energy producing quadratic energy contours.
- We derived the exact nonlinear period:
  
  $$
    T(\theta_0) = 4\sqrt{\frac{l}{g}}
    K\left(\sin\frac{\theta_0}{2}\right)
  $$

  and showed how to compute $K(k)$ efficiently via the AGM.
- Finally, we connected the divergence of $T(\theta_0)$ as $\theta_0\to\pi$ to the separatrix: the system approaches the upright equilibrium asymptotically, requiring infinite time in the ideal model.
