# Simple Pendulum

This directory contains a small, self-contained simple pendulum simulator plus Matplotlib-based animation utilities used to generate videos and figures for the Jekyll site.

The implementation supports:

- An analytical small-angle approximation (SAA) solution.
- A nonlinear numerical solution using a fourth-order Runge-Kutta (RK4) integrator.
- Animation workflows for single runs and side-by-side comparisons.
- A period-vs-amplitude figure using the complete elliptic integral computed via the AGM method.

## Files

### `simple_pendulum.py`

Core simulation code.

- `SimplePendulum`
  - Stores physical parameters (`l`, `g`) and derived natural frequency `omega_n = sqrt(g/l)`.
  - Provides `simulate(...)` to generate a time series of angle and angular velocity.
  - Supports `algorithm="saa"` (small-angle analytical solution) and `algorithm="rk4"` (numerical nonlinear integration).
  - Estimates the period (if possible) from zero-crossings of `omega(t)`.

- `SimplePendulumResult`
  - Immutable container for simulation outputs.
  - Stores internal state in radians and seconds; convenience properties provide degree outputs.

### `simple_pendulum_animator.py`

Matplotlib animation helpers.

- `PendulumAnimator`
  - Centralizes theme settings (dark styling) and shared plotting helpers.

- `SimplePendulumAnimator`
  - `animate_single_simulation(...)` renders:
    - Pendulum geometry
    - Phase portrait (omega vs theta)
    - Angle vs time
  - `animate_comparison_simulations(...)` renders:
    - Phase portrait for simulation 1 with an inset pendulum overlay
    - Phase portrait for simulation 2 with an inset pendulum overlay
    - Shared angle vs time panel

Output is saved as an `.mp4` via `FFMpegWriter`.

### `simulate_and_animate.py`

Convenience script to generate site assets (videos and figures).

Provides routines to:

- Generate single SAA pendulum animations (`make_saa_pendulums`)
- Generate SAA vs RK4 comparison animations across amplitudes (`make_comparisons`)
- Generate a separatrix example (`make_separatrix_pendulum`)
- Generate an open-trajectory example (`make_open_pendulum`)
- Generate a period-vs-amplitude and percent-error figure (`plot_period_and_percent_error_vs_angle`)

## Output locations

This script writes videos into:

- `assets/media/physics/pendulums/simple-pendulum/videos/`

And (by default) writes the period figure to:

- `assets/media/physics/pendulums/simple-pendulum/figures/period_vs_amplitude.png`

If the `figures/` directory does not already exist, create it.

Recommended layout:

- `assets/`
  - `media/`
    - `physics/`
      - `pendulums/`
        - `simple-pendulum/`
          - `videos/`
          - `figures/`

## Dependencies

Runtime requirements for this directory:

- Python 3.11+
- NumPy
- Matplotlib

Additionally required for saving `.mp4` videos:

- FFmpeg installed and available on PATH

Quick check:

```bash
ffmpeg -version
```
