# fractall

A portable fractal explorer written in C using SDL.

![License](https://img.shields.io/badge/license-GPL--2.0-blue.svg)
![Version](https://img.shields.io/badge/version-1.0-green.svg)

## Features

- **17 fractal types**: Von Koch, Dragon, Mandelbrot, Julia, Newton, Phoenix, Burning Ship, Tricorn, Mandelbulb, Buddhabrot, Lyapunov, and more
- **5 smooth color palettes**: SmoothFire, SmoothOcean, SmoothForest, SmoothViolet, SmoothRainbow
- **Smooth coloring**: Eliminates visible color bands using continuous iteration interpolation
- **Interactive zoom**: Click to zoom in/out, explore fractal details
- **Status bar**: Displays fractal type, palette, zoom level, coordinates, and render time
- **Optional GMP support**: High-precision arithmetic for deep zooms

## Building

```bash
./autogen.sh      # If configure doesn't exist
./configure
make
./src/fractall
```

### Dependencies

- SDL 1.2.0+
- Standard math library (`-lm`)
- GMP (optional): For high-precision arithmetic (`./configure --with-gmp`)

## Usage

```bash
fractall [OPTIONS]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-x<n>` | Window width | 800 |
| `-y<n>` | Window height | 600 |
| `-f` | Fullscreen mode | - |
| `-g<n>` | GUI height | 51 |
| `-nogui` | Disable GUI | - |
| `-h` | Show help | - |

## Controls

| Key/Action | Effect |
|------------|--------|
| **F1-F8** | Von Koch, Dragon, Mandelbrot, Julia, Julia Sin, Newton, Phoenix, Sierpinski |
| **F9-F12** | Burning Ship, Tricorn, Mandelbulb, Buddhabrot |
| **GUI buttons** | Select any of the 17 fractal types |
| **C** | Cycle color palette (5 palettes) |
| **S** | Screenshot (Screenshot.bmp) |
| **Q** / **ESC** | Quit |
| **Left click** | Zoom in / +1 iteration (vector fractals) |
| **Right click** | Zoom out / -1 iteration (vector fractals) |

## Fractal Types

### Vector Fractals (F1-F2)
- **Von Koch** - Snowflake fractal
- **Dragon** - Dragon curve

### Escape-Time Fractals (F3-F12 + GUI)
- **Mandelbrot** - The classic z = z² + c
- **Julia** - Julia set with configurable seed
- **Julia Sin** - z = c * sin(z)
- **Newton** - Newton's method fractal
- **Phoenix** - Phoenix fractal
- **Sierpinski** - Sierpinski gasket
- **Barnsley J/M** - Barnsley fractal variants
- **Magnet J/M** - Magnetic fractal variants
- **Burning Ship** - Uses absolute values
- **Tricorn** - Uses complex conjugate
- **Mandelbulb** - 2D slice of 3D fractal (power 8)
- **Buddhabrot** - Density-based rendering of escape trajectories
- **Lyapunov** - Lyapunov exponent of logistic map

## Color Palettes

| Palette | Description |
|---------|-------------|
| SmoothFire | Black → Red → Yellow → White (smooth gradients, 4× repetition) |
| SmoothOcean | Black → Blue → Cyan → White (smooth gradients, 4× repetition) |
| SmoothForest | Black → Green → Yellow → White (smooth gradients, 4× repetition) |
| SmoothViolet | Black → Violet → Pink → White (smooth gradients, 4× repetition) |
| SmoothRainbow | Full rainbow spectrum (Red → Orange → Yellow → Green → Cyan → Blue → Violet) |

All palettes use continuous iteration interpolation based on |z| to eliminate visible color banding. The palettes alternate forward/reverse to avoid abrupt transitions.

**Note**: Buddhabrot and Lyapunov fractals have their own specialized coloring algorithms.

## Screenshots

Press **S** to save the current view as `Screenshot.bmp`.

## License

GPL-2.0 - See LICENSE file for details.

## Version History

- **1.0** (2025) - Added Lyapunov fractal, expanded color palettes, improved rendering
- **0.5.1** (2025) - Added smooth coloring, status bar, bug fixes
- **0.5** (2003) - Initial release

## Author

Arnaud VERHILLE (2001-2026)
