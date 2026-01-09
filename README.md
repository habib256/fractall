# fractall

A portable fractal explorer written in C using SDL.

![License](https://img.shields.io/badge/license-GPL--2.0-blue.svg)

## Features

- **15 fractal types**: Von Koch, Dragon, Mandelbrot, Julia, Newton, Phoenix, Burning Ship, Tricorn, and more
- **7 color palettes**: Normal, Monochrome, Fire, Ocean, Rainbow, Smooth Fire, Smooth Ocean
- **Smooth coloring**: Eliminates visible color bands using continuous iteration interpolation
- **Interactive zoom**: Click to zoom in/out, explore fractal details
- **Status bar**: Displays fractal type, palette, zoom level, coordinates, and render time

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
| **F1-F11** | Select fractal type |
| **C** | Cycle color palette |
| **S** | Screenshot (Screenshot.bmp) |
| **Q** / **ESC** | Quit |
| **Left click** | Zoom in / +1 iteration (vector fractals) |
| **Right click** | Zoom out / -1 iteration (vector fractals) |

## Fractal Types

### Vector Fractals (F1-F2)
- **Von Koch** - Snowflake fractal
- **Dragon** - Dragon curve

### Escape-Time Fractals (F3-F11)
- **Mandelbrot** - The classic z = z² + c
- **Julia** - Julia set with configurable seed
- **Julia Sin** - z = c * sin(z)
- **Newton** - Newton's method fractal
- **Phoenix** - Phoenix fractal
- **Sierpinski** - Sierpinski gasket
- **Burning Ship** - Uses absolute values
- **Tricorn** - Uses complex conjugate
- **Mandelbulb** - 2D slice of 3D fractal (power 8)

## Color Palettes

| Palette | Description |
|---------|-------------|
| Normal | Rainbow gradient (R, G, B cycling) |
| Mono | Grayscale |
| Fire | Black → Red → Yellow → White |
| Ocean | Black → Blue → Cyan → White |
| Rainbow | HSV rainbow with smooth coloring |
| SmoothFire | Fire with smooth gradients |
| SmoothOcean | Ocean with smooth gradients |

The "Smooth" palettes use continuous iteration interpolation based on |z| to eliminate visible color banding.

## Screenshots

Press **S** to save the current view as `Screenshot.bmp`.

## License

GPL-2.0 - See LICENSE file for details.

## Author

Arnaud VERHILLE (2001-2003)
