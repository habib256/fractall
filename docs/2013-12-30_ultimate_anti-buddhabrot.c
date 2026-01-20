// Ultimate Anti-Buddhabrot (c) 2013 Claude Heiland-Allen
// http://mathr.co.uk/blog/2013-12-30_ultimate_anti-buddhabrot.html
// gcc -std=c99 -Wall -pedantic -Wextra -O3 -fopenmp -lm -o anti anti.c
// ./anti > anti.ppm

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


typedef unsigned int N;
typedef unsigned long int NN;
typedef int Z;
typedef long int ZZ;
typedef double R;
typedef double complex C;
typedef float F;


static const R pi = 3.141592653589793;


static const Z zwidth = 4096;
static const Z zheight = 4096;
static const C zpixel_size = 5.0 / 4096;
static const C z0 = 0;


static F count[4096][4096][4];
static ZZ total = 0;


R cabs2(C z) { return creal(z) * creal(z) + cimag(z) * cimag(z); }


static void hsv2rgb(F h, F s, F v, F *rp, F *gp, F *bp) {
  F i, f, p, q, t, r, g, b;
  Z ii;
  if (s == 0.0) {
    r = v;
    g = v;
    b = v;
  } else {
    h = h - floor(h);
    h = h * 6.0;
    i = floor(h);
    ii = i;
    f = h - i;
    p = v * (1 - s);
    q = v * (1 - (s * f));
    t = v * (1 - (s * (1 - f)));
    switch(ii) {
      case 0:  r = v; g = t; b = p; break;
      case 1:  r = q; g = v; b = p; break;
      case 2:  r = p; g = v; b = t; break;
      case 3:  r = p; g = q; b = v; break;
      case 4:  r = t; g = p; b = v; break;
      default: r = v; g = p; b = q; break;
    }
  }
  *rp = r;
  *gp = g;
  *bp = b;
}


static void plot_iterates(C z1, C c, N period, R r, R g, R b) {
  for (N j = 0; j < 2; ++j) {
    C z = z1;
    for (N i = 0; i < period; ++i) {
      C pz0 = (zwidth / 2.0 - 0.5) + I * (zheight / 2.0 - 0.5);
      C pz = (z - z0) / zpixel_size + pz0;
      Z x = creal(pz);
      Z y = cimag(pz);
      if (0 <= x && x < zwidth && 0 <= y && y < zheight) {
        #pragma omp atomic update
        count[y][x][0] += r;
        #pragma omp atomic update
        count[y][x][1] += g;
        #pragma omp atomic update
        count[y][x][2] += b;
        #pragma omp atomic update
        count[y][x][3] += 1;
        #pragma omp atomic update
        total++;
      }
      z = z * z + c;
    }
    z1 = conj(z1);
    c = conj(c);
  }
}


static void clear_image() {
  memset(&count[0][0][0], 0, zwidth * zheight * 4 * sizeof(F));
  total = 0;
}


static void output_image() {
  R s = zwidth * zheight / (R) total;
  printf("P6\n%u %u\n255\n", zwidth, zheight);
  for (Z y = 0; y < zheight; ++y) {
    for (Z x = 0; x < zwidth; ++x) {
      if (count[y][x][3]) {
        R v = log2(1 + count[y][x][3] * s);
        for (Z c = 0; c < 3; ++c) {
          R u = count[y][x][c] / count[y][x][3];
          Z o = fminf(fmaxf(54 * v * u, 0), 255);
          putchar(o);
        }
      } else {
        putchar(0);
        putchar(0);
        putchar(0);
      }
    }
  }
}



int wucleus(C *z0, C c, N period) {
  R eps = 1e-12;
  R er2 = 16;
  C z = *z0;
  for (N j = 0; j < 256; ++j) {
    C dz = 1.0;
    for (N k = 0; k < period; ++k) {
      dz = 2.0 * dz * z;
      z = z * z + c;
    }
    R z2 = cabs(z);
    if (! (z2 < er2)) {
      break;
    }
    z = *z0 - (z - *z0) / (dz - 1.0);
    R e = cabs(z - *z0);
    *z0 = z;
    if (e < eps) {
      return 1;
    }
  }
  return 0;
}

R interior_distance(C *w, C c, N period, R pixel_size) {
  if (wucleus(w, c, period)) {
    C z = *w;
    C dz = 1.0;
    C dzdz = 0.0;
    C dc = 0.0;
    C dcdz = 0.0;
    for (N j = 0; j < period; ++j) {
      dcdz = 2.0 * (z * dcdz + dz * dc);
      dc = 2.0 * z * dc + 1.0;
      dzdz = 2.0 * (dz * dz + z * dzdz);
      dz = 2.0 * z * dz;
      z = z * z + c;
    }
    return (1.0 - cabs2(dz)) / (cabs(dcdz + dzdz * dc / (1.0 - dz)) * pixel_size);
  }
  return -1.0;
}

static void render_recursive_interior(F red, F grn, F blu, N period, C z, C c, R grid_spacing, N depth) {
  if (depth == 0) { return; }
  C z1 = z;
  R de = interior_distance(&z1, c, period, grid_spacing);
  if (de > 0) {
    plot_iterates(z1, c, period, red, grn, blu);
    render_recursive_interior(red, grn, blu, period, z1, c + 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
    render_recursive_interior(red, grn, blu, period, z1, c - 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
    render_recursive_interior(red, grn, blu, period, z1, c + I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
    render_recursive_interior(red, grn, blu, period, z1, c - I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
  }
}

static void render_recursive_unknown(C c, R grid_spacing, N depth) {
  if (depth == 0) { return; }
  C z = 0;
  C dz = 0;
  R mz2 = 1.0/0.0;
  N maxiters = 1 << 12;
  for (N i = 1; i < maxiters; ++i) {
    dz = 2 * z * dz + 1;
    z = z * z + c;
    R z2 = cabs2(z);
    if (! (z2 < 65536)) {
      R de = sqrt(z2) * log(z2) / (cabs(dz) * grid_spacing);
      if (de < 4) {
        render_recursive_unknown(c + 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
        render_recursive_unknown(c - 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
        render_recursive_unknown(c + I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
        render_recursive_unknown(c - I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
      }
      break;
    }
    if (z2 < mz2) {
      mz2 = z2;
      C z1 = z;
      R de = interior_distance(&z1, c, i, grid_spacing);
      if (de > 0) {
        R phi5 = pow((sqrt(5) + 1) / 2, 5);
        R hue = i / phi5 ;
        F red, grn, blu;
        hsv2rgb(hue, 1, 1, &red, &grn, &blu);
        plot_iterates(z1, c, i, red, grn, blu);
        if (de > 0.25) {
          render_recursive_interior(red, grn, blu, i, z1, c + 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_interior(red, grn, blu, i, z1, c - 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_interior(red, grn, blu, i, z1, c + I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_interior(red, grn, blu, i, z1, c - I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
        } else {
          render_recursive_unknown(c + 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_unknown(c - 1 * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_unknown(c + I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
          render_recursive_unknown(c - I * 0.5 * grid_spacing, 0.5 * grid_spacing, depth - 1);
        }
        break;
      }
    }
  }
}

static void render() {
  clear_image();
  N grid = 1024;
  R grid_spacing = 5.0 / grid;
  #pragma omp parallel for schedule(dynamic, 1)
  for (N y = 0; y < grid; ++y) {
    for (N x = 0; x < grid; ++x) {
      C c = grid_spacing * ((x + 0.5 - grid/2.0) + I * (y + 0.5 - grid/2.0));
      if (cimag(c) >= 0) {
        render_recursive_unknown(c, grid_spacing, 9);
      }
    }
    #pragma omp critical
    fprintf(stderr, "%8d\r", y);
  }
  output_image();
}


extern int main(int argc, char **argv) {
  (void) argc;
  (void) argv;
  render();
  return 0;
}
