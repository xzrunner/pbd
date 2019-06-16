#pragma once

#include <stdlib.h>

namespace pbd
{

static const float RENDER_SCALE = 50.0f;

static const double H  = 2;
static const double H2 = 4;
static const double H6 = 64;
static const double H9 = 512;

static const double PARTICLE_RAD = .25;
static const double PARTICLE_DIAM = .5;

static const double EPSILON = .0001;
static const double M_PI    = 3.14159265358979323846;   // pi

#ifndef D2R
#define D2R(d) (d * M_PI / 180)
#endif // D2R

#ifndef R2D
#define R2D(r) (r * 180 / M_PI)
#endif // R2D

// Generally helpful functions
inline double frand() { return (double)rand() / (double)RAND_MAX; }
inline double urand(double a, double b) { return a + (b - a) * frand(); }

}