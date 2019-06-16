#pragma once

#include "pbd/config.h"

#include <glm/vec2.hpp>
#include <glm/geometric.hpp>

#include <vector>
#include <memory>

namespace pbd
{

// Phase of mass for particles
enum Phase
{
    SOLID,
    FLUID,
    GAS,
    NUM_PHASES
};

struct Body;
struct SDFData;

// Individual particle representation
struct Particle
{
    glm::dvec2 p, ep, v, f; // position, guess position, and velocity
    double imass, tmass, sFriction, kFriction, t; // inverse mass, temporary height-scaled mass, coeffs of friction
    int bod; // body (if any) this particle belongs to, for disabling collisions
    Phase ph; // phase of this particle

    Particle()
        : p(glm::dvec2()), v(glm::dvec2()), ph(NUM_PHASES) { Init(0); }

    Particle(glm::dvec2 pos, double mass, Phase phase = SOLID)
        : p(pos), v(glm::dvec2()), ph(phase) { Init(mass); }

    Particle(glm::dvec2 pos, glm::dvec2 vel, double mass, Phase phase)
        : p(pos), v(vel), ph(phase) { Init(mass); }

    void Init(double mass) {
        t = 4.;
        ep = glm::dvec2();
        bod = -1;

        if (mass <= 0) {
            imass = -mass;
        } else {
            imass = 1. / mass;
        }
        tmass = imass;

        f = glm::dvec2();
        sFriction = 0;
        kFriction = 0; // usually smaller the coefficient of static friction
    }

    inline void SetStatic() { imass = 0.; }

    inline glm::dvec2 Guess(double seconds) {
        return imass == 0. ? p : p + seconds * v;
    }

    inline void ConfirmGuess() {
        if (glm::length(ep - p) < EPSILON) {
            v = glm::dvec2(0,0); return;
        }
        p = ep;
    }

    void ScaleMass() {
        if (imass != 0.0) {
            tmass = 1. / ((1. / imass) * exp(-p.y));
        } else {
            tmass = 0.0;
        }
    }

    // Used for stabilization-related constraints
    inline glm::dvec2 GetP(bool stabile) { return stabile ? p : ep; }

    SDFData getSDFData(const std::vector<std::shared_ptr<Body>>& bodies, int idx);
};

}