#pragma once

#include <glm/vec2.hpp>
#include <glm/gtx/rotate_vector.hpp>

namespace pbd
{

// Signed distance field data for rigid-body collisions
struct SDFData
{
    SDFData()
        : gradient(glm::dvec2()), distance(-1.0) {}

    SDFData(const glm::dvec2 grad, double dist)
        : gradient(grad), distance(dist) {}

    inline void rotate(double angle) { gradient = glm::rotate(gradient, angle); }

    glm::dvec2 gradient;
    double distance;
};

}