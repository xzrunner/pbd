#pragma once

#include "pbd/SDFData.h"

#include <glm/vec2.hpp>

#include <vector>
#include <unordered_map>
#include <memory>

namespace pbd
{

struct Particle;
class Constraint;

// A single rigid body
struct Body
{
    std::vector<int> particles; // index into global particles list
    std::unordered_map<int, glm::dvec2> rs; // map from global particles index to r vector
    std::unordered_map<int, SDFData> sdf; // map from global particles index to SDF data
    std::shared_ptr<Constraint> shape = nullptr;
    glm::dvec2 center; // center of mass
    double imass = 0, angle = 0; // total inverse mass

    void UpdateCOM(const std::vector<std::unique_ptr<Particle>>& estimates, bool useEstimates = true);
    void ComputeRs(const std::vector<std::unique_ptr<Particle>>& estimates);

}; // Body

}