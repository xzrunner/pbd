#pragma once

#include <glm/vec2.hpp>

#include <vector>
#include <memory>

namespace ur { class Device; class Context; }

namespace pbd
{

struct Particle;

// Groups of constraints to be solved together
enum ConstraintGroup
{
    STABILIZATION,
    CONTACT,
    STANDARD,
    SHAPE,
    NUM_CONSTRAINT_GROUPS
};

// Abstract superclass of all constraint types
class Constraint
{
public:
    Constraint() {}
    virtual ~Constraint() {}

    virtual void Draw(const ur::Device& dev, ur::Context& ctx,
        const std::vector<std::unique_ptr<Particle>>& particles) = 0;

    // For iterative solving of constraints
    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) = 0;

    // For matrix-oriented solving of constraints
    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) = 0;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) = 0;
    virtual void UpdateCounts(std::vector<int>& counts) = 0;

protected:
    double m_stiffness = 1;

}; // Constraint

}