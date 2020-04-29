#pragma once

#include "pbd/Constraint.h"

namespace pbd
{
namespace constraint
{

// Collision with the boundaries of the world
class BoundaryCS : public Constraint
{
public:
    BoundaryCS(int index, double val, bool x_boundary, bool greater, bool st = false);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur::Device& dev, ur::Context& ctx, 
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

private:
    int    m_idx = 0;
    double m_value = 0;

    bool m_is_x = false, m_is_greater_than = false, m_stabile = false;

}; // BoundaryCS

}
}
