#pragma once

#include "pbd/Constraint.h"

namespace pbd
{
namespace constraint
{

// Contact between two particles where AT LEAST ONE is not a solid
class ContactCS : public Constraint
{
public:
    ContactCS(int first, int second, bool st = false);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur2::Device& dev, ur2::Context& ctx,
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

private:
    int m_i1 = 0, m_i2 = 0;
    bool m_stabile = false;

}; // ContactCS

}
}