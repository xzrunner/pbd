#pragma once

#include "pbd/Constraint.h"

namespace pbd
{
namespace constraint
{

// Two particles must be exactly a certain distance away
class DistanceCS : public Constraint
{
public:
    DistanceCS(double distance, int first, int second, bool st = false);
    DistanceCS(int first, int second, const std::vector<std::unique_ptr<Particle>>& particles);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur::Device& dev, ur::Context& ctx,
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

private:
    double m_d = 0;
    int m_i1 = 0, m_i2 = 0;
    bool m_stabile = false;

}; // DistanceCS

}
}