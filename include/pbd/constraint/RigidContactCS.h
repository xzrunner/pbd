#pragma once

#include "pbd/Constraint.h"

namespace pbd
{

struct Body;

namespace constraint
{

class RigidContactCS : public Constraint
{
public:
    RigidContactCS(int first, int second, const std::vector<std::shared_ptr<Body>>& bodies, bool st = false);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur::Device& dev, ur::Context& ctx, 
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

    bool InitBoundary(const std::unique_ptr<Particle>& p1, const std::unique_ptr<Particle>& p2);

private:
    const std::vector<std::shared_ptr<Body>>& m_bods;
    glm::dvec2 m_n;
    double m_d = 0;
    int m_i1 = 0, m_i2 = 0;
    bool m_stabile = false;

}; // RigidContactCS

}
}