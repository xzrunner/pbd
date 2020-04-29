#pragma once

#include "pbd/Constraint.h"

namespace pbd
{

struct Body;

namespace constraint
{

class TotalShapeCS : public Constraint
{
public:
    TotalShapeCS(const std::shared_ptr<Body>& bod, double stiff = 1.0);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur::Device& dev, ur::Context& ctx,
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

    glm::dvec2 Guess(int idx);

private:
    std::shared_ptr<Body> m_body = nullptr;

}; // TotalShapeCS

}
}