#pragma once

#include "pbd/Constraint.h"

#include <unordered_map>

namespace pbd
{
namespace constraint
{

class GasCS : public Constraint
{
public:
    GasCS(double density, const std::vector<int>& particles, bool open);

    virtual void Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts) override;
    virtual void Draw(const ur2::Device& dev, ur2::Context& ctx, 
        const std::vector<std::unique_ptr<Particle>>& particles) override;

    virtual double Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates) override;
    virtual glm::dvec2 Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect) override;
    virtual void UpdateCounts(std::vector<int>& counts) override;

    double Poly6(double rlen);
    glm::dvec2 SpikyGrad(const glm::dvec2 &r, double rlen);
    glm::dvec2 Grad(const std::vector<std::unique_ptr<Particle>>& estimates, int k, int j);

    void AddParticle(int index);

private:
    double m_p0;
    std::vector<int> m_ps;
    std::vector<std::unique_ptr<std::vector<int>>> m_neighbors;
    int m_num_particles;
    std::vector<glm::dvec2> m_deltas;
    std::unordered_map<int, double> m_lambdas;
    bool m_open;

}; // GasCS

}
}