#pragma once

#include <glm/vec2.hpp>

#include <memory>
#include <vector>

namespace pbd
{

struct Particle;
namespace constraint { class GasCS; }

class OpenSmokeEmitter
{
public:
    OpenSmokeEmitter(glm::dvec2 posn, double particlesPerSec,
        const std::shared_ptr<constraint::GasCS>& gs);

    void Tick(std::vector<std::unique_ptr<Particle>>& estimates, double secs);
    auto& GetParticles() { return m_particles; }
    auto& GetPosn() const { return m_posn; }

private:
    double Poly6(double r2);
    glm::dvec2 SpikyGrad(const glm::dvec2 &r, double rlen2);

private:
    glm::dvec2 m_posn;
    double m_particlesPerSec = 0;
    std::vector<std::unique_ptr<Particle>> m_particles;
    double timer = 0;
    std::shared_ptr<constraint::GasCS> m_gs = nullptr;

}; // OpenSmokeEmitter

}