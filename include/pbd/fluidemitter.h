#pragma once

#include <glm/vec2.hpp>

#include <memory>
#include <vector>

namespace pbd
{

struct Particle;
namespace constraint { class TotalFluidCS; }

class FluidEmitter
{
public:
    FluidEmitter(glm::dvec2 posn, double particlesPerSec,
        const std::shared_ptr<constraint::TotalFluidCS>& fs);

    void Tick(std::vector<std::unique_ptr<Particle>>& estimates, double secs);
    auto& GetParticles() { return m_grains; }
    auto& GetPosn() const { return m_posn; }

private:
    glm::dvec2 m_posn;
    double m_particlesPerSec = 0;
    double m_timer = 0;
    double m_tot_timer = 0;
    std::shared_ptr<constraint::TotalFluidCS> m_fs = nullptr;
    std::vector<std::unique_ptr<Particle>> m_grains;

}; // FluidEmitter

}
