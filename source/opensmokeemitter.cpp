#include "pbd/OpenSmokeEmitter.h"
#include "pbd/Particle.h"
#include "pbd/config.h"
#include "pbd/constraint/GasCS.h"

namespace pbd
{

OpenSmokeEmitter::OpenSmokeEmitter(glm::dvec2 posn, double particlesPerSec,
                                   const std::shared_ptr<constraint::GasCS>& gs)
    : m_posn(posn)
    , m_particlesPerSec(particlesPerSec)
    , m_gs(gs)
{
    timer = 0;
}

void OpenSmokeEmitter::Tick(std::vector<std::unique_ptr<Particle>>& estimates, double secs) {
    timer += secs;
    while(timer >= 1./m_particlesPerSec)
    {
        timer -= 1./m_particlesPerSec;
        auto p = std::make_unique<Particle>(m_posn, .1, GAS);
        m_particles.push_back(std::move(p));
        if(m_gs != NULL)
        {
            p = std::make_unique<Particle>(m_posn, 1, GAS);
            m_gs->AddParticle(estimates.size());
            estimates.push_back(std::move(p));
        }
    }
    for(auto& p: m_particles)
    {
        if(p->ph == FLUID || p->ph == GAS)
        {
            p->v = glm::dvec2();
            double sum = 0;
            for(auto& n: estimates)
            {
                glm::dvec2 r = p->p - n->p;
                double p6 = Poly6(glm::dot(r,r));
                p->v += n->v * p6;
                sum += p6;
            }

            if (sum > 0) {
                p->p += p->v * secs / sum;
            }
        }
    }
}

double OpenSmokeEmitter::Poly6(double r2)
{
    if(r2 >= H2) return 0;
    double term2 = (H2 - r2);
    return (315. / (64. * M_PI * H9)) * (term2 * term2 * term2);
}

glm::dvec2 OpenSmokeEmitter::SpikyGrad(const glm::dvec2 &r, double rlen2)
{
    if(rlen2 >= H) return glm::dvec2();
    if(rlen2 == 0) return glm::dvec2();
    return -glm::normalize(r) * (45. / (M_PI * H6)) * (H - rlen2) * (H - rlen2);
//    return -r / (H*H*rlen);
}

}