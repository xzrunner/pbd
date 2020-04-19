#include "pbd/constraint/BoundaryCS.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include <algorithm>

namespace pbd
{
namespace constraint
{

BoundaryCS::BoundaryCS(int index, double val, bool x_boundary, bool greater, bool st)
    : m_idx(index)
    , m_value(val)
    , m_is_x(x_boundary)
    , m_is_greater_than(greater)
    , m_stabile(st)
{
}

void BoundaryCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    auto& p = estimates[m_idx];

    // Add a little random jitter for fluids and gases so particles do not become trapped on boundaries
    double extra = p->ph == FLUID || p->ph == GAS ? frand() * .003 : 0;
    double d = (PARTICLE_RAD + extra);
    glm::dvec2 n = glm::dvec2();

    // Move the particle back into a valid spot (if necessary)
    if (m_is_greater_than)
    {
        if (m_is_x)
        {
            // Quit if no longer valid
            if (p->ep.x >= m_value + PARTICLE_RAD) {
                return;
            }
            p->ep.x = m_value + d;
            if (m_stabile) {
                p->p.x = m_value + d;
            }
            n = glm::dvec2(1,0);
        }
        else
        {
            // Quit if no longer valid
            if (p->ep.y >= m_value + PARTICLE_RAD) {
                return;
            }
            p->ep.y = m_value + d;
            if (m_stabile) {
                p->p.y = m_value + d;
            }
            n = glm::dvec2(0,1);
        }
    }
    else
    {
        if (m_is_x)
        {
            // Quit if no longer valid
            if (p->ep.x <= m_value - PARTICLE_RAD) {
                return;
            }
            p->ep.x = m_value - d;
            if (m_stabile) {
                p->p.x = m_value - d;
            }
            n = glm::dvec2(-1,0);
        }
        else
        {
            // Quit if no longer valid
            if (p->ep.y <= m_value - PARTICLE_RAD) {
                return;
            }
            p->ep.y = m_value - d;
            if (m_stabile) {
                p->p.y = m_value - d;
            }
            n = glm::dvec2(0,-1);
        }
    }

    if (m_stabile) {
        return;
    }

    // Apply friction - boundaries have a coefficient of friction of 1
    glm::dvec2 dp = (p->ep - p->p) / (double)counts[m_idx],
               dpt = dp - glm::dot(dp, n) * n;
    double ldpt = glm::length(dpt);

    if (ldpt < EPSILON) {
        return;
    }

    // Choose between static and kinetic friction
    if (ldpt < sqrt(p->sFriction) * d) {
        p->ep -= dpt;
    } else {
        p->ep -= dpt * std::min(sqrt(p->kFriction )* d / ldpt, 1.);
    }
}

void BoundaryCS::Draw(const ur2::Device& dev, ur2::Context& ctx,
                      const std::vector<std::unique_ptr<Particle>>& particles)
{
}

double BoundaryCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    auto& p = estimates[m_idx];
    if (m_is_greater_than)
    {
        if (m_is_x) {
            return (m_value + PARTICLE_RAD) - p->GetP(m_stabile).x;
        } else {
            return (m_value + PARTICLE_RAD) - p->GetP(m_stabile).y;
        }
    }
    else
    {
        if (m_is_x) {
            return p->GetP(m_stabile).x - (m_value - PARTICLE_RAD);
        } else {
            return p->GetP(m_stabile).y - (m_value - PARTICLE_RAD);
        }
    }
}

glm::dvec2 BoundaryCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    if (respect != m_idx) {
        return glm::dvec2();
    }

    if (m_is_greater_than)
    {
        if (m_is_x) {
            return glm::dvec2(-1,0);
        } else {
            return glm::dvec2(0,-1);
        }
    }
    else
    {
        if (m_is_x) {
            return glm::dvec2(1,0);
        } else {
            return glm::dvec2(0,1);
        }
    }
}

void BoundaryCS::UpdateCounts(std::vector<int>& counts)
{
    counts[m_idx]++;
}

}
}