#include "pbd/Body.h"
#include "pbd/Particle.h"

namespace pbd
{

void Body::UpdateCOM(const std::vector<std::unique_ptr<Particle>>& estimates, bool useEstimates)
{
    // Recompute center of mass
    glm::dvec2 total(0, 0);
    for (size_t i = 0; i < particles.size(); i++) {
        auto& p = estimates[particles[i]];
        total += (useEstimates ? p->ep : p->p) / p->imass;
    }
    center = total * imass;

    // Recompute angle delta guess
    angle = 0.0;
    double prev = 0.0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        int index = particles[i];
        glm::dvec2 q = rs[index];
        if (glm::dot(q,q) == 0) {
            continue;
        }
        auto& p = estimates[index];
        glm::dvec2 r = p->ep - center;

        double cos = r.x * q.x + r.y * q.y,
               sin = r.y * q.x - r.x * q.y,
               next = atan2(sin, cos);

        // Ensure all guesses are close to each other
        if (i > 0) {
            if (prev - next >= M_PI) {
                next += 2 * M_PI;
            }
        } else {
            if (next < 0) {
                next += 2 * M_PI;
            }
        }

        prev = next;
        next /= p->imass;
        angle += next;
    }
    angle *= imass;
}

void Body::ComputeRs(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    imass = 0.0;
    for (size_t i = 0; i < particles.size(); i++)
    {
        int idx = particles[i];
        auto& p = estimates[idx];
        glm::dvec2 r = p->p - center;
        rs[idx] = r;

        if (glm::dot(r,r) != 0) {
            imass += (1.0 / p->imass);
        }
    }
    imass = 1.0 / imass;
}

}