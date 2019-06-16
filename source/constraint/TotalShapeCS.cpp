#include "pbd/constraint/TotalShapeCS.h"
#include "pbd/Particle.h"
#include "pbd/Body.h"
#include "pbd/config.h"

#include <tessellation/Painter.h>
#include <painting2/RenderSystem.h>

namespace pbd
{
namespace constraint
{

TotalShapeCS::TotalShapeCS(const std::shared_ptr<Body>& bod, double stiff)
    : m_body(bod)
{
    m_stiffness = stiff;
}

void TotalShapeCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    m_body->UpdateCOM(estimates);

    // implemented using http://labs.byhook.com/2010/06/29/particle-based-rigid-bodies-using-shape-matching/
    for (size_t i = 0; i < m_body->particles.size(); i++) {
        int idx = m_body->particles[i];
        auto& p = estimates[idx];
        p->ep += (Guess(idx) - p->ep) * m_stiffness;
    }
}

void TotalShapeCS::Draw(const std::vector<std::unique_ptr<Particle>>& particles)
{
    tess::Painter pt;
    pt.SetAntiAliased(false);

    sm::mat4 mt = sm::mat4::Scaled(RENDER_SCALE, RENDER_SCALE, RENDER_SCALE);
    const sm::vec2 center = mt * sm::vec2(static_cast<float>(m_body->center.x), static_cast<float>(m_body->center.y));
    for (size_t i = 0; i < m_body->particles.size(); i++)
    {
        int idx = m_body->particles[i];
        auto& p = particles[idx];
        auto pos = mt * sm::vec2(static_cast<float>(p->p.x), static_cast<float>(p->p.y));
        pt.AddLine(pos, center, 0xff00ff00);
    }

    const float half_point_size = 1.5f;
    pt.AddRectFilled(center, half_point_size, 0xffff00ff);
    for (size_t i = 0; i < m_body->particles.size(); i++)
    {
        int idx = m_body->particles[i];
        auto& p = particles[idx];
        auto pos = mt * sm::vec2(static_cast<float>(p->p.x), static_cast<float>(p->p.y));
        pt.AddRectFilled(pos,half_point_size, 0xffff00ff);
    }
}

double TotalShapeCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    (void) estimates;
    return 0;
}

glm::dvec2 TotalShapeCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
//    if (body->rs.contains(respect)) {
//        Particle *p = estimates->at(respect);
//        glm::dvec2 out = Guess(respect) - p->ep;
//        if (out == glm::dvec2()) {
//            return glm::dvec2(0,0);
//        }
//        return -glm::normalize(out);
//    }
    (void) estimates;
    (void) respect;
    return glm::dvec2();
}

void TotalShapeCS::UpdateCounts(std::vector<int>& counts)
{
    for (size_t i = 0; i < m_body->particles.size(); i++) {
        counts[m_body->particles[i]]++;
    }
}

glm::dvec2 TotalShapeCS::Guess(int idx)
{
    double c = cos(m_body->angle), s = sin(m_body->angle);

    glm::dvec2 q = m_body->rs[idx],
               d = glm::dvec2(c * q.x - s * q.y, s * q.x + c * q.y);
    return d + m_body->center;
}

}
}