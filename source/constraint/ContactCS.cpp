#include "pbd/constraint/ContactCS.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include <tessellation/Painter.h>
#include <painting2/RenderSystem.h>

namespace pbd
{
namespace constraint
{

ContactCS::ContactCS(int first, int second, bool st)
    : m_i1(first)
    , m_i2(second)
    , m_stabile(st)
{
}

void ContactCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    if (p1->tmass == 0.f && p2->tmass == 0.f) {
        return;
    }

    glm::dvec2 diff = p1->GetP(m_stabile) - p2->GetP(m_stabile);
    double wSum = p1->tmass + p2->tmass,
            dist = glm::length(diff),
            mag = dist - PARTICLE_DIAM;

    // Previous iterations have moved particles out of collision
    if (mag > 0) {
        return;
    }

    double scale = mag / wSum;
    glm::dvec2 dp = (scale / dist) * diff,
            dp1 = -p1->tmass * dp / (double)counts[m_i1],
              dp2 = p2->tmass * dp / (double)counts[m_i2];

    p1->ep += dp1;
    p2->ep += dp2;

    if (m_stabile) {
        p1->p += dp1;
        p2->p += dp2;
    }
}

void ContactCS::Draw(const std::vector<std::unique_ptr<Particle>>& particles)
{
    auto& p1 = particles[m_i1];
    auto& p2 = particles[m_i2];

    tess::Painter pt;
    pt.SetAntiAliased(false);

    sm::mat4 mt = sm::mat4::Scaled(RENDER_SCALE, RENDER_SCALE, RENDER_SCALE);
    const sm::vec2 pos1 = mt * sm::vec2(static_cast<float>(p1->p.x), static_cast<float>(p1->p.y));
    const sm::vec2 pos2 = mt * sm::vec2(static_cast<float>(p2->p.x), static_cast<float>(p2->p.y));

    pt.AddLine(pos1, pos2, 0xff00ffff);

    const float half_point_size = 1.5f;
    pt.AddRectFilled(pos1, half_point_size, 0xffff00ff);
    pt.AddRectFilled(pos2, half_point_size, 0xffff00ff);

    pt2::RenderSystem::DrawPainter(pt);
}

double ContactCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    double dist = glm::length(p1->GetP(m_stabile) - p2->GetP(m_stabile));
    return dist > PARTICLE_DIAM ? 0 : dist - PARTICLE_DIAM;
}

glm::dvec2 ContactCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    if (!(respect == m_i1 || respect == m_i2)) {
        return glm::dvec2();
    }

    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    glm::dvec2 diff = p1->GetP(m_stabile) - p2->GetP(m_stabile);
    double dist = glm::length(diff);

    if (dist > PARTICLE_DIAM) {
        return glm::dvec2();
    }

    glm::dvec2 n = diff / dist;
    if (respect == m_i1) {
        return n;
    } else {
        return -n;
    }
}

void ContactCS::UpdateCounts(std::vector<int>& counts)
{
    counts[m_i1]++;
    counts[m_i2]++;
}

}
}