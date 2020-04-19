#include "pbd/constraint/DistanceCS.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include <tessellation/Painter.h>
#include <unirender2/RenderState.h>
#include <painting2/RenderSystem.h>

namespace pbd
{
namespace constraint
{

DistanceCS::DistanceCS(double distance, int first, int second, bool st)
    : m_d(distance)
    , m_i1(first)
    , m_i2(second)
    , m_stabile(st)
{
}

DistanceCS::DistanceCS(int first, int second, const std::vector<std::unique_ptr<Particle>>& particles)
    : m_d(0.0)
    , m_i1(first)
    , m_i2(second)
{
    m_d = glm::length(particles[m_i1]->p - particles[m_i2]->p);
}

void DistanceCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    if (p1->imass == 0.f && p2->imass == 0.f) {
        return;
    }

    glm::dvec2 diff = p1->ep - p2->ep;
    double wSum = p1->imass + p2->imass,
            dist = glm::length(diff),
            mag = dist - m_d,
            scale = mag / wSum;

    glm::dvec2 dp = (scale / dist) * diff,
              dp1 = -p1->imass * dp / (double)counts[m_i1],
              dp2 = p2->imass * dp / (double)counts[m_i2];

    p1->ep += dp1;
    p2->ep += dp2;
}

void DistanceCS::Draw(const ur2::Device& dev, ur2::Context& ctx,
                      const std::vector<std::unique_ptr<Particle>>& particles)
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

    ur2::RenderState rs;
    pt2::RenderSystem::DrawPainter(dev, ctx, rs, pt);
}

double DistanceCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    return glm::length(p1->GetP(m_stabile) - p2->GetP(m_stabile)) - m_d;
}

glm::dvec2 DistanceCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    if (!(respect == m_i1 || respect == m_i2)) {
        return glm::dvec2();
    }

    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    glm::dvec2 n = glm::normalize(p1->GetP(m_stabile) - p2->GetP(m_stabile));
    if (respect == m_i1) {
        return n;
    } else {
        return -n;
    }
}

void DistanceCS::UpdateCounts(std::vector<int>& counts)
{
    counts[m_i1]++;
    counts[m_i2]++;
}

}
}