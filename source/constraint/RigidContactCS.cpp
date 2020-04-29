#include "pbd/constraint/RigidContactCS.h"
#include "pbd/Particle.h"
#include "pbd/SDFData.h"
#include "pbd/config.h"

#include <algorithm>

namespace pbd
{
namespace constraint
{

RigidContactCS::RigidContactCS(int first, int second, const std::vector<std::shared_ptr<Body>>& bodies, bool st)
    : m_d(0.0)
    , m_i1(first)
    , m_i2(second)
    , m_stabile(st)
    , m_bods(bodies)
{

}

void RigidContactCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    SDFData dat1 = p1->getSDFData(m_bods, m_i1), dat2 = p2->getSDFData(m_bods, m_i2);

    if (dat1.distance < 0 || dat2.distance < 0) {
        glm::dvec2 x12 = p2->GetP(m_stabile) - p1->GetP(m_stabile);
        double len = glm::length(x12);
        m_d = PARTICLE_DIAM - len;
        if (m_d < EPSILON) return;
        m_n = x12 / len;
    } else {
        if (dat1.distance < dat2.distance) {
            m_d = dat1.distance;
            m_n = dat1.gradient;
        } else {
            m_d = dat2.distance;
            m_n = -dat2.gradient;
        }

        if (m_d < PARTICLE_DIAM + EPSILON) {
            if (InitBoundary(p1, p2)) {
                return;
            }
        }
    }

    double wSum = p1->tmass + p2->tmass;
    glm::dvec2 dp = (1.0 / wSum) * m_d * m_n,
              dp1 = -p1->tmass * dp  / (double)counts[m_i1],
              dp2 = p2->tmass * dp / (double)counts[m_i2];

    if (!m_stabile) {
        p1->ep += dp1;
        p2->ep += dp2;
    } else {
        p1->p += dp1;
        p2->p += dp2;
    }

    // Apply friction
    glm::dvec2 nf = glm::normalize(m_n);
    glm::dvec2 dpf = (p1->ep - p1->p) - (p2->ep - p2->p),
               dpt = dpf - glm::dot(dpf, nf) * nf;
    double ldpt = glm::length(dpt);
    if (ldpt < EPSILON) {
        return;
    }
    double sFric = sqrt(p1->sFriction * p2->sFriction),
            kFric = sqrt(p1->kFriction * p2->kFriction);

    if (ldpt < sFric * m_d) {
        if (m_stabile) {
            p1->p -= dpt * p1->tmass / wSum;
            p2->p += dpt * p2->tmass / wSum;
        }
        p1->ep -= dpt * p1->tmass / wSum;
        p2->ep += dpt * p2->tmass / wSum;
    } else {
        glm::dvec2 delta = dpt * std::min(kFric * m_d / ldpt, 1.);
        if (m_stabile) {
            p1->p -= delta * p1->tmass / wSum;
            p2->p += delta * p2->tmass / wSum;
        }
        p1->ep -= delta * p1->tmass / wSum;
        p2->ep += delta * p2->tmass / wSum;
    }
}

void RigidContactCS::Draw(const ur::Device& dev, ur::Context& ctx,
                          const std::vector<std::unique_ptr<Particle>>& particles)
{

}

double RigidContactCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    auto& p1 = estimates[m_i1];
    auto& p2 = estimates[m_i2];
    SDFData dat1 = p1->getSDFData(m_bods, m_i1), dat2 = p2->getSDFData(m_bods, m_i2);

    if (dat1.distance < 0 || dat2.distance < 0) {
        glm::dvec2 x12 = p2->GetP(m_stabile) - p1->GetP(m_stabile);
        double len = glm::length(x12);
        m_d = PARTICLE_DIAM - len;
        m_n = len > EPSILON ? -x12 / len : glm::dvec2(0,1);
    } else {
        if (dat1.distance < dat2.distance) {
            m_d = dat1.distance;
            m_n = dat1.gradient;
        } else {
            m_d = dat2.distance;
            m_n = -dat2.gradient;
        }

        if (m_d < PARTICLE_DIAM + EPSILON) {
            InitBoundary(p1, p2);
        }
    }

    return m_d;
}

glm::dvec2 RigidContactCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    if (respect == m_i1) {
        return -m_n;
    }

    if (respect == m_i2) {
        return m_n;
    }

    return glm::dvec2();
}

void RigidContactCS::UpdateCounts(std::vector<int>& counts)
{
    counts[m_i1]++;
    counts[m_i2]++;
}

bool RigidContactCS::InitBoundary(const std::unique_ptr<Particle>& p1, const std::unique_ptr<Particle>& p2)
{
    glm::dvec2 x12 = p1->GetP(m_stabile) - p2->GetP(m_stabile);
    double len = glm::length(x12);
    m_d = PARTICLE_DIAM - len;
    if (m_d < EPSILON) return true;
    x12 = len > EPSILON ? x12 / len : glm::dvec2(0,1);
    double dp = glm::dot(x12, m_n);
    if (dp < 0) {
        m_n = x12 - 2.0 * dp * m_n;
    } else {
        m_n = x12;
    }
    return false;
}

}
}