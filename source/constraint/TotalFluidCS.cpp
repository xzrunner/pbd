#include "pbd/constraint/TotalFluidCS.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include <iostream>

namespace
{

//#define H 3.5
//#define H6 1838.265625
//#define H9 78815.6386719

//#define H 5.
//#define H6 15625.
//#define H9 1953125.

//#define H 4.
//#define H6 4096.
//#define H9 262144.

// Epsilon in gamma correction denominator
const double RELAXATION = .01;

// Pressure terms
const double K_P  = .1;
const double E_P  = 4;
const double DQ_P = .2;

// Fluid-solid coupling constant
const double S_SOLID = 0.;

}

namespace pbd
{
namespace constraint
{

TotalFluidCS::TotalFluidCS(double density, const std::vector<int>& particles)
    : m_p0(density)
{
    m_num_particles = particles.size();

    m_neighbors.resize(m_num_particles);
    for (int i = 0; i < m_num_particles; ++i) {
        m_neighbors[i] = std::make_unique<std::vector<int>>();
    }

    m_deltas.resize(m_num_particles);

    for (int i = 0; i < m_num_particles; i++) {
        m_ps.push_back(particles[i]);
    }
}

void TotalFluidCS::AddParticle(int index)
{
    assert(m_neighbors.size() == m_num_particles);
    m_neighbors.push_back(std::make_unique<std::vector<int>>());

    assert(m_deltas.size() == m_num_particles);
    m_deltas.push_back(glm::dvec2());

    // Clear
    for (int i = 0; i < m_num_particles; ++i) {
        m_neighbors[i]->clear();
        m_deltas[i] = glm::dvec2();
    }

    m_ps.push_back(index);

    m_num_particles++;
}

void TotalFluidCS::RemoveParticle(int index)
{
//    if(ps.contains(index)) {

        assert(m_neighbors.size() == m_num_particles);
        m_neighbors.pop_back();

        assert(m_deltas.size() == m_num_particles);
        m_deltas.pop_back();

        m_num_particles--;

        // Clear
        for (int i = 0; i < m_num_particles; ++i) {
            m_neighbors[i]->clear();
            m_deltas[i] = glm::dvec2();
        }

        m_ps.erase(m_ps.begin() + index);
//    }
}

void TotalFluidCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{
    // Find neighboring particles and estimate pi for each particle
    m_lambdas.clear();
    for (size_t k = 0; k < m_ps.size(); k++)
    {
        m_neighbors[k]->clear();
        int i = m_ps[k];
        auto& p_i = estimates[i];
        double pi = 0., denom = 0.;

        // Find neighbors
        for (size_t j = 0; j < estimates.size(); j++)
        {
            // Check if the next particle is actually this particle
            if (j != i)
            {
                auto& p_j = estimates[j];

                // Ignore fixed particles
                if (p_j->imass == 0) continue;
                glm::dvec2 r = p_i->ep - p_j->ep;
                double rlen2 = glm::dot(r, r);
                if (rlen2 < H2)
                {
                    // Found a neighbor! Remember it and add to pi and the gamma denominator
                    m_neighbors[k]->push_back(j);
                    double incr = Poly6(rlen2) / p_j->imass;
                    if (p_j->ph == SOLID) {
                        incr *= S_SOLID;
                    }
                    pi += incr;

                    glm::dvec2 gr = Grad(estimates, k, j);
                    denom += glm::dot(gr, gr);
                }
            // If it is, cut to the chase
            }
            else
            {
                m_neighbors[k]->push_back(j);
                pi += Poly6(0) / p_i->imass;
            }
        }

        glm::dvec2 gr = Grad(estimates, k, i);
        denom += glm::dot(gr, gr);

        // Compute the gamma value
//        cout << i << " estimated " << pi << endl;
        double lambda = -((pi / m_p0) - 1.) / (denom + RELAXATION);
        m_lambdas[i] = lambda;
    }

    // Compute actual deltas
    for (size_t k = 0; k < m_ps.size(); k++)
    {
        glm::dvec2 delta = glm::dvec2();
        int i = m_ps[k];
        auto& p_i = estimates[i];

        for (size_t x = 0; x < m_neighbors[k]->size(); x++)
        {
            int j = (*m_neighbors[k])[x];
            if (i == j) continue;
            auto& p_j = estimates[j];
            glm::dvec2 r = p_i->ep - p_j->ep;
            double rlen = glm::length(r);
            glm::dvec2 sg = SpikyGrad(r, rlen);
            double lambdaCorr = -K_P * pow((Poly6(rlen * rlen) / Poly6(DQ_P * DQ_P * H * H)), E_P);
            delta += (m_lambdas[i] + m_lambdas[j] + lambdaCorr) * sg;
        }
        m_deltas[k] = (delta / m_p0);
    }

    for (size_t k = 0; k < m_ps.size(); k++)
    {
        int i = m_ps[k];
        auto& p_i = estimates[i];
        p_i->ep += m_deltas[k] / ((double) m_neighbors[k]->size() + counts[i]);
    }
}

void TotalFluidCS::Draw(const ur::Device& dev, ur::Context& ctx,
                        const std::vector<std::unique_ptr<Particle>>& particles)
{
}

double TotalFluidCS::Poly6(double r2)
{
    if(r2 >= H2) return 0;
    double term2 = (H2 - r2);
    return (315. / (64. * M_PI * H9)) * (term2 * term2 * term2);
//    return (H-r) / (H*H);
}

glm::dvec2 TotalFluidCS::SpikyGrad(const glm::dvec2 &r, double rlen2)
{
    if(rlen2 >= H) return glm::dvec2();
    if(rlen2 == 0) return glm::dvec2();
    return -glm::normalize(r) * (45. / (M_PI * H6)) * (H - rlen2) * (H - rlen2);
//    return -r / (H*H*rlen);
}

glm::dvec2 TotalFluidCS::Grad(const std::vector<std::unique_ptr<Particle>>& estimates, int k, int j)
{
    int i = m_ps[k];
    auto& p_i = estimates[i];
    auto& p_j = estimates[j];
    glm::dvec2 r = p_i->ep - p_j->ep;
    double rlen = glm::length(r);
    if (p_i != p_j) {
        return -SpikyGrad(r, rlen) / (m_p0);
    }

    glm::dvec2 out = glm::dvec2();
    for (size_t x = 0; x < m_neighbors[k]->size(); x++)
    {
        auto& p_j = estimates[(*m_neighbors[k])[x]];
        r = p_i->ep - p_j->ep;
        rlen = glm::length(r);
        out += (p_j->ph == SOLID ? S_SOLID : 1.) * SpikyGrad(r, rlen);
    }

    return out / (m_p0);
}

double TotalFluidCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    std::cout << "You shouldn't be calling evaluate on fluids" << std::endl;
    exit(1);
}

glm::dvec2 TotalFluidCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    std::cout << "You shouldn't be calling gradient on fluids" << std::endl;
    exit(1);
}

void TotalFluidCS::UpdateCounts(std::vector<int>& counts)
{
}

}
}