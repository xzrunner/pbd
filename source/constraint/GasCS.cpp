#include "pbd/constraint/GasCS.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include <glm/vec3.hpp>

#include <iostream>

namespace
{

// USE H = 4, density = .5, look for the lattice

// Epsilon in gamma correction denominator
const double RELAXATION = .01;

// Pressure terms
const double K_P  = .2;
const double E_P  = 4;
const double DQ_P = .25;

// Fluid-solid coupling constant
const double S_SOLID = .5;

}

namespace pbd
{
namespace constraint
{

GasCS::GasCS(double density, const std::vector<int>& particles, bool open)
    : Constraint(), m_p0(density), m_open(open)
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

void GasCS::AddParticle(int index)
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

void GasCS::Project(const std::vector<std::unique_ptr<Particle>>& estimates, const std::vector<int>& counts)
{

    // Find neighboring particles and estimate pi for each particle
    m_lambdas.clear();
    for (size_t k = 0; k < m_ps.size(); k++) {
        m_neighbors[k]->clear();
        int i = m_ps[k];
        auto& p_i = estimates[i];
        double pi = 0., denom = 0.;

        // Find neighbors
        for (size_t j = 0; j < estimates.size(); j++) {

            // Check if the next particle is actually this particle
            if (j != i) {
                auto& p_j = estimates[j];

                // Ignore fixed particles
                if (p_j->imass == 0) continue;
                glm::dvec2 r = p_i->ep - p_j->ep;
                double rlen2 = glm::dot(r, r);
                if (rlen2 < H2) {

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
            } else {
                m_neighbors[k]->push_back(j);
                pi += Poly6(0) / p_i->imass;
            }
        }

        glm::dvec2 gr = Grad(estimates, k, i);
        denom += glm::dot(gr, gr);

        // Compute the gamma value
//        cout << i << " estimated " << pi << endl;

        double p_rat = (pi/m_p0);
        if(m_open) p_i->f += p_i->v * (1.-p_rat) * -50.;
//        if(p_rat < 1) p_rat = 1;
        double lambda = -(p_rat - 1.) / (denom + RELAXATION);
        m_lambdas[i] = lambda;
    }

    // Compute actual deltas
    for (size_t k = 0; k < m_ps.size(); k++) {
        glm::dvec2 delta = glm::dvec2();
        glm::dvec2 f_vort = glm::dvec2();
        int i = m_ps[k];
        auto& p_i = estimates[i];

        for (size_t x = 0; x < m_neighbors[k]->size(); x++) {
            int j = (*m_neighbors[k])[x];
            if (i == j) continue;
            auto& p_j = estimates[j];
            glm::dvec2 r = p_i->ep - p_j->ep;
            double rlen = glm::length(r);
            glm::dvec2 sg = SpikyGrad(r, rlen);
            double lambdaCorr = -K_P * pow((Poly6(rlen * rlen) / Poly6(DQ_P * DQ_P * H * H)), E_P);
            delta += (m_lambdas[i] + m_lambdas[j] + lambdaCorr) * sg;
//            vorticity
            glm::dvec2 gradient = SpikyGrad(r, glm::dot(r,r));
            glm::dvec2 w = gradient * p_j->v;
            glm::dvec3 cross = glm::cross(glm::dvec3(0,0,glm::length(w)), glm::dvec3(r.x, r.y, 0));
            f_vort += glm::dvec2(cross.x, cross.y) * Poly6(glm::dot(r,r));
        }
        m_deltas[k] = (delta / m_p0);
        p_i->f += f_vort;
    }

    for (size_t k = 0; k < m_ps.size(); k++) {
        int i = m_ps[k];
        auto& p_i = estimates[i];
        p_i->ep += m_deltas[k] / ((double) m_neighbors[k]->size() + counts[i]);
    }

//    // Find neighboring particles and estimate pi for each particle
//    lambdas.Clear();
//    for (size_t k = 0; k < ps.size(); k++) {
//        neighbors[k]->Clear();
//        int i = ps[k];
//        Particle *p_i = estimates->at(i);
//        double pi = 0., denom = 0.;

//        // Find neighbors and calculate forces
//        for (size_t j = 0; j < estimates->size(); j++) {

//            // Check if the next particle is actually this particle
//            if (j != i) {
//                Particle *p_j = estimates->at(j);
//                glm::dvec2 r = p_i->ep - p_j->ep;
//                double rlen2 = glm::dot(r, r);
//                if (rlen2 < H2) {

//                    // Found a neighbor! Remember it and add to pi and the gamma denominator
//                    neighbors[k]->push_back(j);
//                    double incr = Poly6(rlen2) / p_j->imass;
//                    if (p_j->ph == SOLID) {
//                        incr *= S_SOLID;
//                    }
//                    pi += incr;

//                    glm::dvec2 gr = Grad(estimates, k, j);
//                    denom += glm::dot(gr, gr);
//                }

//            // If it is, cut to the chase
//            } else {
//                neighbors[k]->push_back(j);
//                pi += Poly6(0) / p_i->imass;
//            }
//        }

//        glm::dvec2 gr = Grad(estimates, k, i);
//        denom += glm::dot(gr, gr);

//        // Compute the gamma value
////        cout << i << " estimated " << pi << endl;
//        double p_rat = (pi/p0);
////        if(m_open) p_i->f += p_i->v * (1.-p_rat) * -50.;
////        if(p_rat < 1) p_rat = 1;
//        double lambda = -(p_rat - 1.) / (denom + RELAXATION);
//        lambdas[i] = lambda;
//    }

//    // Compute actual deltas
//    for (size_t k = 0; k < ps.size(); k++) {
//        glm::dvec2 delta = glm::dvec2();
//        glm::dvec2 f_vort = glm::dvec2();
//        int i = ps[k];
//        Particle *p_i = estimates->at(i);

//        for (size_t x = 0; x < neighbors[k]->size(); x++) {
//            int j = neighbors[k][x];
//            if (i == j) continue;
//            Particle *p_j = estimates->at(j);
//            glm::dvec2 r = p_i->ep - p_j->ep;
//            double rlen = glm::length(r);
//            glm::dvec2 sg = SpikyGrad(r, rlen);
//            double lambdaCorr = -K_P * pow((Poly6(rlen) / Poly6(DQ_P * H)), E_P);
//            delta += (lambdas[i] + lambdas[j] + lambdaCorr) * sg;

//            // vorticity
////            glm::dvec2 gradient = SpikyGrad(r, glm::dot(r,r));
////            glm::dvec2 w = gradient * p_j->v;
////            glm::dvec3 cross = glm::cross(glm::dvec3(0,0,glm::length(w)), glm::dvec3(r.x, r.y, 0));
////            f_vort += glm::dvec2(cross.x, cross.y) * Poly6(glm::dot(r,r));
//        }
//        deltas[k] = (delta / p0);
////        p_i->f += f_vort;
//    }

//    for (size_t k = 0; k < ps.size(); k++) {
//        int i = ps[k];
//        Particle *p_i = estimates->at(i);
//        p_i->ep += deltas[k] / ((double) neighbors[k]->size() + counts[i]);
//    }
}

void GasCS::Draw(const std::vector<std::unique_ptr<Particle>>& particles)
{

}

double GasCS::Poly6(double r2)
{
    if(r2 >= H2) return 0;
    double term2 = (H2 - r2);
    return (315. / (64. * M_PI * H9)) * (term2 * term2 * term2);
//    return (H-r) / (H*H);
}

glm::dvec2 GasCS::SpikyGrad(const glm::dvec2 &r, double rlen)
{
    if(rlen >= H) return glm::dvec2();
    if(rlen == 0) return glm::dvec2();
    return -glm::normalize(r) * (45. / (M_PI * H6)) * (H - rlen) * (H - rlen);
//    return -r / (H*H*rlen);
}

glm::dvec2 GasCS::Grad(const std::vector<std::unique_ptr<Particle>>& estimates, int k, int j)
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
    for (size_t x = 0; x < m_neighbors[k]->size(); x++) {
        r = p_i->ep - estimates[(*m_neighbors[k])[x]]->ep;
        rlen = glm::length(r);
        out += SpikyGrad(r, rlen);
    }

    return out / (m_p0);
}

double GasCS::Evaluate(const std::vector<std::unique_ptr<Particle>>& estimates)
{
    std::cout << "You shouldn't be calling evaluate on fluids" << std::endl;
    exit(1);
}

glm::dvec2 GasCS::Gradient(const std::vector<std::unique_ptr<Particle>>& estimates, int respect)
{
    std::cout << "You shouldn't be calling gradient on fluids" << std::endl;
    exit(1);
}

void GasCS::UpdateCounts(std::vector<int>& counts)
{
}

}
}