#include "pbd/Solver.h"
#include "pbd/Constraint.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

namespace
{

const double RELAXATION_PARAMETER = 1.;

}

namespace pbd
{

Solver::Solver()
{
    m_b = Eigen::VectorXd(2);
    m_gamma = Eigen::VectorXd(2);
    m_n_cons = -1;

    m_dp = Eigen::VectorXd(2);
    m_counts.resize(2, 0);
    m_n_parts = -1;
}

int Solver::GetCount(int idx)
{
    return m_counts[idx];
}

void Solver::SetupM(const std::vector<std::unique_ptr<Particle>>& particles, bool contact)
{
    m_invM.resize(particles.size() * 2, particles.size() * 2);
    for (size_t i = 0; i < particles.size(); i++) {

        // Down the diagonal
        auto& p = particles[i];
        m_invM.insert(2 * i, 2 * i) = contact ? p->tmass : p->imass;
        m_invM.insert(2 * i + 1, 2 * i + 1) = contact ? p->tmass : p->imass;
    }
}

void Solver::SetupSizes(int numParts, const std::vector<std::shared_ptr<Constraint>>& constraints)
{
    bool change = true;
    int numCons = constraints.size();

    // Only update some things if the number of particles changed
    if (m_n_parts != numParts) {
        m_n_parts = numParts;
        m_dp = Eigen::VectorXd(m_n_parts * 2);
        m_counts.resize(m_n_parts, 0);
        change = true;
    }

    // Update how many constraints affect each particle
    for (int i = 0; i < m_n_parts; i++) {
        m_counts[i] = 0;
    }
    for (int i = 0; i < numCons; i++) {
        constraints[i]->UpdateCounts(m_counts);
    }

    if (m_n_cons != numCons) {
        m_n_cons = numCons;
        m_b = Eigen::VectorXd(m_n_cons);
        m_gamma = Eigen::VectorXd(m_n_cons);
        change = true;
    }

    if (change) {
        m_JT.resize(m_n_parts * 2, m_n_cons);
    }
}

void Solver::SolveAndUpdate(const std::vector<std::unique_ptr<Particle>>& particles, const std::vector<std::shared_ptr<Constraint>>& constraints, bool stabile)
{
    if (constraints.size() == 0) {
        return;
    }

    // Reset J^T and b
    bool updated = false;

    // Loop!
    for (size_t i = 0; i < particles.size(); i++) {
        for (size_t j = 0; j < constraints.size(); j++) {
            auto& cons = constraints[j];

            // Update b
            if (!updated) {
                m_b[j] = -cons->Evaluate(particles);
            }
            glm::vec2 grad_ji = cons->Gradient(particles, i);
            m_JT.coeffRef(2 * i, j) = grad_ji.x;
            m_JT.coeffRef(2 * i + 1, j) = grad_ji.y;
        }
        updated = true;
    }

    auto temp = m_invM * m_JT;
    m_A = m_JT.transpose() * temp;
//    m_JT.printMatrix(4,false);

    Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> chol(m_A);  // performs a Cholesky factorization of A
    m_gamma = chol.solve(m_b);

    //m_eq.setA(&m_A);
//    cout << endl;
//    for (size_t i = 0; i < particles.size(); i++) {
//        printf("%.4f\n", m_b[i]);
//    }
//    bool result = m_eq.solve(m_b, m_gamma);
//    cout << result << endl;
//    for (size_t i = 0; i < particles.size(); i++) {
//        printf("%.4f\n", m_gamma[i]);
//    }
//    cout << endl;
    m_dp = temp * m_gamma;

    for (size_t i = 0; i < particles.size(); i++) {
        auto& p = particles[i];
        int n = m_counts[i];
        double mult = n > 0 ? (RELAXATION_PARAMETER / (double)n) : 0.,
               dx = m_dp[2 * i] * mult,
               dy = m_dp[2 * i + 1] * mult;
//        cout << dx << " " << dy << endl;

        p->ep.x += (fabs(dx) > EPSILON ? dx : 0);
        p->ep.y += (fabs(dy) > EPSILON ? dy : 0);

        if (stabile) {
            p->p.x += (fabs(dx) > EPSILON ? dx : 0);
            p->p.y += (fabs(dy) > EPSILON ? dy : 0);
        }
    }
}

}