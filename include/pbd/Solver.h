#pragma once

#include <Eigen/Sparse>

#include <vector>
#include <memory>

namespace pbd
{

class Constraint;
struct Particle;

class Solver
{
public:
    Solver();

    Eigen::SparseMatrix<double> m_invM, m_JT, m_A;
    Eigen::VectorXd m_b, m_gamma, m_dp;
    std::vector<int> m_counts;
    int m_n_parts, m_n_cons;

    int GetCount(int idx);

    void SetupM(const std::vector<std::unique_ptr<Particle>>& particles, bool contact = false);
    void SetupSizes(int numParts, const std::vector<std::shared_ptr<Constraint>>& constraints);
    void SolveAndUpdate(const std::vector<std::unique_ptr<Particle>>& particles, const std::vector<std::shared_ptr<Constraint>>& constraints, bool stabile = false);

}; // Solver

}