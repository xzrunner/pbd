#pragma once

#include "pbd/Solver.h"
#include "pbd/Constraint.h"
#include "pbd/SDFData.h"

#include <SM_Matrix.h>

#include <unordered_map>

namespace ur { class Device; class Context; }
namespace tess { class Painter; }

namespace pbd
{

class OpenSmokeEmitter;
class FluidEmitter;
struct Particle;
struct Body;

namespace constraint
{
    class GasCS;
    class TotalFluidCS;
}

// Built-in simulation scenes
enum SimulationType
{
    FRICTION_TEST,
    SDF_TEST,
    GRANULAR_TEST,
    STACKS_TEST,
    WALL_TEST,
    PENDULUM_TEST,
    ROPE_TEST,
    FLUID_TEST,
    FLUID_SOLID_TEST,
    GAS_ROPE_TEST,
    WATER_BALLOON_TEST,
    CRADLE_TEST,
    NUM_SIMULATION_TYPES,
    SMOKE_OPEN_TEST,
    SMOKE_CLOSED_TEST,
    VOLCANO_TEST,
    WRECKING_BALL,
    DEBUG_SCENE,
};

// The basic simulation, implementing the "main solve loop" from the paper.
class Simulation
{
public:
    Simulation();
    ~Simulation();

    void Init(SimulationType type);

    // Initializers for test scenes
    void InitFriction();
    void InitSdf();
    void InitGranular();
    void InitBoxes();
    void InitPendulum();
    void InitRope();
    void InitFluid();
    void InitFluidSolid();
    void InitWall();
    void InitGas();
    void InitWaterBalloon();
    void InitNewtonsCradle();
    void InitSmokeOpen();
    void InitSmokeClosed();
    void InitRopeGas();
    void InitVolcano();
    void InitWreckingBall();
    void InitDebugScene();

    // Basic interaction events
    void Tick(double seconds);
    void Draw(const ur::Device& dev, ur::Context& ctx);
    void Resize(const glm::ivec2 &dim);
    void MousePressed(const glm::dvec2 &p);

    // Debug information and flags
    int GetNumParticles();
    double GetKineticEnergy();

private:
    // Reset the simulation
    void Clear();

    // Creation functions for different types of matter
    void CreateRigidBody(std::vector<std::unique_ptr<Particle>>& verts,
        const std::vector<SDFData>& sdfData);
    std::shared_ptr<constraint::TotalFluidCS>
        CreateFluid(std::vector<std::unique_ptr<Particle>>& particles, double density);
    std::shared_ptr<constraint::GasCS>
        CreateGas(std::vector<std::unique_ptr<Particle>>& particles, double density, bool open);
    void CreateSmokeEmitter(glm::dvec2 posn, double particlesPerSec,
        const std::shared_ptr<constraint::GasCS>& gs);
    void CreateFluidEmitter(glm::dvec2 posn, double particlesPerSec,
        const std::shared_ptr<constraint::TotalFluidCS>& fs);

    // Simple drawing routines
    void DrawGrid(tess::Painter& pt, const sm::mat4& mt);
    void DrawParticles(tess::Painter& pt, const sm::mat4& mt);
    void DrawBodies(const ur::Device& dev, ur::Context& ctx,
        tess::Painter& pt, const sm::mat4& mt);
    void DrawGlobals(const ur::Device& dev, ur::Context& ctx,
        tess::Painter& pt, const sm::mat4& mt);
    void DrawSmoke(tess::Painter& pt, const sm::mat4& mt);

    uint32_t GetColor(int body, float alpha);

public:
    bool m_debug = false;

private:
    // Counts for iterative particle solver
    std::vector<int> m_counts;

    // Storage of global particles, rigid bodies, and general constraints
    std::vector<std::unique_ptr<Particle>> m_particles;
    std::vector<std::shared_ptr<Body>>     m_bodies;
    std::vector<std::unique_ptr<OpenSmokeEmitter>> m_smoke_emitters;
    std::vector<std::unique_ptr<FluidEmitter>>     m_fluid_emitters;
    std::unordered_map<ConstraintGroup, std::vector<std::shared_ptr<Constraint>>> m_global_constraints;

    // Solvers for regular and contact constraints
    Solver m_standard_solver;
    Solver m_contact_solver;

    // Drawing and boundary information
    glm::ivec2 m_dimensions;
    glm::dvec2 m_x_boundaries, m_y_boundaries;
    glm::dvec2 m_gravity;
    glm::dvec2 m_point;
};

}