#include "pbd/Simulation.h"
#include "pbd/OpenSmokeEmitter.h"
#include "pbd/FluidEmitter.h"
#include "pbd/Body.h"
#include "pbd/Particle.h"
#include "pbd/config.h"

#include "pbd/constraint/DistanceCS.h"
#include "pbd/constraint/TotalShapeCS.h"
#include "pbd/constraint/BoundaryCS.h"
#include "pbd/constraint/ContactCS.h"
#include "pbd/constraint/RigidContactCS.h"
#include "pbd/constraint/TotalShapeCS.h"
#include "pbd/constraint/GasCS.h"
#include "pbd/constraint/TotalFluidCS.h"

#include <SM_Calc.h>
#include <tessellation/Painter.h>
#include <unirender/RenderState.h>
#include <painting2/RenderSystem.h>

#include <iostream>

// Iterative or matrix solve
#define ITERATIVE

// Use stabilization pass or not, and if so how many iterations
//#define USE_STABILIZATION
//#define STABILIZATION_ITERATIONS 2

namespace
{

// Number of solver iterations per timestep
const int SOLVER_ITERATIONS = 3;

// Gravity scaling factor for gases
const double ALPHA = -.2;

}

namespace pbd
{

Simulation::Simulation()
{
    //Init(WRECKING_BALL);
    Init(DEBUG_SCENE);
//    m_debug = true;
}

Simulation::~Simulation()
{
}

void Simulation::Clear()
{
    m_particles.clear();
    m_smoke_emitters.clear();
    m_fluid_emitters.clear();
    m_bodies.clear();
//    for (int i = 0; i < NUM_CONSTRAINT_GROUPS; i++) {
//        if(m_globalConstraints.find((ConstraintGroup) i) != m_globalConstraints.end()) {
//            auto& group = m_globalConstraints[(ConstraintGroup) i];
//            for (size_t j = group.size()-1; j >=0; j--) {
//                auto& c = group[j];
//                for (int k = 0; k < NUM_CONSTRAINT_GROUPS; k++) {
//                    if(m_globalConstraints.find((ConstraintGroup) k) != m_globalConstraints.end()) {
////                        m_globalConstraints[(ConstraintGroup) k].removeAll(c);
//                        for (auto& itr = m_globalConstraints[(ConstraintGroup)k].begin(); itr != m_globalConstraints[(ConstraintGroup)k].end(); ) {
//                            if (*itr == c) {
//                                itr = m_globalConstraints[(ConstraintGroup)k].erase(itr);
//                            } else {
//                                ++itr;
//                            }
//                        }
//                    }
//                }
//                delete(c);
//            }
//        }
//    }
    m_global_constraints.clear();

    m_counts.clear();
}

void Simulation::Init(SimulationType type)
{
    this->Clear();

    // Default gravity value
    m_gravity = glm::dvec2(0,-9.8);

    m_point = glm::dvec2(0, 0);

    switch (type) {
    case FRICTION_TEST:
        InitFriction(); break;
    case SDF_TEST:
        InitSdf(); break;
    case GRANULAR_TEST:
        InitGranular(); break;
    case STACKS_TEST:
        InitBoxes(); break;
    case WALL_TEST:
        InitWall(); break;
    case PENDULUM_TEST:
        InitPendulum(); break;
    case ROPE_TEST:
        InitRope(); break;
    case FLUID_TEST:
        InitFluid(); break;
    case FLUID_SOLID_TEST:
        InitFluidSolid(); break;
    case GAS_ROPE_TEST:
        InitRopeGas(); break;
    case WATER_BALLOON_TEST:
        InitWaterBalloon(); break;
    case CRADLE_TEST:
        InitNewtonsCradle(); break;
    case SMOKE_OPEN_TEST:
        InitSmokeOpen(); break;
    case SMOKE_CLOSED_TEST:
        InitSmokeClosed(); break;
    case VOLCANO_TEST:
        InitVolcano(); break;
    case WRECKING_BALL:
        InitWreckingBall(); break;
    case DEBUG_SCENE:
        InitDebugScene(); break;
        break;
    default:
        InitBoxes(); break;
    }

    // Set up the M^-1 matrix
    m_standard_solver.SetupM(m_particles);

    m_counts.resize(m_particles.size(), 0);
}

// (#) in the main simulation loop refer to lines from the main loop in the paper
void Simulation::Tick(double seconds)
{
    std::unordered_map<ConstraintGroup, std::vector<std::shared_ptr<Constraint>>> constraints;

    // Add all rigid body shape constraints
    for (size_t i = 0; i < m_bodies.size(); i++) {
        auto& b = m_bodies[i];
        if (auto c = std::dynamic_pointer_cast<constraint::TotalShapeCS>(b->shape)) {
            constraints[SHAPE].push_back(c);
        } else {
            std::cout << "Rigid body's attached constraint was not a shape constraint." << std::endl;
            exit(1);
        }
    }

    // Add all other global constraints
    for (size_t i = 0; i < m_global_constraints.size(); i++) {
        auto& group = m_global_constraints[(ConstraintGroup) i];
        for (size_t j = 0; j < group.size(); j++) {
            constraints[(ConstraintGroup) i].push_back(group[j]);
        }
    }

    // (1) For all particles
    for (size_t i = 0; i < m_particles.size(); i++) {
        auto& p = m_particles[i];

        // (2) Apply forces
        glm::dvec2 myGravity = m_gravity;
        if(p->ph == GAS) myGravity *= ALPHA;
//        for(OpenSmokeEmitter *e: m_emitters) {
//            for(auto& p : m_particles) {
//                if(glm::distance(p->p, e->GetPosn()) < 1) {
//                    p->f += glm::dvec2(0,.03);
//                }
//            }
//        }
        p->v = p->v + seconds * myGravity + seconds * p->f;
        p->f = glm::dvec2();

        // (3) Predict positions, reset n
        p->ep = p->Guess(seconds);
        m_counts[i] = 0;

        // (4) Apply mass scaling (used by certain constraints)
        p->ScaleMass();
    }
    // (5) End for

    m_contact_solver.SetupM(m_particles, true);

    // (6) For all particles
    for (size_t i = 0; i < m_particles.size(); i++) {
        auto& p = m_particles[i];

        // (7) Find neighboring particles and solid contacts, naive solution
        for (size_t j = i + 1; j < m_particles.size(); j++) {
            auto& p2 = m_particles[j];

            // Skip collision between two immovables
            if (p->imass == 0 && p2->imass == 0) {
                continue;

            // Skip collisions betwee particles in the same rigid body
            } else if (p->ph == SOLID && p2->ph == SOLID && p->bod == p2->bod && p->bod != -1) {
                continue;
            } else {

                // Collision happens when circles overlap
                double dist = glm::distance(p->ep, p2->ep);
                if (dist < PARTICLE_DIAM - EPSILON) {

                    // Rigid contact constraints (which include friction) apply to solid-solid contact
                    if (p->ph == SOLID && p2->ph == SOLID) {
                        constraints[CONTACT].push_back(
                            std::make_shared<constraint::RigidContactCS>(i, j, m_bodies)
                        );
#ifdef USE_STABILIZATION
                        constraints[STABILIZATION].push_back(
                            std::make_shared<constraint::RigidContactCS>(i, j, m_bodies, true)
                        );
#endif
                    // Regular contact constraints (which have no friction) apply to other solid-other contact
                    } else if (p->ph == SOLID || p2->ph == SOLID) {
                        constraints[CONTACT].push_back(
                            std::make_shared<constraint::ContactCS>(i, j)
                        );
                    }
                }
            }
        }

        // (8) Find solid boundary contacts
        if (p->ep.x < m_x_boundaries.x + PARTICLE_RAD) {
            constraints[CONTACT].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_x_boundaries.x, true, true)
            );
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_x_boundaries.x, true, true, true)
            );
#endif
        } else if (p->ep.x > m_x_boundaries.y - PARTICLE_RAD) {
            constraints[CONTACT].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_x_boundaries.y, true, false)
            );
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_x_boundaries.y, true, false, true)
            );
#endif
        }

        if (p->ep.y < m_y_boundaries.x + PARTICLE_RAD) {
            constraints[CONTACT].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_y_boundaries.x, false, true)
            );
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_y_boundaries.x, false, true, true)
            );
#endif
        } else if (p->ep.y > m_y_boundaries.y - PARTICLE_RAD) {
            constraints[CONTACT].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_y_boundaries.y, false, false));
#ifdef USE_STABILIZATION
            constraints[STABILIZATION].push_back(
                std::make_shared<constraint::BoundaryCS>(i, m_y_boundaries.y, false, false, true));
#endif
        }
    }
    // (9) End for

    m_contact_solver.SetupSizes(m_particles.size(), constraints[STABILIZATION]);

#ifdef ITERATIVE

    // (17) For constraint group
    for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
        ConstraintGroup g = (ConstraintGroup) j;

        // Skip the stabilization constraints
        if (g == STABILIZATION) {
            continue;
        }

        //  (18, 19, 20) Update n based on constraints in g
        for (size_t k = 0; k < constraints[g].size(); k++) {
            constraints[g][k]->UpdateCounts(m_counts);
        }
    }

#endif

#ifdef USE_STABILIZATION

    // (10) For stabilization iterations
    for (int i = 0; i < STABILIZATION_ITERATIONS; i++) {

#ifdef ITERATIVE
        // (11, 12, 13, 14) Solve contact constraints and update p, ep, and n
        for (int k = 0; k < constraints[STABILIZATION].size(); k++) {
            constraints[STABILIZATION][k]->Project(m_particles, m_counts);
        }
#else
        // (11, 12, 13, 14) Solve contact constraints and update p, ep, and n
        if (constraints[STABILIZATION].size() > 0) {
            m_contact_solver.SolveAndUpdate(m_particles, constraints[STABILIZATION], true);
        } else {
            break;
        }
#endif

    }
    // (15) End for

#endif

#ifdef ITERATIVE

    // (16) For solver iterations
    for (int i = 0; i < SOLVER_ITERATIONS; i++) {

        // (17) For constraint group
        for (int j = 0; j < (int) NUM_CONSTRAINT_GROUPS; j++) {
            ConstraintGroup g = (ConstraintGroup) j;

            // Skip the stabilization constraints
            if (g == STABILIZATION) {
                continue;
            }

            //  (18, 19, 20) Solve constraints in g and update ep
            for (size_t k = 0; k < constraints[g].size(); k++) {
                constraints[g][k]->Project(m_particles, m_counts);
            }
        }
    }

#else

    m_standard_solver.SetupSizes(m_particles.size(), constraints[STANDARD]);
    m_contact_solver.SetupSizes(m_particles.size(), constraints[CONTACT]);

    // (16) For solver iterations
    for (int i = 0; i < SOLVER_ITERATIONS; i++) {

        // (17, 18, 19, 20) for constraint group, solve constraints and update ep
        if (constraints[CONTACT].size() > 0) {
            m_contact_solver.SolveAndUpdate(m_particles, constraints[CONTACT]);
        }

        if (constraints[STANDARD].size() > 0) {
            m_standard_solver.SolveAndUpdate(m_particles, constraints[STANDARD]);
        }

        if (constraints[SHAPE].size() > 0) {
            for (size_t j = 0; j < constraints[SHAPE].size(); j++) {
                constraints[SHAPE][j]->Project(m_particles, m_counts);
            }
        }
        // (21) End for
    }
    // (22) End for
#endif

    // (23) For all particles
    for (size_t i = 0; i < m_particles.size(); i++) {
        auto& p = m_particles[i];

        // (24) Update velocities
        p->v = (p->ep - p->p) / seconds;

        // (25, 26) Advect diffuse particles, apply internal forces
        /// TODO

        // (27) Update positions or apply sleeping
        p->ConfirmGuess();
    }
    // (28) End for

    // Delete temporary conact constraints
    constraints[CONTACT].clear();
    constraints[STABILIZATION].clear();

    for(auto& e : m_smoke_emitters) {
        e->Tick(m_particles, seconds);
        // (8) Find solid boundary contacts
        for(auto& p : e->GetParticles()) {
            if (p->p.x < m_x_boundaries.x) {
                p->p.x = m_x_boundaries.x;
            } else if (p->p.x > m_x_boundaries.y) {
                p->p.x = m_x_boundaries.y;
            }
            if (p->p.y < m_y_boundaries.x) {
                p->p.y = m_y_boundaries.x;
            } else if (p->p.y > m_y_boundaries.y) {
                p->p.y = m_y_boundaries.y;
            }
        }
    }
    for(auto& e: m_fluid_emitters) {
        e->Tick(m_particles, seconds);
    }
    m_counts.resize(m_particles.size(), 0);
}

void Simulation::CreateRigidBody(std::vector<std::unique_ptr<Particle>>& verts,
                                 const std::vector<SDFData>& sdfData)
{
    if(verts.size() <= 1) {
        std::cout << "Rigid bodies must be at least 2 points." << std::endl;
        exit(1);
    }

    // Compute the total mass, add all the particles to the system and the body
    auto body = std::make_shared<Body>();
    int offset = m_particles.size(), bodyIdx = m_bodies.size();
    double totalMass = 0.0;
    for (size_t i = 0; i < verts.size(); i++) {
        auto& p = verts[i];
        p->bod = bodyIdx;
        p->ph = SOLID;

        if (p->imass == 0.0) {
            std::cout << "A rigid body cannot have a point of infinite mass." << std::endl;
            exit(1);
        }

        totalMass += (1.0 / p->imass);

        m_particles.push_back(std::move(p));
        body->particles.push_back(i + offset);
        body->sdf[i + offset] = sdfData[i];
    }

    // Update the body's global properties, including initial r_i vectors
    body->imass = 1.0 / totalMass;
    body->UpdateCOM(m_particles, false);
    body->ComputeRs(m_particles);
    body->shape = std::make_shared<constraint::TotalShapeCS>(body);

    m_bodies.push_back(body);

    verts.clear();
}

std::shared_ptr<constraint::GasCS>
Simulation::CreateGas(std::vector<std::unique_ptr<Particle>>& verts, double density, bool open=false)
{
    int offset = m_particles.size();
    int bod = static_cast<int>(100 * frand());
    std::vector<int> indices;
    for (size_t i = 0; i < verts.size(); i++) {
        auto& p = verts[i];
        p->ph = GAS;
        p->bod = bod;

        if (p->imass == 0.0) {
            std::cout << "A fluid cannot have a point of infinite mass." << std::endl;
            exit(1);
        }

        m_particles.push_back(std::move(p));
        indices.push_back(offset + i);
    }
    auto gs = std::make_shared<constraint::GasCS>(density, indices, open);
    m_global_constraints[STANDARD].push_back(gs);

    verts.clear();

    return gs;
}

std::shared_ptr<constraint::TotalFluidCS>
Simulation::CreateFluid(std::vector<std::unique_ptr<Particle>>& verts, double density)
{
    int offset = m_particles.size();
    int bod = static_cast<int>(100 * frand());
    std::vector<int> indices;
    for (size_t i = 0; i < verts.size(); i++) {
        auto& p = verts[i];
        p->ph = FLUID;
        p->bod = bod;

        if (p->imass == 0.0) {
            std::cout << "A fluid cannot have a point of infinite mass." << std::endl;
            exit(1);
        }

        m_particles.push_back(std::move(p));
        indices.push_back(offset + i);
    }
    auto fs = std::make_shared<constraint::TotalFluidCS>(density, indices);
    m_global_constraints[STANDARD].push_back(fs);

    verts.clear();

    return fs;
}

void Simulation::CreateSmokeEmitter(glm::dvec2 posn, double particlesPerSec,
                                    const std::shared_ptr<constraint::GasCS>& gs)
{
    m_smoke_emitters.push_back(
        std::make_unique<OpenSmokeEmitter>(posn, particlesPerSec, gs)
    );
}

void Simulation::CreateFluidEmitter(glm::dvec2 posn, double particlesPerSec,
                                    const std::shared_ptr<constraint::TotalFluidCS>& fs)
{
    m_fluid_emitters.push_back(
        std::make_unique<FluidEmitter>(posn, particlesPerSec, fs)
    );
}

void Simulation::Draw(const ur::Device& dev, ur::Context& ctx)
{
    tess::Painter pt;
    pt.SetAntiAliased(false);

    sm::mat4 mt = sm::mat4::Scaled(RENDER_SCALE, RENDER_SCALE, RENDER_SCALE);

    DrawGrid(pt, mt);
    if (m_debug) {
        DrawParticles(pt, mt);
    }
    DrawBodies(dev, ctx, pt, mt);
    DrawGlobals(dev, ctx, pt, mt);
    DrawSmoke(pt, mt);

    //const float half_point_size = 2.5f;
    //pt.AddRectFilled(sm::vec2(m_point.x, m_point.y), half_point_size, 0xffffffff);

    ur::RenderState rs;
    pt2::RenderSystem::DrawPainter(dev, ctx, rs, pt);
}

void Simulation::Resize(const glm::ivec2 &dim)
{
    m_dimensions = dim;
}

void Simulation::DrawGrid(tess::Painter& pt, const sm::mat4& mt)
{
    const float dx = static_cast<float>(m_dimensions.x);
    const float dy = static_cast<float>(m_dimensions.y);

    for (int x = -m_dimensions.x; x <= m_dimensions.x; x++)
    {
        const float px = static_cast<float>(x);
        pt.AddLine(mt * sm::vec2(px, -dy), mt * sm::vec2(px, dy), 0xff333333);
    }
    for (int y = -m_dimensions.y; y <= m_dimensions.y; y++)
    {
        const float py = static_cast<float>(y);
        pt.AddLine(mt * sm::vec2(-dy, py), mt * sm::vec2(dy, py), 0xff333333);
    }

    pt.AddLine(mt * sm::vec2(-dx, 0), mt * sm::vec2(dx, 0), 0xffffffff);
    pt.AddLine(mt * sm::vec2(0, -dy), mt * sm::vec2(0, dy), 0xffffffff);

    const float line_width = 3;
    auto bmin = mt * sm::vec2(static_cast<float>(m_x_boundaries.x), static_cast<float>(m_y_boundaries.x));
    auto bmax = mt * sm::vec2(static_cast<float>(m_x_boundaries.y), static_cast<float>(m_y_boundaries.y));
    pt.AddLine(sm::vec2(bmin.x, bmin.y), sm::vec2(bmin.x, bmax.y), 0xffffffff, line_width);
    pt.AddLine(sm::vec2(bmax.x, bmin.y), sm::vec2(bmax.x, bmax.y), 0xffffffff, line_width);
    pt.AddLine(sm::vec2(bmin.x, bmin.y), sm::vec2(bmax.x, bmin.y), 0xffffffff, line_width);
    pt.AddLine(sm::vec2(bmin.x, bmax.y), sm::vec2(bmax.x, bmax.y), 0xffffffff, line_width);
}

void Simulation::DrawParticles(tess::Painter& pt, const sm::mat4& mt)
{
    for (size_t i = 0; i < m_particles.size(); i++)
    {
        const auto& p = m_particles[i];

        uint32_t color;
        if (p->imass == 0.f) {
            color = 0xff0000ff;
        } else if (p->ph == FLUID || p->ph == GAS){
            color =
                0xff000000 |
                static_cast<uint32_t>((1 - p->bod / 100.) * 255.0) << 16 |
                static_cast<uint32_t>(p->bod / 100. * 255.0) << 8;
        } else if (p->ph == SOLID) {
            color = GetColor(p->bod, 1);
        } else {
            color = 0xffff0000;
        }

        const float radius = 16.0f * PARTICLE_RAD;
        auto center = mt * sm::vec2(static_cast<float>(p->p.x), static_cast<float>(p->p.y));
        pt.AddCircleFilled(center, radius, color);
    }
}

void Simulation::DrawBodies(const ur::Device& dev, ur::Context& ctx,
                            tess::Painter& pt, const sm::mat4& mt)
{
    for (size_t i = 0; i < m_bodies.size(); i++)
    {
        auto& b = m_bodies[i];
        if (m_debug)
        {
            b->shape->Draw(dev, ctx, m_particles);
        }
        else
        {
            for (size_t i = 0; i < b->particles.size(); i++)
            {
                auto& p = m_particles[(b->particles[i])];

                sm::vec2 rect[4];

                const float r = static_cast<float>(PARTICLE_RAD);
                const sm::vec2 c = mt * sm::vec2(static_cast<float>(p->p.x), static_cast<float>(p->p.y));
                const float a = static_cast<float>(b->angle);
                rect[0] = mt * sm::rotate_vector(sm::vec2(-r, -r), a) + c;
                rect[1] = mt * sm::rotate_vector(sm::vec2(-r,  r), a) + c;
                rect[2] = mt * sm::rotate_vector(sm::vec2( r,  r), a) + c;
                rect[3] = mt * sm::rotate_vector(sm::vec2( r, -r), a) + c;

                pt.AddPolygonFilled(&rect[0], 4, GetColor(p->bod, 0.6f));
//                pt.AddPolyline(&rect[0], 4, 0xff000000);
            }
        }
    }
}

void Simulation::DrawGlobals(const ur::Device& dev, ur::Context& ctx,
                             tess::Painter& pt, const sm::mat4& mt)
{
    for (size_t i = 0; i < m_global_constraints.size(); i++) {
        for (size_t j = 0; j < m_global_constraints[(ConstraintGroup) i].size(); j++) {
            m_global_constraints[(ConstraintGroup)i][j]->Draw(dev, ctx, m_particles);
        }
    }
}

void Simulation::DrawSmoke(tess::Painter& pt, const sm::mat4& mt)
{
    float rad = static_cast<float>(PARTICLE_RAD / 7.0);
    for (size_t i = 0; i < m_smoke_emitters.size(); i++)
    {
        auto& particles = m_smoke_emitters[i]->GetParticles();
        for(size_t j = 0; j < particles.size(); j++)
        {
            auto& p = particles[j];
            auto center = mt * sm::vec2(static_cast<float>(p->p.x), static_cast<float>(p->p.y));
            pt.AddRectFilled(center, rad, 0xffffffff);
//            glPushMatrix();
//            glTranslatef(p->p.x, p->p.y, 0);
//            glScalef(PARTICLE_RAD/7., PARTICLE_RAD/7., 0);
//            DrawCircle();
//            glPopMatrix();
        }
    }
}

uint32_t Simulation::GetColor(int body, float alpha)
{
//    uint32_t rgb = 0;
    uint32_t bgr = 0;
    uint32_t a = static_cast<uint32_t>(alpha * 255) << 24;

    int choice = abs(body) % 5;
    if (choice == 0) {
        bgr = 0x00b3ff;
    } else if (choice == 1) {
        bgr = 0xf3c05a;
    } else if (choice == 2) {
        bgr = 0xcc8dff;
    } else if (choice == 3) {
        bgr = 0xe6d91a;
    } else {
        bgr = 0x99e61a;
    }

    return bgr + a;
}

void Simulation::InitFriction()
{
    m_x_boundaries = glm::dvec2(-20,20);
    m_y_boundaries = glm::dvec2(0,1000000);

    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    glm::ivec2 dim = glm::ivec2(3,2);
    for (int x = 0; x < dim.x; x++) {
        double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
        for (int y = 0; y < dim.y; y++) {
            double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
            auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), (x == 0 && y == 0 ? 1 : 1.));
            part->v.x = 5;
            part->kFriction = .01;
            part->sFriction = .1;
            vertices.push_back(std::move(part));
        }
    }
    CreateRigidBody(vertices, data);
}

void Simulation::InitGranular()
{
    m_x_boundaries = glm::dvec2(-100,100);
    m_y_boundaries = glm::dvec2(-5, 1000);
    m_gravity = glm::dvec2(0,-9.8);

    for (int i = -15; i <= 15; i++) {
        for (int j = 0; j < 30; j++) {
            glm::dvec2 pos = glm::dvec2(i * (PARTICLE_DIAM + EPSILON), pow(j,1.2) * (PARTICLE_DIAM) + PARTICLE_RAD + m_y_boundaries.x);
            auto part= std::make_unique<Particle>(pos, 1, SOLID);
            part->sFriction = .35;
            part->kFriction = .3;
            m_particles.push_back(std::move(part));
        }
    }

    auto jerk = std::make_unique<Particle>(glm::dvec2(-25.55, 40), 100.f, SOLID);
    jerk->v.x = 8.5;
    m_particles.push_back(std::move(jerk));
}

void Simulation::InitSdf()
{
    m_x_boundaries = glm::dvec2(-20,20);
    m_y_boundaries = glm::dvec2(0,1000000);

    int numBoxes = 2;
    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,0)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,0)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    glm::ivec2 dim = glm::ivec2(2,3);
    for (int i = numBoxes - 1; i >= 0; i--) {
        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2) + i * PARTICLE_RAD;
            for (int y = 0; y < dim.y; y++) {
                double yVal = ((40 * i) * dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), 4.);
                if (i > 0) part->v.y = -120;
                vertices.push_back(std::move(part));
            }
        }
        CreateRigidBody(vertices, data);
    }
}

void Simulation::InitBoxes()
{
    m_x_boundaries = glm::dvec2(-20,20);
    m_y_boundaries = glm::dvec2(0,1000000);

    int numBoxes = 25, numColumns = 2;
    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -numColumns; j <= numColumns; j++) {
        glm::ivec2 dim = glm::ivec2(3,2);
        for (int i = numBoxes - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double xVal = j * 4 + PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
                for (int y = 0; y < dim.y; y++) {
                    double yVal = ((2 * i + 1) * dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                    auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), 4.);
                    part->sFriction = 1.;
                    part->kFriction = 1.;
                    vertices.push_back(std::move(part));
                }
            }
            CreateRigidBody(vertices, data);
        }
    }
}

void Simulation::InitWall()
{
    m_x_boundaries = glm::dvec2(-50,50);
    m_y_boundaries = glm::dvec2(0,1000000);

    glm::dvec2 dim = glm::dvec2(6,2);
    int height = 11, width = 5;
    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));

    for (int i = 0; i < dim.x - 2; i++) {
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    }

    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -width; j <= width; j++) {
        for (int i = height - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double num = (i % 2 == 0 ? 3 : -1);
                double xVal = j * (EPSILON + dim.x / 2.) + PARTICLE_DIAM * (x % (int)dim.x) - num * PARTICLE_RAD;
                for (int y = 0; y < dim.y; y++) {
                    double yVal = (i * dim.y + (y % (int)dim.y) + EPSILON) * PARTICLE_DIAM + PARTICLE_RAD;
                    auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), 1.);
                    part->sFriction = 1;
                    part->kFriction = 0;
                    vertices.push_back(std::move(part));
                }
            }
            CreateRigidBody(vertices, data);
        }
    }
}

void Simulation::InitPendulum()
{
    m_x_boundaries = glm::dvec2(-10,10);
    m_y_boundaries = glm::dvec2(0,1000000);

    int chainLength = 3;
    m_particles.push_back(std::make_unique<Particle>(glm::dvec2(0, chainLength * 3 + 6) * PARTICLE_DIAM + glm::dvec2(0,2), 0, SOLID));

    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD));

    std::vector<std::unique_ptr<Particle>> vertices;
    double xs[6] = {-1,-1,0,0,1,1};

    for (int i = chainLength; i >= 0; i--) {
        for (int j = 0; j < 6; j++) {
            double y = ((i + 1) * 3 + (j % 2)) * PARTICLE_DIAM + 2;
            auto part = std::make_unique<Particle>(glm::dvec2(xs[j] * PARTICLE_DIAM, y), 1.);
            part->v.x = 3;
            vertices.push_back(std::move(part));
        }
        CreateRigidBody(vertices, data);

        if (i < chainLength) {
            int basePrev = 1 + (chainLength - i - 1) * 6, baseCur = basePrev + 6;
            m_global_constraints[STANDARD].push_back(
                std::make_shared<constraint::DistanceCS>(baseCur + 1, basePrev, m_particles)
            );
            m_global_constraints[STANDARD].push_back(
                std::make_shared<constraint::DistanceCS>(baseCur + 5, basePrev + 4, m_particles)
            );
        }
    }

    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(0, 4, m_particles)
    );
}

void Simulation::InitRope()
{
    double scale = 5.;
    m_x_boundaries = glm::dvec2(-scale,scale);
    m_y_boundaries = glm::dvec2(0,1000000);

    double top = 6, dist = PARTICLE_RAD;

    auto e1 = std::make_unique<Particle>(glm::dvec2(m_x_boundaries.x, top), 0, SOLID);
    e1->bod = -2;
    m_particles.push_back(std::move(e1));

    for (double i = m_x_boundaries.x + dist; i < m_x_boundaries.y - dist; i += dist) {
        auto part = std::make_unique<Particle>(glm::dvec2(i, top), 1., SOLID);
        part->bod = -2;
        m_particles.push_back(std::move(part));
        m_global_constraints[STANDARD].push_back(
            std::make_shared<constraint::DistanceCS>(dist, m_particles.size() - 2, m_particles.size() - 1)
        );
    }

    auto e2 = std::make_unique<Particle>(glm::dvec2(m_x_boundaries.y, top), 0, SOLID);
    e2->bod = -2;
    m_particles.push_back(std::move(e2));

    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(dist, m_particles.size() - 2, m_particles.size() - 1)
    );

    double delta = .7;
    std::vector<std::unique_ptr<Particle>> particles;

    for(double x = -scale; x < scale; x += delta) {
        for(double y = 10; y < 10 + scale; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    CreateFluid(particles, 1.75);
}

void Simulation::InitFluid()
{
    double scale = 4., delta = .7;
    m_gravity = glm::dvec2(0,-9.8);
    m_x_boundaries = glm::dvec2(-2 * scale,2 * scale);
    m_y_boundaries = glm::dvec2(-2 * scale, 10 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    double num = 2.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < scale; y += delta) {
                particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        CreateFluid(particles, 1 + .75 * d);
    }
}

void Simulation::InitFluidSolid()
{
    double scale = 3., delta = .7;
    m_gravity = glm::dvec2(0, -9.8);
    m_x_boundaries = glm::dvec2(-2 * scale,2 * scale);
    m_y_boundaries = glm::dvec2(-2 * scale, 100 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    double num = 1.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y + 3) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        CreateFluid(particles, 1. + 1.25 * (d + 1));
    }

    if(true) {
        particles.clear();
        std::vector<SDFData> data;
        double root2 = sqrt(2);
        glm::ivec2 dim = glm::ivec2(5,2);
        data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
        data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
        for (int i = 0; i < dim.x - 2; i++) {
            data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
            data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
        }
        data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
        data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
            for (int y = 0; y < dim.y; y++) {
                double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                particles.push_back(std::make_unique<Particle>(glm::dvec2(xVal-3, yVal + 10), 2));
            }
        }
        CreateRigidBody(particles, data);
    }

    if(true) {
        particles.clear();
        std::vector<SDFData> data;
        double root2 = sqrt(2);
        glm::ivec2 dim = glm::ivec2(5,2);
        data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
        data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));
        for (int i = 0; i < dim.x - 2; i++) {
            data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
            data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
        }
        data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
        data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

        for (int x = 0; x < dim.x; x++) {
            double xVal = PARTICLE_DIAM * ((x % dim.x) - dim.x / 2);
            for (int y = 0; y < dim.y; y++) {
                double yVal = (dim.y + (y % dim.y) + 1) * PARTICLE_DIAM;
                particles.push_back(std::make_unique<Particle>(glm::dvec2(xVal+3, yVal + 10), .2));
            }
        }
        CreateRigidBody(particles, data);
    }
}


void Simulation::InitGas()
{
    double scale = 2., delta = .7;
    m_gravity = glm::dvec2(0, -9.8);
    m_x_boundaries = glm::dvec2(-2  * scale,2 * scale);
    m_y_boundaries = glm::dvec2(-2  * scale, 10 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    double num = 2.;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        CreateGas(particles, .75 + 3*(d));
    }

    scale = 3;
    for (int d = 0; d < num; d++) {
        double start = -2 * scale + 4 * scale * (d / num);
        for(double x = start; x < start + (4 * scale / num); x += delta) {
            for(double y = -2 * scale; y < 2 * scale; y += delta) {
                particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y+10) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
            }
        }
        CreateFluid(particles, 4. + .75 * (d + 1));
    }
}

void Simulation::InitWaterBalloon()
{
    double scale = 10.;
    m_x_boundaries = glm::dvec2(-scale,scale);
    m_y_boundaries = glm::dvec2(-10,1000000);

    double samples = 60, da = 360. / samples;

    for (int i = 0; i < samples; i++) {
        double angle = D2R(i * da);
        auto part = std::make_unique<Particle>(glm::dvec2(sin(angle), cos(angle)) * 3., 1);
        part->bod = -2;
        int idx = m_particles.size();
        m_particles.push_back(std::move(part));

        if (i > 0) {
            m_global_constraints[STANDARD].push_back(
                std::make_shared<constraint::DistanceCS>(idx, idx - 1, m_particles)
            );
        }
    }
    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(0, m_particles.size() - 1, m_particles)
    );
    int idk = m_particles.size();

    for (int i = 0; i < samples; i++) {
        double angle = D2R(i * da);
        auto part = std::make_unique<Particle>(glm::dvec2(sin(angle), cos(angle) + 3) * 3., 1);
        part->bod = -3;
        int idx = m_particles.size();
        m_particles.push_back(std::move(part));

        if (i > 0) {
            m_global_constraints[STANDARD].push_back(
                std::make_shared<constraint::DistanceCS>(idx, idx - 1, m_particles)
            );
        }
    }
    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(idk, m_particles.size() - 1, m_particles)
    );

    double delta = 1.5 * PARTICLE_RAD;
    std::vector<std::unique_ptr<Particle>> particles;

    for(double x = -2; x <= 2; x += delta) {
        for(double y = -2; y <= 2; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    CreateFluid(particles, 1.75);

    for(double x = -2; x <= 2; x += delta) {
        for(double y = -2; y <= 2; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y + 9) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    CreateFluid(particles, 1.75);
}

void Simulation::InitNewtonsCradle()
{
    m_x_boundaries = glm::dvec2(-10,10);
    m_y_boundaries = glm::dvec2(-5,1000000);

    int n = 2;

    for (int i = -n; i <= n; i++) {
        int idx = m_particles.size();
        m_particles.push_back(std::make_unique<Particle>(glm::dvec2(i * PARTICLE_DIAM, 0), 0.f));
        if (i != -n) {
            m_particles.push_back(std::make_unique<Particle>(glm::dvec2(i * PARTICLE_DIAM, -3), 1.f));
        } else {
            m_particles.push_back(std::make_unique<Particle>(glm::dvec2(i * PARTICLE_DIAM - 3, 0), 1.f));
        }
        m_global_constraints[STANDARD].push_back(
            std::make_shared<constraint::DistanceCS>(idx, idx+1, m_particles)
        );
    }
}

void Simulation::InitSmokeOpen()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_x_boundaries = glm::dvec2(-3  * scale,3 * scale);
    m_y_boundaries = glm::dvec2(-2  * scale,100 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    double start = -2 * scale;
    for(double x = start; x < start + (4 * scale); x += delta) {
        for(double y = -2 * scale; y < 2 * scale; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }

    auto gs = CreateGas(particles, 1.5, true);
    CreateSmokeEmitter(glm::dvec2(0,-2*scale+1), 15, gs);
}

void Simulation::InitSmokeClosed()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_x_boundaries = glm::dvec2(-2  * scale,2 * scale);
    m_y_boundaries = glm::dvec2(-2  * scale,2 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    double start = -2 * scale;
    for(double x = start; x < start + (4 * scale); x += delta) {
        for(double y = -2 * scale; y < 2 * scale; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }

    auto gs = CreateGas(particles, 1.5, false);
    CreateSmokeEmitter(glm::dvec2(0,-2*scale+1), 15, NULL);
}

void Simulation::InitRopeGas()
{
    double scale = 2., delta = .63;
    m_gravity = glm::dvec2(0, -9.8);
    m_x_boundaries = glm::dvec2(-4  * scale,4 * scale);
    m_y_boundaries = glm::dvec2(-2  * scale,100 * scale);

    double top = 12, dist = PARTICLE_RAD;

    auto e1 = std::make_unique<Particle>(glm::dvec2(0, top), 0, SOLID);
    e1->bod = -2;
    m_particles.push_back(std::move(e1));

    for (double i = 0 + dist; i < 4*scale - dist; i += dist) {
        auto part = std::make_unique<Particle>(glm::dvec2(i, top), 2, SOLID);
        part->bod = -2;
        m_particles.push_back(std::move(part));
        m_global_constraints[STANDARD].push_back(
            std::make_shared<constraint::DistanceCS>(dist, m_particles.size() - 2, m_particles.size() - 1)
        );
    }

//    auto e2 = std::make_unique<Particle>(glm::dvec2(2*scale, top), 0, SOLID);
//    e2->bod = -2;
//    m_particles.push_back(e2);

    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(dist, m_particles.size() - 2, m_particles.size() - 1)
    );

    std::vector<std::unique_ptr<Particle>> particles;

    double start = -.5 * scale;
    for(double x = start; x < start + (1 * scale); x += delta) {
        for(double y = -.5 * scale; y < .5 * scale; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    auto gs = CreateGas(particles, 1.5, true);
    CreateSmokeEmitter(glm::dvec2(0,0), 15, gs);
}

void Simulation::InitVolcano()
{
    double scale = 10., delta = .2;

    for(double x = 1.; x <= scale; x+=delta) {
        m_particles.push_back(std::make_unique<Particle>(glm::dvec2(-x,scale-x), 0));
        m_particles.push_back(std::make_unique<Particle>(glm::dvec2(x,scale-x), 0));
    }

    m_gravity = glm::dvec2(0,-9.8);
    m_x_boundaries = glm::dvec2(-2 * scale,2 * scale);
    m_y_boundaries = glm::dvec2(0, 10 * scale);
    std::vector<std::unique_ptr<Particle>> particles;

    delta = .8;
    for(double y = 0.; y < scale-1.; y+=delta) {
        for(double x = 0.; x < scale-y-1; x += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1.1));
            particles.push_back(std::make_unique<Particle>(glm::dvec2(-x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1.1));
        }
    }

    auto fs = CreateFluid(particles, 1);
    CreateFluidEmitter(glm::dvec2(0,0), scale*4, fs);

//    double top = scale-.5, dist = PARTICLE_RAD;

//    auto e1 = std::make_unique<Particle>(glm::dvec2(-1-dist, top), 0, SOLID);
//    e1->bod = -2;
//    m_particles.push_back(e1);

//    for (double i = -1; i <= 2; i += dist) {
//        auto part = std::make_unique<Particle>(glm::dvec2(i, top), 1, SOLID);
//        part->bod = -2;
//        m_particles.push_back(part);
//        m_globalConstraints[STANDARD].push_back(
//                    new DistanceCS(dist, m_particles.size() - 2, m_particles.size() - 1));
//    }
}

void Simulation::InitWreckingBall()
{
    m_x_boundaries = glm::dvec2(-15,100);
    m_y_boundaries = glm::dvec2(0,1000000);

    glm::dvec2 dim = glm::dvec2(6,2);
    int height = 8, width = 2;
    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));

    for (int i = 0; i < dim.x - 2; i++) {
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    }

    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -width; j <= width; j++) {
        for (int i = height - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double num = (i % 2 == 0 ? 3 : -1);
                double xVal = j * (EPSILON + dim.x / 2.) + PARTICLE_DIAM * (x % (int)dim.x) - num * PARTICLE_RAD;
                for (int y = 0; y < dim.y; y++) {
                    double yVal = (i * dim.y + (y % (int)dim.y) + EPSILON) * PARTICLE_DIAM + PARTICLE_RAD;
                    auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), 30.);
                    part->sFriction = 1;
                    part->kFriction = 1;
                    vertices.push_back(std::move(part));
                }
            }
            CreateRigidBody(vertices, data);
        }
    }

    double scale = 6., delta = .4;
    std::vector<std::unique_ptr<Particle>> particles;

    double num = 1.;
    double start = m_x_boundaries.x + 1;
    for(double x = start; x < start + (scale / num); x += delta) {
        for(double y = 0; y < 1.2 * scale; y += delta) {
            particles.push_back(std::make_unique<Particle>(glm::dvec2(x,y) + .2 * glm::dvec2(frand() - .5, frand() - .5), 1));
        }
    }
    CreateFluid(particles, 2.5);

    int idx = m_particles.size();
    m_particles.push_back(std::make_unique<Particle>(glm::dvec2(10, 50), 0));
    data.clear();

    glm::dvec2 base = glm::dvec2(57, 50);
    particles.push_back(std::make_unique<Particle>(base, 1000));
    for (double a = 0; a <= 360; a+=30) {
        glm::dvec2 vec = glm::dvec2(cos(D2R(a)), sin(D2R(a)));
        particles.push_back(std::make_unique<Particle>(vec * PARTICLE_RAD + base, 1000));
        data.push_back(SDFData(vec, PARTICLE_RAD * 1.5));
    }
    data.push_back(SDFData());
    CreateRigidBody(particles, data);

    m_global_constraints[STANDARD].push_back(
        std::make_shared<constraint::DistanceCS>(idx, idx + 1, m_particles)
    );
}

void Simulation::InitDebugScene()
{
    m_x_boundaries = glm::dvec2(-5, 5);
    m_y_boundaries = glm::dvec2(0,1000000);

    glm::dvec2 dim = glm::dvec2(6,2);
//    int height = 8, width = 2;
//    int height = 4, width = 1;
    int height = 2, width = 1;
    double root2 = sqrt(2);
    std::vector<std::unique_ptr<Particle>> vertices;
    std::vector<SDFData> data;
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(-1,1)), PARTICLE_RAD * root2));

    for (int i = 0; i < dim.x - 2; i++) {
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,-1)), PARTICLE_RAD));
        data.push_back(SDFData(glm::normalize(glm::dvec2(0,1)), PARTICLE_RAD));
    }

    data.push_back(SDFData(glm::normalize(glm::dvec2(1,-1)), PARTICLE_RAD * root2));
    data.push_back(SDFData(glm::normalize(glm::dvec2(1,1)), PARTICLE_RAD * root2));

    for (int j = -width; j <= width; j++) {
        for (int i = height - 1; i >= 0; i--) {
            for (int x = 0; x < dim.x; x++) {
                double num = (i % 2 == 0 ? 3 : -1);
                double xVal = j * (EPSILON + dim.x / 2.) + PARTICLE_DIAM * (x % (int)dim.x) - num * PARTICLE_RAD;
                for (int y = 0; y < dim.y; y++) {
                    double yVal = (i * dim.y + (y % (int)dim.y) + EPSILON) * PARTICLE_DIAM + PARTICLE_RAD;
                    auto part = std::make_unique<Particle>(glm::dvec2(xVal, yVal), 30.);
                    part->sFriction = 1;
                    part->kFriction = 1;
                    vertices.push_back(std::move(part));
                }
            }
            CreateRigidBody(vertices, data);
        }
        break;
    }
}

int Simulation::GetNumParticles()
{
    return m_particles.size();
}

double Simulation::GetKineticEnergy()
{
    double energy = 0;
    for (size_t i = 0; i < m_particles.size(); i++) {
        auto& p = m_particles[i];
        if (p->imass != 0.) {
            energy += .5 * glm::dot(p->v, p->v) / p->imass;
        }
    }
    return energy;
}

void Simulation::MousePressed(const glm::dvec2 &p)
{
    for (size_t i = 0; i < m_particles.size(); i++) {
        auto& part = m_particles[i];

        glm::dvec2 to = glm::normalize(p - part->p);
        part->v += 7. * to;
    }
    m_point = p;
}

}