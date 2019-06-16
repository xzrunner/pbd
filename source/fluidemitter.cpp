#include "pbd/FluidEmitter.h"
#include "pbd/Particle.h"
#include "pbd/constraint/TotalFluidCS.h"

namespace pbd
{

FluidEmitter::FluidEmitter(glm::dvec2 posn, double particlesPerSec,
                           const std::shared_ptr<constraint::TotalFluidCS>& fs)
    : m_posn(posn)
    , m_particlesPerSec(particlesPerSec)
    , m_fs(fs)
{
    m_timer = 0;
    m_tot_timer = 0;
}

void FluidEmitter::Tick(std::vector<std::unique_ptr<Particle>>& estimates, double secs)
{
    for(int i = m_fs->m_ps.size()-1; i >= 0; i--)
    {
        auto& p = estimates[m_fs->m_ps[i]];
            //            std::cout << p << std::endl;
//            double lambda = m_fs->lambdas[i];
            //            std::cout << lambda << std::endl;
        //            if(lambda >= -.1 && glm::length(p->v) < .05 && glm::length(p->p - p->ep) < .05) {
        //            if(p->p.y >= 10 || fabs(p->p.x) >= 10 ) {
        if(glm::length(p->v) < .06 && p->p.y <= 5)
        {
            if(m_fs->m_lambdas[i] <= 0)
            {
                p->t -= 1;
                if(p->t <= 0)
                {
                    p->t = 0;

                    //                p->ph = SOLID;
                    //                Particle *newP = new Particle(p->p, 0, SOLID);
                    //                newP->v = p->v;
                    //            p->imass -= secs;
                    //            if(p->imass == 0)

                    p->imass = 0;
                    p->ph = SOLID;
                    p->ep = p->p;
                    p->v = glm::dvec2();
                    p->f = glm::dvec2();
                    //                m_grains.push_back(newP);
                    //                estimates->append(newP);
                    //                estimates->removeAt(m_fs->ps.at(i));
                    //                if(m_fs->ps.contains(i))
                    m_fs->RemoveParticle(i);
                    //                delete p;
                    //                p->imass = 1;
                    //                p->ph = SOLID;
                    //                m_fs->ps.removeAt(i);
                }
            }
            else
            {
                p->t += secs;
                if (p->t > 3) {
                    p->t = 3;
                }
            }
            //        }
        }
    }


    //    for(int i=0; i<m_grains.size(); i++) {
    //        Particle *p = m_grains.at(i);
//        if(glm::length(p->v) <= .02) {
//            m_grains.removeAt(i);
//            p->imass = 0;
//            p->v = glm::dvec2();
//            p->f = glm::dvec2();
//            p->ep = p->p;
//        }
//    }

    m_timer += secs;
    m_tot_timer += secs;
    while(m_tot_timer < 5 && m_timer >= 1./m_particlesPerSec)
    {
        m_timer -= 1./m_particlesPerSec;
        if(m_fs != NULL)
        {
            auto p = std::make_unique<Particle>(m_posn, 1, FLUID);
            p->v = glm::dvec2(frand(),1);
            m_fs->AddParticle(estimates.size());
            estimates.push_back(std::move(p));
        }
    }
}

}