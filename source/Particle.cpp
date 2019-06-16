#include "pbd/Particle.h"
#include "pbd/SDFData.h"
#include "pbd/Body.h"

#include <math.h>

namespace pbd
{

SDFData Particle::getSDFData(const std::vector<std::shared_ptr<Body>>& bodies, int idx)
{
    if (ph != SOLID || bod < 0) {
        return SDFData();
    }

    auto& body = bodies[bod];
    SDFData out = body->sdf[idx];
    out.rotate(body->angle);
    return out;
}

}