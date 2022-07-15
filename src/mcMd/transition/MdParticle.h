#ifndef MCMD_MD_PARTICLE_H
#define MCMD_MD_PARTICLE_H
/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/space/Vector.h>

namespace McMd 
{

    using namespace Util;

    struct MdParticle 
    {
       Vector position;
       Vector velocity;

       MdParticle();
    };

}
#endif