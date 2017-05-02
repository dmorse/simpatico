#ifndef MCMD_GROUP_REBRIDGE_BASE_H
#define MCMD_GROUP_REBRIDGE_BASE_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mcMoves/SystemMove.h>  // base class
#include <util/boundary/Boundary.h>
#include <mcMd/mcSimulation/McSystem.h>
#ifdef SIMP_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#include <util/space/Vector.h>

namespace McMd
{

   using namespace Util;

   class Atom;
   class McSystem;

   /**
   * Base class for rebridging a group of atoms forming a tetrahedron.
   *
   * \ingroup McMd_McMove_Module MD_Module
   */
   class GroupRebridgeBase : public SystemMove 
   {
   
   public:
   
      /**
      * Constructor. 
      */
      GroupRebridgeBase(McSystem& system);
   
      /**
      * Destructor. 
      */
      virtual ~GroupRebridgeBase();

      /**
      * Calculate the energy cost for rebriding a tetra group.
      *
      * \param aPtr     ptr to a atom
      * \param bPtr     ptr to b atom
      * \param cPtr     ptr to c atom
      * \param dPtr     ptr to d atom
      * \param bondType bondType
      * \param energy   energy cost of rebridging
      */
      void tetraEnergy(Atom* aPtr, Atom* bPtr, Atom* cPtr, Atom* dPtr,
                       int bondType, double &energy);

      #ifdef SIMP_ANGLE
      /**
      * Calculate the angle energy for a bead triple.
      *
      *      a -- b -- c
      *
      * \param a     atom a
      * \param b     atom b
      * \param c     atom c
      * \param type  angle type
      */
      double angleEnergy(const Atom &a, const Atom &b, const Atom &c, int type);
      #endif

      /**
      * Calculate the energy cost for rebriding an octa group.
      *
      * \param aPtr     ptr to a atom
      * \param bPtr     ptr to b atom
      * \param cPtr     ptr to c atom
      * \param dPtr     ptr to d atom
      * \param mPtr     ptr to d atom
      * \param nPtr     ptr to d atom
      * \param bondType bondType
      * \param energy   energy cost of rebridging
      */
      void octaEnergy(Atom* aPtr, Atom* bPtr, Atom* cPtr, Atom* dPtr,
                      Atom* mPtr, Atom* nPtr, int bondType, double &energy);
  
   };

   #ifdef SIMP_ANGLE
   /*
   */
   inline
   double GroupRebridgeBase::angleEnergy(const Atom &a, const Atom &b,
                       const Atom &c, int type)
   {
      Vector r1, r2;
      double cosTheta;
      system().boundary().distanceSq(b.position(), a.position(), r1);
      system().boundary().distanceSq(c.position(), b.position(), r2);
      r1 /= r1.abs();
      r2 /= r2.abs();
      cosTheta = r1.dot(r2);
      return system().anglePotential().energy(cosTheta, type);
   }
   #endif

}      
#endif
