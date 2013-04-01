#ifndef MCMD_EXTERNAL_POTENTIAL_H
#define MCMD_EXTERNAL_POTENTIAL_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>   // base class

#include <string>

namespace Util
{
   class Vector;
   class Tensor;
   class Random;
}

namespace McMd
{

   using namespace Util;

   class Atom;


   /**
   * Abstract External Potential class.
   *
   * \ingroup McMd_External_Module
   */
   class ExternalPotential : public ParamComposite
   {

   public:

      /**
      * Constructor .
      */
      ExternalPotential()
      {  setClassName("ExternalPotential"); }

      /**
      * Destructor (does nothing)
      */
      virtual ~ExternalPotential()
      {}

      //@{ External Interaction Interface

      /**
      * Returns external potential energy of a single particle. 
      *
      * \param position atomic position Vector
      * \param i        atom type.
      * \return external potential energy
      */
      virtual double energy(const Vector& position, int i) const = 0;

      /**
      * Returns force caused by the external potential.
      *
      * \param position  atom position
      * \param type      atom type id
      * \param force     force on the atom (on output)
      */
      virtual void getForce(const Vector& position, int type, Vector& force) const = 0;

      //@}

      /**
      * Add external force of an Atom to the total force acting on it.
      */
      virtual void addForces() = 0;

      /**
      * Return total external energy of this System.
      */
      virtual double energy() const = 0;

      /**
      * Calculate the external energy for one Atom.
      */
      virtual double atomEnergy(const Atom& atom) const = 0;

      /**
      * Return name of external interaction class (e.g., "LamellarOrderingExternal").
      */
      virtual std::string interactionClassName() const = 0;
   
   };

}
#endif
