#ifndef MCMD_NVE_VV_INTEGRATOR_H
#define MCMD_NVE_VV_INTEGRATOR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/mdIntegrators/MdIntegrator.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   /**
   * An NVE Verlet molecular dynamics integrator.
   *
   * \ingroup McMd_MdIntegrator_Module
   */
   class NveVvIntegrator : public MdIntegrator
   {
   
   public:

      /// Constructor. 
      NveVvIntegrator(MdSystem& system);
 
      /// Destructor.   
      virtual ~NveVvIntegrator();

      /**
      * Read parameters from file and initialize this MdSystem.
      *
      * \param in input file stream.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load the internal state to an archive.
      *
      * \param ar archive object.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save the internal state to an archive.
      *
      * \param ar archive object.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Setup private variables before main loop.
      */
      virtual void setup();

      /**
      * Take a complete NVE MD integration step.
      */
      virtual void step();

   private:

      /// Factors of 0.5*dt/mass for different atom types.
      DArray<double> prefactors_;

   }; 

} 
#endif
