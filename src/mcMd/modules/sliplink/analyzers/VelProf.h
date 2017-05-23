#ifndef VEL_PROF_H
#define VEL_PROF_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2014, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>  // base class template
#include <mcMd/simulation/System.h>             // class template parameter
#include <util/accumulators/Average.h>          // member
#include <util/containers/DArray.h>             // member template
#include <util/space/Vector.h>                   // member template parameter

#include <cstdio> 
#include <cstring> 

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

   /**
   *
   * Evaluates x-velocity profile as a function of z
   *
   * \ingroup Analyzer_Module
   */
   class VelProf : public SystemAnalyzer<System>
   {
   
   public:

      /**
      * Constructor.
      *
      * \param system reference to parent System object
      */
      VelProf(System &system);
   
      /**
      * Read parameters from file, and allocate data arrays.
      *
      * \param in parameter input stream
      */
      virtual void readParameters(std::istream& in);
   
      /** 
      * Empty.
      */
      virtual void setup();
   
      /**
      * Evaluate and output x-velocity profile as a function of z.
      *
      * \param iStep step counter
      */
      virtual void sample(long iStep);
   
      /**
      * Output param file at end of simulation.
      */
      virtual void output();

   private:

      /// Output file stream
      std::ofstream outputFile_;

      /// number of slabs in which Lz is divided
      int nSlabs_;
      
      /// Number of samples per block average output.
      int nSamplePerBlock_;
      
      /// Array of Average objects - statistical accumulators
      DArray<Average>  accumulator_;
      
      /// Total transferred momentum.
      double   totalMomentum_;
      
      /// velx_[i] contains the mean x-velocity for slab i
      DArray<double> velx_;
      
      /// nAtomi_[i] contains the number of atoms in slab i
      DArray<int> nAtomi_;    
         
      int nSpec_;
      
      double slabWidth_;
      
      int    iBlock_;       
      
   };

}
#endif
