#ifdef  MCMD_PERTURB
#ifdef MCMD_EXTERNAL
#ifndef MC_EXTERNAL_PERTURBATION_H
#define MC_EXTERNAL_PERTURBATION_H


/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/perturb/LinearPerturbation.h>      // base class
#include <mcMd/potentials/external/ExternalPotentialImpl.h>   

namespace McMd
{

   using namespace Util;

   class McSystem;

   /**
   * A Perturbation in the external potential external parameter.
   *
   * An McExternalPerturbation describes a Hamiltonian with a variable
   * external parameter.
   *
   * \ingroup Perturb_Module
   */
   class McExternalPerturbation : public LinearPerturbation<McSystem>
   {

   public:
  
      /**
      * Constructor. 
      */
      McExternalPerturbation(McSystem& system);

      /**
      * Destructor
      */
      virtual ~McExternalPerturbation();

      /**
      * Read external parameter from file.
      *
      * \param in input stream (file or std in).
      */ 
      virtual void readParam(std::istream& in);
 
      /**
      * Set external parameter for this System.
      */
      virtual void setParameter();

      /**
      * Return the external parameter for this System.
      *
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double parameter(int i) const;

      /**
      * Return external potential energy / ( kT *external parameter )
      *
      * \param i index of the perturbation parameter. 
      * In this case, i is always 0 since nParameters = 1.
      */
      virtual double derivative(int i) const;
      
   private: 

      /*
      Number of perturbation parameters associated with a System.
      nParameters = 1 for McExternalPerturbation.
      */
      int nParameters_;
   };

}
#endif
#endif  // #ifdef MCMD_EXTERNAL
#endif // #ifdef MCMD_PERTURB 
