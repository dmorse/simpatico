#ifndef DDMD_POTENTIAL_H
#define DDMD_POTENTIAL_H

#include <util/param/ParamComposite.h>  // base class
#include <util/util/Setable.h>          // template for members
#include <util/space/Tensor.h>          // parameter for member

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

namespace DdMd
{

   using namespace Util;

   /**
   * A Potential represents a potential energy contribution.
   *
   * A Potential provides methods to calculate potential energy,
   * and associated atomic forces and stress.
   *  
   *  \ingroup DdMd_Potential_Module
   */
   class Potential : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Potential();

      /**
      * Destructor.
      */
      virtual ~Potential();

      /**
      * Set flag to identify if reverse communication is enabled.
      *
      * \param reverseUpdateFlag true if reverse communication is enabled.
      */
      void setReverseUpdateFlag(bool reverseUpdateFlag);

      /**
      * Get flag to identify if reverse communication is enabled.
      */
      bool reverseUpdateFlag() const;

      /// \name Total Energy, Force and Stress 
      //@{

      /**
      * Add force contributions to all atomic forces.
      */
      virtual void addForces() = 0;

      /**
      * Compute potential energy on all processors.
      *
      * This method must be called on all processors. The result is
      * stored on the master processor, and may be retrieved by 
      * calling energy() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computeEnergy(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeEnergy() = 0;
      #endif

      /**
      * Return the total potential, from all processors.
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeEnergy.
      */
      double energy() const;

      /**
      * Compute stress on all processors.
      *
      * This method must be called on all processors. The result 
      * is stored on the master processor, and may be retrieved 
      * by calling energy() on this processor.
      */
      #ifdef UTIL_MPI
      virtual void computeStress(MPI::Intracomm& communicator) = 0;
      #else
      virtual void computeStress() = 0;
      #endif

      /**
      * Return the stress tensor.
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeStress.
      */
      Util::Tensor stress() const;

      /**
      * Return the pressure. 
      *
      * This method should only be called on the master (rank 0) 
      * processor, after a previous call to computeStress.
      */
      double pressure() const;

      //@}

   protected:

      /**
      * Set a value for the total energy.
      */
      void setEnergy(double energy);

      /**
      * Mark the energy as unknown (nullify).
      */
      void unsetEnergy();

      /**
      * Set a value for the total stress.
      */
      void setStress(const Tensor& stress);

      /**
      * Mark the stress as unknown (nullify).
      */
      void unsetStress();

      /*
      * Add a pair contribution to the stress tensor.
      */
      void incrementPairStress(const Vector& f, const Vector& dr, 
                               Tensor& stress) const;
   private:

      /// Total stress.
      Setable<Tensor> stress_;

      /// Total energy.
      Setable<double> energy_;

      /// Is reverse update communication enabled?
      bool reverseUpdateFlag_;

   };

   inline bool Potential::reverseUpdateFlag() const
   {  return reverseUpdateFlag_; }

   /*
   * Add a pair contribution to the virial tensor (protected).
   */
   inline void 
   Potential::incrementPairStress(const Vector& f, const Vector& dr, 
                                  Tensor& stress) const
   {
      int i, j;
      for (i = 0; i < Dimension; ++i) {
         for (j = 0; j < Dimension; ++j) {
            stress(i, j) += f[i]*dr[j];
         }
      }
   }

}
#endif
