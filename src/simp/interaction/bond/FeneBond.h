#ifndef SIMP_FENE_BOND_H
#define SIMP_FENE_BOND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/global.h>
#include <util/param/ParamComposite.h>  // base class
#include <util/random/Random.h>         // Util namespace 

#include <cmath>
#include <sstream>

namespace Simp
{

   using namespace Util;

   /**
   * A finitely-extensible nonlinear element (FENE) bond potential.
   *
   * Finitely extensible nonlinear element (FENE) bond potential 
   * V(r) = -0.5*kappa*r0^2 log( 1 - (r^2/r0^2)) with a spring 
   * constant kappa and maximum bond length r0.
   *
   * This implementation uses the FENE potential for distances less
   * than a cutoff rl at which the FENE potentail reaches a parameter 
   * energyCutoff, and a linear potential for all r > rl.
   *
   * \sa \ref simp_interaction_bond_FeneBond_page "Parameter file format"
   * \sa \ref simp_interaction_bond_interface_page
   *
   * \ingroup Simp_Interaction_Bond_Module
   */
   class FeneBond : public ParamComposite 
   {
   
   public:

      /**
      * Default constructor.
      */
      FeneBond();
   
      /**
      * Copy constructor.
      */
      FeneBond(const FeneBond& other);
   
      /**
      * Assignment.
      */
      FeneBond& operator = (const FeneBond& other);
  
      // Default destructor.
 
      /**
      * Set the number of bond types.
      *
      * \param nBondType number of bond types
      */
      void setNBondType(int nBondType);
       
      /**
      * Read bond interaction parameters from input stream.
      *
      * Format:
      *     kappa        CArray<double>
      *     length       CArray<double>
      *     energyCutoff double
      *
      * \pre nBondType must be set, by calling setNBondType().
      *
      * \param in input stream
      */
      void readParameters(std::istream &in);

      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive &ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive &ar);

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param type   bond type index 
      * \param value  new value of parameter
      */
      void set(std::string name, int type, double value);

      /**
      * Returns potential energy for one bond.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      */
      double energy(double rSq, int type) const;
   
      /**
      * Returns force/distance for one bond, for use in MD.
      *
      * A positive return value represents a repulsive radial force.
      *
      * \param rSq  square of distance between bonded particles.
      * \param type type of bond.
      * \return     scalar repulsive force divided by distance.
      */
      double forceOverR(double rSq, int type) const;

      /**
      * Return bond length chosen from equilibrium distribution.
      *
      * This function returns a bond length chosen from the Boltzmann
      * distribution of lengths for bonds of random orientation. The
      * distribution P(l) of values of the length l is proportional 
      * to l*l*exp[-beta*phi(l) ], where phi(l) is the bond energy. 
      *
      * \param random pointer to random number generator object.
      * \param beta   inverse temperature
      * \param type   bond type
      * \return random bond length chosen from equilibrium distribution.
      */
      double randomBondLength(Random *random, double beta, int type) const;
   
      /**
      * Get a parameter value, identified by a string.
      *
      * \param name parameter name
      * \param type bond type index 
      */
      double get(std::string name, int type) const;

   private:
   
      /// Maximum possible number of bond types
      static const int MaxNBondType = 2;
  
      /// Spring constant.
      double kappa_[MaxNBondType];

      /// Maximum FENE bond length.
      double r0_[MaxNBondType];  

      /// Maximum FENE bond length squared.
      double r0Sq_[MaxNBondType];  

      /// Inverse square maximum length.
      double r0SqInv_[MaxNBondType];  

      /// Prefactor in energy.
      double ce_[MaxNBondType];  

      /// Square of value of r where force is forceCutoff_
      double rlSq_[MaxNBondType];

      /// Value of r where Fene force is forceCutoff_
      double rl_[MaxNBondType];
      
      /// Value of Fene potential at rl.
      double  energyCutoff_[MaxNBondType];  
      
      /// Internal constant
      // double a_[MaxNBondType];
      
      /// Value of Fene force at rl.
      double  forceCutoff_;
      
      /// Number of bond types.
      int    nBondType_;

   };
   
   // Inline method definitions
   
   /* 
   * Return interaction energy for one bond
   */
   inline double FeneBond::energy(double rSq, int type) const
   {
     if (rSq < rlSq_[type]){     
       double g = 1.0 - rSq*r0SqInv_[type];
       return ce_[type]*log(g);
     } else {
       //double b = rl_[type]*kappa_[type]/a_[type];
       return energyCutoff_[type] + forceCutoff_*(sqrt(rSq) - rl_[type]);  
     } 
   }

   /* 
   * Return force / distance for one bond, for use in MD
   */
   inline 
   double FeneBond::forceOverR(double rSq, int type) const
   {
     if (rSq < rlSq_[type]){ 
       double g = 1.0 - rSq*r0SqInv_[type];
       return -kappa_[type]/g;
     } else {
       if (rSq > r0Sq_[type]) {
          std::cout << "Warning: Long FENE bond: Rsq = " << rSq << std::endl;
       }
       //return -rl_[type]*kappa_[type]/(a_[type]*sqrt(rSq));  
       return -forceCutoff_/sqrt(rSq);  
     }              
   }

} 
#endif
