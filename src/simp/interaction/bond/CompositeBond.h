#ifndef SIMP_COMPOSITE_BOND_H
#define SIMP_COMPOSITE_BOND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/param/Begin.h>
#include <util/param/End.h>
#include <util/global.h>

#include <math.h>

namespace Util 
{  class Random; }

namespace Simp
{

   using namespace Util;

   /**
   * Composite bond potential.
   *
   * This class template creates a bond potential energy which is the 
   * sum of a BareBond bond potential and a BarePair pair potential.
   * This is intended for use in models in which the nonbonded pair 
   * interactions between bonded atoms are masked, to allow the pair
   * potential to be explicitly added back in.
   * 
   * \ingroup Simp_Interaction_Bond_Module
   */
   template <class BareBond, class BarePair>
   class CompositeBond : public ParamComposite 
   {
   
   public:
   
      /**
      * Constructor.
      */
      CompositeBond();

      // Mutators

      /**  
      * Set nAtomType value.
      *
      * \param nBondType number of atom types.
      */
      void setNBondType(int nBondType);

      /**
      * Read epsilon and sigma, initialize other variables.
      *
      * \pre nAtomType must be set, by calling setNAtomType().
      *
      * \param in  input stream 
      */
      void readParameters(std::istream &in);
      
      /**
      * Modify a parameter, identified by a string.
      *
      * \param name  parameter name
      * \param type  bond type index 
      * \param value new value of parameter
      */
      void set(std::string name, int type, double value);

      // Accessors
      
      /**
      * Returns interaction energy for a single pair of particles. 
      *
      * \param rsq  square of distance between pair of atoms
      * \param type type index for bond.
      * \return     bond interaction energy
      */
      double energy(double rsq, int type) const;
   
      /**
      * Returns ratio of bond force to atom separation.
      *
      * Multiply this quantity by the separation vector to obtain
      * the force vector. A positive value for the return value
      * represents a repulsive bond force.
      *
      * \param rsq        square of distance between particles
      * \param bondTypeId type index for bond.
      * \return force divided by distance 
      */
      double forceOverR(double rsq, int bondTypeId) const;

      /**
      * Throws exception if called.
      * 
      * \param randomPtr pointer to a random number generator
      * \param beta      inverse absolute temperature (inverse energy)
      * \param type      type index for bond.
      */  
      double randomBondLength(Random* randomPtr, double beta, int type) 
             const;

      /**
      * Get a parameter value, identified by a string.
      *
      * \param  name parameter name
      * \param  type bond type index
      * \return value of parameter
      */
      double get(std::string name, int type) const;

   private:

      // Bare bond potential object.
      BareBond bond_;

      // Bare pair potential object.
      BarePair pair_;
      
      /**
      * Copy constructor.
      */
      CompositeBond(const CompositeBond<BareBond, BarePair>& other);
 
   };
     
 
   /* 
   * Constructor.
   */
   template <class BareBond, class BarePair>
   CompositeBond<BareBond, BarePair>::CompositeBond() 
    : bond_(),
      pair_()
   {
      std::string name("CompositeBond");
      name += "<";
      name += bond_.className();
      name += ",";
      name += pair_.className();
      name += ">";
      setClassName(name.c_str());
   }
   
   
   /* 
   * Copy constructor.
   */   
   template <class BareBond, class BarePair>
   CompositeBond<BareBond, BarePair>::CompositeBond(const CompositeBond<BareBond, BarePair>& other) 
    : bond_(other.bond_),
      pair_(other.pair_)
   {}
     
   /* 
   * Read potential parameters from file.
   */  
   template <class BareBond, class BarePair>
   void CompositeBond<BareBond, BarePair>::readParameters(std::istream &in)
   {

      readParamComposite(in, bond_);
      readParamComposite(in, pair_);

      #if 0
      std::string label;
      Begin* beginPtr = 0;
      End*   endPtr = 0;
      // Read bare bond potential
      beginPtr = &addBegin(bond_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParameters(in);
      readParamComposite(in, bond_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParameters(in);

      // Read bare pair potential
      beginPtr = &addBegin(pair_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParameters(in);
      readParamComposite(in, pair_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParameters(in);
      #endif

   }
   
   /* 
   * Set nAtomType
   */
   template <class BareBond, class BarePair>
   void CompositeBond<BareBond, BarePair>::setNBondType(int nBondType) 
   {  
      bond_.setNBondType(nBondType);
      pair_.setNAtomType(1);
   }   

   /* 
   * Calculate interaction energy for a bond.
   */
   template <class BareBond, class BarePair>
   inline double 
   CompositeBond<BareBond, BarePair>::energy(double rsq, int type)
   const 
   {
      double total = bond_.energy(rsq, type);
      total += pair_.energy(rsq, 0, 0);
      return total;
   }
  
   /* 
   * Compute force/distance for a bond.
   */
   template <class BareBond, class BarePair>
   inline double 
   CompositeBond<BareBond, BarePair>::forceOverR(double rsq, int type)
   const 
   {
      double total = bond_.forceOverR(rsq, type);
      total += pair_.forceOverR(rsq, 0, 0);
      return total;
   }

   /*
   * Modify a parameter, identified by a string.
   */
   template <class BareBond, class BarePair>
   void CompositeBond<BareBond, BarePair>
        ::set(std::string name, int type, double value)
   {
      UTIL_THROW("Unrecognized parameter name");
   }

   /*
   * Get a parameter value, identified by a string.
   */
   template <class BareBond, class BarePair>
   double CompositeBond<BareBond, BarePair>::
          get(std::string name, int type) const
   {
      UTIL_THROW("Unrecognized parameter name");
      return 0.0;
   }

   /*
   * Throws exception if called.
   */  
   template <class BareBond, class BarePair>
   double CompositeBond<BareBond, BarePair>
          ::randomBondLength(Random* randomPtr, double beta, int type) const
   {  
      UTIL_THROW("Unimplemented function"); 
      return 0.0; // To avoid compiler warnings.
   }

}
#endif
