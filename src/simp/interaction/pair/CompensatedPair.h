#ifndef SIMP_COMPENSATED_PAIR_H
#define SIMP_COMPENSATED_PAIR_H

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

namespace Simp
{

   using namespace Util;

   /**
   * Compensated pair potential for transient gel.
   *
   * This class is intended for use in a model with transient crosslinks
   * that can be created or destroyed at fixed chemical potential.
   * 
   * \ingroup Simp_Interaction_Potential_Module
   */
   template <class BarePair, class LinkPotential>
   class CompensatedPair : public ParamComposite 
   {
   
   public:
   
      /**
      * Constructor.
      */
      CompensatedPair();

      /// \name Mutators

      /**
      * Read epsilon and sigma, initialize other variables.
      *
      * \pre nAtomType must be set, by calling setNAtomType().
      *
      * \param in  input stream 
      */
      void readParameters(std::istream &in);
      
      /**  
      * Set nAtomType value.
      *
      * \param nAtomType number of atom types.
      */
      void setNAtomType(int nAtomType);
      
      /**
      * Set LJ interaction energy for a specific pair of Atom types.
      *
      * \param i        type of Atom 1
      * \param j        type of Atom 2
      * \param epsilon  LJ energy parameter
      */
      void setEpsilon(int i, int j, double epsilon);     

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value);

      /// \name Accessors
      
      /**
      * Returns interaction energy for a single pair of particles. 
      *
      * \param rsq square of distance between particles
      * \param i   type of particle 1
      * \param j   type of particle 2
      * \return    pair interaction energy
      */
      double energy(double rsq, int i, int j) const;
   
      /**
      * Returns ratio of scalar pair interaction force to pair separation.
      *
      * Multiply this quantity by the components of the separation vector
      * to obtain the force vector. A positive value for the return value
      * represents a repulsive force between a pair of particles.
      *
      * \param rsq square of distance between particles
      * \param i type of particle 1
      * \param j type of particle 2
      * \return  force divided by distance 
      */
      double forceOverR(double rsq, int i, int j) const;
   
      /**
      * Get cutoff for pair potential, for a particular type pair.
      */
      double cutoffSq(int i, int j) const;
 
      /**
      * Get maximum of pair cutoff distance, for all atom type pairs.
      */
      double maxPairCutoff() const;      
 
      /// \name Other Accessors
      
      /**
      * Get LJ interaction energy for a specific pair of Atom types.
      *
      * \param i   type of Atom 1
      * \param j   type of Atom 2
      * \return    epsilon_[i][j]
      */
      double epsilon(int i, int j) const;      

      /**
      * Get a parameter value, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      */
      double get(std::string name, int i, int j) const;
 
   private:
      
      /// Maximum allowed link length.
      double maxLinkLength_;
           
      /// Square of maximum link length.
      double maxLinkLengthSq_;

      /**
      * Maximum pair potential cutoff radius, for all monomer type pairs.
      *
      * Used in construction of a cell list or Verlet pair list.
      */
      double maxPairCutoff_;

      /**
      * Square of maximum pair potential cutoff radius.
      */
      double maxPairCutoffSq_;
           
      BarePair        pair_;
      LinkPotential   link_;
      double          temp_;
      double          mu_;
      double          activity_;
      double          beta_;

      /**
      * Copy constructor.
      */
      CompensatedPair(BarePair barePot);
 
      /**
      * Copy constructor.
      */
      CompensatedPair(const CompensatedPair<BarePair, LinkPotential>& other);

   };
  
      
   /* 
   * Constructor.
   */
   template <class BarePair, class LinkPotential>
   CompensatedPair<BarePair, LinkPotential>::CompensatedPair() 
    : maxPairCutoff_(0.0),
      maxPairCutoffSq_(0.0),
      temp_(0.0),
      mu_(0.0),
      activity_(0.0),
      beta_(0.0)
   {
      std::string name("CompensatedPair");
      name += "<";
      name += pair_.className();
      name += ",";
      name += link_.className();
      name += ">";
      setClassName(name.c_str());
   }
   
   /* 
   * Copy constructor.
   */   
   template <class BarePair, class LinkPotential>
   CompensatedPair<BarePair, LinkPotential>::CompensatedPair(BarePair barePot) 
    : maxPairCutoff_(0.0),
      pair_(barePot),
      temp_(0.0),
      mu_(0.0),      
      activity_(0.0),
      beta_(0.0)
   {}
     
   /* 
   * Read potential parameters from file.
   */  
   template <class BarePair, class LinkPotential>
   void CompensatedPair<BarePair, LinkPotential>::readParameters(std::istream &in)
   {
      #if 0
      std::string label;
      Begin* beginPtr = 0;
      End*   endPtr = 0;

      // Read bare pair potential
      //label = pair_.className();
      //label += "{";
      beginPtr = &addBegin(pair_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParameters(in);
      readParamComposite(in, pair_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParameters(in);

      // Read bare link potential
      //label = link_.className();
      //label += "{";
      beginPtr = &addBegin(link_.className().c_str());
      beginPtr->setIndent(*this);
      beginPtr->readParameters(in);
      readParamComposite(in, link_);
      endPtr = &addEnd();
      endPtr->setIndent(*this);
      endPtr->readParameters(in);
      #endif

      readParamComposite(in, pair_);
      readParamComposite(in, link_);

      read<double>(in, "maxLinkLength", maxLinkLength_);	
      read<double>(in, "temperature", temp_);
      read<double>(in, "mu", mu_);

      beta_ = 1.0/temp_;
      activity_ = exp(mu_*beta_);

      maxLinkLengthSq_ = maxLinkLength_*maxLinkLength_;
      maxPairCutoff_ = maxLinkLength_;
      if (pair_.maxPairCutoff() > maxLinkLength_) {
         maxPairCutoff_ = pair_.maxPairCutoff();
      }
      maxPairCutoffSq_ = maxPairCutoff_*maxPairCutoff_;
   }
   
   /* 
   * Set nAtomType
   */
   template <class BarePair, class LinkPotential>   
   void CompensatedPair<BarePair, LinkPotential>::setNAtomType(int nAtomType) 
   {  
      pair_.setNAtomType(nAtomType);
      link_.setNBondType(1);
   }   

   /* 
   * SetEpsilon
   */
   template <class BarePair, class LinkPotential>   
   void CompensatedPair<BarePair, LinkPotential>::setEpsilon(int i, int j, double epsilon) 
   {  pair_.setEpsilon(i, j, epsilon); }   


   /* 
   * Get maximum of pair cutoff distance, for all atom type pairs.
   */
   template <class BarePair, class LinkPotential>  
   double CompensatedPair<BarePair, LinkPotential>::maxPairCutoff() const
   { return maxPairCutoff_; }


   /* 
   * Calculate interaction energy for a pair, as function of squared distance.
   */
   template <class BarePair, class LinkPotential>
   inline double CompensatedPair<BarePair, LinkPotential>::energy(double rsq, int i, int j) const 
   {
        double total = pair_.energy(rsq, i, j);
        if (rsq < maxLinkLengthSq_) {
           total += activity_ * temp_ * exp(-beta_ * link_.energy(rsq, 0));
        }
        return total;
   }
 
  
   /* 
   * Calculate force/distance for a pair as function of squared distance.
   */
   template <class BarePair, class LinkPotential>
   inline double 
   CompensatedPair<BarePair, LinkPotential>::forceOverR(double rsq, int i, int j) const
   {
        double total = 0.0;
        if (rsq < pair_.cutoffSq(i, j)) {
           total += pair_.forceOverR(rsq, i, j);
        }
        if (rsq < maxLinkLengthSq_) {
           total -= activity_ * exp(-beta_ * link_.energy(rsq, 0)) * link_.forceOverR(rsq, 0);
        }
        return total;
   }

   /*
   * Get square of maximum pair cutoff distance.
   */
   template <class BarePair, class LinkPotential>
   inline 
   double CompensatedPair<BarePair, LinkPotential>::cutoffSq(int i, int j) const
   {  return maxPairCutoffSq_; }
 
   /*
   * Modify a parameter, identified by a string.
   */
   template <class BarePair, class LinkPotential>
   void CompensatedPair<BarePair, LinkPotential>
        ::set(std::string name, int i, int j, double value)
   {
      if (name == "epsilon") {
         pair_.setEpsilon(i, j, value);
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
   }

   /*
   * Get a parameter value, identified by a string.
   */
   template <class BarePair, class LinkPotential>
   double CompensatedPair<BarePair, LinkPotential>
          ::get(std::string name, int i, int j) const
   {
      double value = 0.0;
      if (name == "epsilon") {
         value = pair_.epsilon(i, j);
      } else {
         UTIL_THROW("Unrecognized parameter name");
      }
      return value;
   }

}
#endif
