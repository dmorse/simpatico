#ifndef MCMD_INTRA_PAIR_AUTO_CORR_H
#define MCMD_INTRA_PAIR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010, David Morse (morse@cems.umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/diagnostics/SystemDiagnostic.h>   // base class template
#include <mcMd/simulation/System.h>              // base class template parameter
#include <util/accumulators/AutoCorrArray.h>     // base class template parameter
#include <util/space/Vector.h>                   // template parameter
#include <util/containers/DArray.h>              // member template

namespace McMd
{

   using namespace Util;

   class Species;

   /**
   * Autocorrelation for vector separation of any two monomers on a molecule.
   * 
   * This class evaluates the autocorrelation function for the vector that
   * separates two atoms with indices atom1Id and atom2Id on a molecule of
   * species speciesId. It can be used to evaluate the autocorrelation
   * function for individual bonds or for the end-to-end vector of a linear
   * polymer.
   *
   * If atom1Id = i and atom2Id = j, let \f$ V(t) = R_i(t) - R_j(t) \f$. 
   * This class calculates and outputs the function  
   * \f[
   *     F(\tau) \equiv \langle V(t) \cdot V(t-\tau) \rangle
   * \f] 
   * where \f$ \tau \f$ is a time separation, and where 
   * \f$ \langle \cdots \rangle\f$ denotes an * average over values 
   * of the time t and over all molecules in a species.  Here, 
   * \f$ V(t) \cdot V(t') \f$ denotes a dot product of values of the 
   * vector \f$ V(t) \f$  at different times.
   *
   * \ingroup McMd_Diagnostic_Module
   */
   class IntraPairAutoCorr 
    : public SystemDiagnostic<System>
   {
   
   public:
  
      /**
      * Constructor.
      *
      * \param system reference to parent System
      */
      IntraPairAutoCorr(System &system);
  
      /** 
      * Read parameters from file.
      *
      * Input format:
      *
      *   - int    interval        sampling interval
      *   - string outputFileName  output file name
      *   - int    speciesId       integer Id for molecular Species
      *   - int    atom1Id         local atom Id for 1st atom
      *   - int    atom2Id         local atom Id for 2nd atom
      *   - int    capacity        capacity of array of previous values.
      *
      * \param in input parameter stream
      */
      virtual void readParam(std::istream& in);
  
      /** 
      * Set number of molecules and clear accumulator.
      */
      virtual void initialize();
   
      /** 
      * Evaluate separation vector for all chains, add to ensemble.
      *
      * \param iStep counter for number of steps
      */
      virtual void sample(long iStep);
   
      /**
      * Output results after simulation is completed.
      */
      virtual void output();

      /**
      * Save state to binary file archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchiveType& ar);

      /**
      * Load state from a binary file archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void load(Serializable::IArchiveType& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

   private:
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Statistical accumulator
      AutoCorrArray<Vector, double> accumulator_;

      /// Array of separation vectors, one per molecule of species.
      DArray<Vector> data_;
   
      /// Pointer to relevant Species.
      Species *speciesPtr_;
   
      /// Index of relevant Species.
      int     speciesId_;
   
      /// Number of molecules in the species (must not change).
      int     nMolecule_;
   
      /// Local index of atom1
      int     atom1Id_;
   
      /// Local index of atom2 (must be greater than atom1Id_).
      int     atom2Id_;
   
      /// Maximum length of each sequence in AutoCorrArray.
      int     capacity_;
   
      /// Has readParam been called?
      bool    isInitialized_;

   };

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void IntraPairAutoCorr::serialize(Archive& ar, const unsigned int version)
   {
      if (!isInitialized_) {
         UTIL_THROW("Error: Object not initialized.");
      }

      // Data not set by readParam
      ar & accumulator_;
      ar & data_;
      ar & nMolecule_;

      // Data set by readParam (check consistency).
      serializeCheck(ar, speciesId_, "speciesId");
      serializeCheck(ar, atom1Id_, "atom1Id");
      serializeCheck(ar, atom2Id_, "atom2Id");
      serializeCheck(ar, capacity_, "capacity");

   }

}
#endif
