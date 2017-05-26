#ifndef MCMD_INTRA_PAIR_AUTO_CORR_H
#define MCMD_INTRA_PAIR_AUTO_CORR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/SystemAnalyzer.h>   // base class template
#include <mcMd/simulation/System.h>              // base class template parameter
#include <util/accumulators/AutoCorrArray.h>     // base class template parameter
#include <util/space/Vector.h>                   // template parameter
#include <util/containers/DArray.h>              // member template

namespace Simp {
   class Species;
}

namespace McMd
{

   using namespace Util;
   using namespace Simp;

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
   * \sa \ref mcMd_analyzer_IntraPairAutoCorr_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class IntraPairAutoCorr 
    : public SystemAnalyzer<System>
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
      virtual void readParameters(std::istream& in);
  
      /**
      * Load state from an archive.
      *
      * \param ar loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to archive.
      *
      * \param ar saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive. 
      *
      * \param ar      saving or loading archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      /** 
      * Set number of molecules and clear accumulator.
      */
      virtual void setup();
   
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

   private:
   
      /// Output file stream
      std::ofstream outputFile_;

      /// Statistical accumulator
      AutoCorrArray<Vector, double> accumulator_;

      /// Array of separation vectors, one per molecule (temporary)
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
      Analyzer::serialize(ar, version);
      ar & speciesId_;
      ar & atom1Id_;
      ar & atom2Id_;
      ar & capacity_;

      ar & nMolecule_;
      ar & accumulator_;
   }

}
#endif
