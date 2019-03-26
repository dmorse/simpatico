#ifndef MCMD_DIHEDRAL_DIST_H
#define MCMD_DIHEDRAL_DIST_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/analyzers/base/DistributionAnalyzer.h> // base class template
#include <mcMd/simulation/System.h>                   // base class param
#include <util/accumulators/Distribution.h>           // member
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   /**
   * DihedralDist evaluates the distribution of dihedral angles.
   *
   * \sa \ref mcMd_analyzer_DihedralDist_page "parameter file format"
   *
   * \ingroup McMd_Analyzer_McMd_Module
   */
   class DihedralDist : public DistributionAnalyzer<System>
   {
   
   public:

      /**	
      * Constructor.
      *
      * \param system reference to parent System object
      */
      DihedralDist(System &system);

      /**
      * Read parameters from file.  
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
      * Save state to an archive.
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
      * Add particle pairs to DihedralDist histogram.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

   private:

      /// Index of relevant Species (-1 for all)
      int speciesId_;
  
      /// Index of dihedral type (-1 for all)
      int typeId_;
 
      /**
      * Sample dihedrals for a single molecular species.
      *
      * \param is species index
      */ 
      void sampleSpecies(int is);

   };

   // Member function template.

   /*
   * Serialize to/from an archive. 
   */
   template <class Archive>
   void DihedralDist::serialize(Archive& ar, const unsigned int version)
   { 
      DistributionAnalyzer<System>::serialize(ar, version);
      ar & speciesId_;
      ar & typeId_;
   }

}
#endif
