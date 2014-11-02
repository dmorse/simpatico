#ifndef UTIL_AUTOCORRELATION_H
#define UTIL_AUTOCORRELATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "AutoCorrelation.h"          // base class
#include <util/containers/GArray.h>   // member
#include <util/global.h>

namespace Util
{

   /**
   * Hierarchical auto-correlation function algorithm.
   *
   * This class represents the primary stage of a linked list of
   * AutoCorrelation objects.
   *
   * \ingroup Accumulators_Module
   */
   template <typename Data, typename Product>
   class AutoCorrelation : public AutoCorrStage
   {

   public:

      /**
      * Constructor
      *
      * \param blockFactor ratio of block sizes of subsequent stages
      */
      AutoCorrelation(int blockFactor = 4);

      /**
      * Destructor.
      *
      * Recursively destroy all descendant stages.
      */
      virtual ~AutoCorrelation();

      /**
      * Set the base output file name.
      */
      void setFileName(int bufferCapacity);

      /**
      * Set buffer capacity, allocate memory and initialize.
      *
      * \param bufferCapacity max. number of values in buffer
      */
      void setCapacity(int bufferCapacity);

      /**
      * Register the creation of a descendant stage.
      *
      * This should be called only by a root stage.
      *
      * \param ptr pointer to a descendant AutoCorrelation.
      */
      virtual void registerDescendant(AutoCorrelation* ptr);

      ///\name Serialization
      //@{

      /**
      * Load state from an archive.
      *
      * \param ar binary loading (input) archive.
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save state to an archive.
      *
      * \param ar binary saving (output) archive.
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Serialize to/from an archive.
      *
      * \param ar      archive
      * \param version archive version id
      */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version);

      //@}
      ///\name Accessors
      //@{

      /**
      * Output the autocorrelation function
      *
      * \param out output stream.
      */
      virtual void output(std::ostream& out);

      //@}

   private:

      // Pointers to descendant AutoCorrStage objects
      GArray<AutoCorrStage*> descendants_;

   };

}
#endif
