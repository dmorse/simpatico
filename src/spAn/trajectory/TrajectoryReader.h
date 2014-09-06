#ifndef SPAN_TRAJECTORY_IO_H
#define SPAN_TRAJECTORY_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace SpAn
{

   class Configuration;
   using namespace Util;

   /**
   * Abstract trajectory file reader.
   *
   * Each concrete subclass of TrajectoryReader implements a specific 
   * trajectory file format.
   *
   * \ingroup SpAn_Trajectory_Module
   */
   class TrajectoryReader  : public ParamComposite
   {

   public:

      /**
      * Default constructor.
      */
      TrajectoryReader();

      /**
      * Constructor.
      *
      * \param configuration parent Configuration object
      */
      TrajectoryReader(Configuration& configuration);

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader();

      /**
      * Read a header (if any).
      *
      * \param file input file 
      */
      virtual void readHeader(std::ifstream& file){}

      /**
      * Read a frame.
      *
      * \param file input file 
      * \return true if a frame was found, false if eof
      */
      virtual bool readFrame(std::ifstream& file) = 0;

   protected:

      Configuration& configuration()
      {  return *configurationPtr_; }

   private:

      Configuration* configurationPtr_;

   };

}
#endif
