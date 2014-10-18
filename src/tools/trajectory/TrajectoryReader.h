#ifndef TOOLS_TRAJECTORY_IO_H
#define TOOLS_TRAJECTORY_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // base class

namespace Tools
{

   class Configuration;
   using namespace Util;

   /**
   * Abstract trajectory file reader.
   *
   * Each concrete subclass of TrajectoryReader implements a specific 
   * trajectory file format.
   *
   * \ingroup Tools_Trajectory_Module
   */
   class TrajectoryReader  : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *
      * \param configuration  parent Configuration object
      * \param isBinary  Is the trajectory file a binary format?
      */
      TrajectoryReader(Configuration& configuration, bool isBinary = false);

      /**
      * Destructor.
      */
      virtual ~TrajectoryReader();

      /**
      * Is the file format binary (true) or text (false)?
      */
      bool isBinary() const;

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

      /**
      * Return parent Configuraiton by reference.
      */
      Configuration& configuration()
      {  return *configurationPtr_; }

   private:

      /// Pointer to parent Tools::Configuration.
      Configuration* configurationPtr_;

      /// Is the file format binary (true) or text (false)?
      bool isBinary_;

   };

   /**
   * Is the file format binary (true) or text (false)?
   */
   inline bool TrajectoryReader::isBinary() const
   {  return isBinary_; }

}
#endif
