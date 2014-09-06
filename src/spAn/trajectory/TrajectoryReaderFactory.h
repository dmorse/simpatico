#ifndef SPAN_TRAJECTORY_READER_FACTORY_H
#define SPAN_TRAJECTORY_READER_FACTORY_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/Factory.h>  
#include <spAn/trajectory/TrajectoryReader.h>

#include <string>

namespace SpAn
{

   using namespace Util;
   class Configuration;

   /**
   * Default Factory for subclasses of TrajectoryReader.
   *
   * \ingroup SpAn_Trajectory_Module
   */
   class TrajectoryReaderFactory : public Factory<TrajectoryReader> 
   {

   public:

      /**
      * Constructor
      *
      * \param configuration parent Configuration object
      */
      TrajectoryReaderFactory(Configuration& configuration);

      /**
      * Create an instance of a specified subclass of TrajectoryReader.
      *
      * If the subclassName is recognized, this method returns a
      * pointer to new object. If the name is not recognized, it
      * returns a null pointer.
      *
      * The new object is created using a constructor that takes
      * the parent SpAn::Configuration as a parameter. The calling
      * function must invoke DdMd::TrajectoryReader::initialize() before
      * using the new TrajectoryReader to read or write configuration
      * files.
      *
      * \param  subclassName  name of a subclass of TrajectoryReader 
      * \return TrajectoryReader*     pointer to new instance of subclassName
      */
      TrajectoryReader* factory(const std::string &subclassName) const;

   private:

      /// Pointer to a parent Configuration.
      Configuration* configurationPtr_;

   };

}
#endif
