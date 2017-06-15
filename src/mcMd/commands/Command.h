#ifndef MCMD_COMMAND_H
#define MCMD_COMMAND_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>
#include <util/random/Random.h>
#include <util/global.h>

namespace McMd
{

   using namespace Util;

   class Simulation;

   /**
   * Command is an object that can be invoked from the command script.
   *
   * \ingroup McMd_Command_Module
   */
   class Command : public ParamComposite
   {

   public:

      /**
      * Constructor.
      */
      Command();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~Command();

      /**
      * Read required parameters from file.
      *
      * Empty default implementation.
      */
      virtual void readParameters(std::istream &in);

      /**
      * Load internal state from an archive.
      *
      * Empty default implementation.
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
      * Read and execute command.
      *
      * \return true if name is recognized, false if not
      */
      virtual bool readCommand(std::string const & name, std::istream& in);

      // Accessor Functions

      /**
      * Output statistics for this move (called at the end of the simulation)
      */
      virtual void output();

   };

}
#endif
