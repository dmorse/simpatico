#ifndef MCMD_CONFIG_IO_H
#define MCMD_CONFIG_IO_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/boundary/Boundary.h>
#include <util/global.h>

#include <iostream>

namespace McMd
{

   using namespace Util;

   class Simulation;
   class System;
   
   /**
   * System configuration file reader and writer.
   *
   * A ConfigIo object is used by a System to read and write the 
   * configuration file for the System. The configuration file contains
   * the Boundary dimensions, the number of molecules of each Species, 
   * and the atomic positions for every Atom in that System.
   *
   * The ConfigIo class defines default implementations of the virtual 
   * read() and write() methods. These define the default config file 
   * format for Simpatico. Other formats may be read and written by using
   * an instance of a subclass of ConfigIo, for which the overridden 
   * read() and write() methods must define the desired file format. 
   *
   * A System holds a pointer to a ConfigIo, which may be set by the
   * System::setConfigIo() method. By default, if no ConfigIo is set
   * explicitly, the System will create and use an instance of 
   * ConfigIo when needed.
   *
   * \ingroup McMd_ConfigIo_Module
   */
   class ConfigIo
   {
   
   public:

      /// Constructor. 
      ConfigIo(System& system);
 
      /// Destructor.   
      virtual ~ConfigIo();
 
      /**
      * Read configuration (particle positions) from file.
      *
      * \param in input file stream.
      */
      virtual void read(std::istream &in) = 0;

      /**
      * Read configuration (particle positions) from file
      * and transform positions from cartesian to 
      * generalized system.
      *
      * \param in input file stream.
      */
      virtual void transformCartConfigToGen(std::istream &in) = 0;
 
      /**
      * Write configuration (particle positions) to file.
      *
      * \param out output file stream.
      */
      virtual void write(std::ostream& out) = 0;

      /**
      * Transform positions from generalized to cartesian
      * system and write configuration (particle positions) 
      * to file.
      *
      * \param out output file stream.
      */
      virtual void transformGenToCartConfig(std::ostream& out) = 0;

   protected:

      /// Get a reference to the parent System. 
      System &system() const;

      /// Get a reference to the parent Simulation. 
      Simulation &simulation() const;

      /// Get the Boundary.
      Boundary &boundary() const;

   private:
   
      /// Boundary object.
      Boundary   *boundaryPtr_;
   
      /// Pointer to parent System;
      System     *systemPtr_;
   
      /// Pointer to parent Simulation.
      Simulation *simulationPtr_;
   
   }; // end class ConfigIo


   // Inline functions 

   /* 
   * Get the parent System. 
   */
   inline System& ConfigIo::system() const
   {
      assert(systemPtr_); 
      return *systemPtr_; 
   }
 
   /* 
   * Get the parent Simulation.
   */
   inline Simulation& ConfigIo::simulation() const
   { 
      assert(simulationPtr_);
      return *simulationPtr_; 
   }

   /* 
   * Get the Boundary.
   */
   inline Boundary& ConfigIo::boundary() const
   {
      assert(boundaryPtr_); 
      return *boundaryPtr_; 
   }


} 
#endif
