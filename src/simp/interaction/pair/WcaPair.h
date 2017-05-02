#ifndef SIMP_WCA_PAIR_H
#define SIMP_WCA_PAIR_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <simp/interaction/pair/LJPair.h>

#include <iostream>

namespace Simp
{

   using namespace Util;

   /**
   * A Weeks-Chandler-Anderson (WCA) repulsive LJ interaction.
   *
   * The WCA potential is given by a Lennard-Jones potential cutoff at
   * the potential minimum at sigma 2^{1/6}, and shifted so that the
   * potential is zero at the cutoff distance.
   *
   * \sa \ref simp_interaction_pair_WcaPair_page "Parameter file format"
   * \sa \ref simp_interaction_pair_interface_page
   * \sa \ref simp_interaction_pair_page
   * 
   * \ingroup Simp_Interaction_Pair_Module
   */
   class WcaPair : public LJPair
   {
   
   public:
   
      /**
      * Read epsilon and sigma, initialize other variables.
      *
      * \pre nAtomType must have been set, by calling setNAtomType().
      *
      * \param in  input stream 
      */
      void readParameters(std::istream &in);

      /**
      * Modify a parameter, identified by a string.
      *
      * \param name   parameter name
      * \param i      atom type index 1
      * \param j      atom type index 2
      * \param value  new value of parameter
      */
      void set(std::string name, int i, int j, double value);

   };

}
#endif  
