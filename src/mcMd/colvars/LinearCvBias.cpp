/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "LinearCvBias.h"

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   LinearCvBias::LinearCvBias() 
   {}

   /*
   * Default destructor.
   */
   LinearCvBias::~LinearCvBias()
   {}

   /*
   * Read parameters from file.
   */
   void LinearCvBias::readParameters(std::istream& in)
   {  read<double>(in, "k", k_); }

   /*
   * Load internal state from archive.
   */
   void LinearCvBias::loadParameters(Serializable::IArchive& ar)
   {  loadParameter<double>(ar, "k", k_); }

   /*
   * Sae internal state to archive.
   */
   void LinearCvBias::save(Serializable::OArchive& ar)
   {  ar << k_; }

   /*
   * Compute and return the bias potential value.
   */
   double LinearCvBias::value(double cv)
   {  return k_*cv; }

   /*
   * Compute and return the derivative of the bias potential.
   */
   double LinearCvBias::derivative(double cv)
   {  return k_; }

}
