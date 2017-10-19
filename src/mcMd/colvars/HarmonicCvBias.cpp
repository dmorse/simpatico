/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "HarmonicCvBias.h"

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   HarmonicCvBias::HarmonicCvBias() 
   {}

   /*
   * Default destructor.
   */
   HarmonicCvBias::~HarmonicCvBias()
   {}

   /*
   * Read parameters from file.
   */
   void HarmonicCvBias::readParameters(std::istream& in)
   {
      read<double>(in, "cv0", cv0_);
      read<double>(in, "k", k_);
   }

   /*
   * Load internal state from archive.
   */
   void HarmonicCvBias::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<double>(ar, "cv0", cv0_);
      loadParameter<double>(ar, "k", k_);
   }

   /*
   * Sae internal state to archive.
   */
   void HarmonicCvBias::save(Serializable::OArchive& ar)
   {
      ar << cv0_;
      ar << k_;
   }

   /*
   * Compute and return the bias potential value.
   */
   double HarmonicCvBias::value(double cv)
   {
      double dcv = cv - cv0_;
      return 0.5*k_*dcv*dcv;
   }

   /*
   * Compute and return the derivative of the bias potential.
   */
   double HarmonicCvBias::derivative(double cv)
   {
      return k_*(cv - cv0_);
   }


}
