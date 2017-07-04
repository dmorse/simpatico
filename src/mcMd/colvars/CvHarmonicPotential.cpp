/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CvHarmonicPotential.h"

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   CvHarmonicPotential::CvHarmonicPotential() 
   {}

   /*
   * Default destructor.
   */
   CvHarmonicPotential::~CvHarmonicPotential()
   {}

   /*
   * Read parameters from file.
   */
   void CvHarmonicPotential::readParameters(std::istream& in)
   {
      read<double>(in, "cv0", cv0_);
      read<double>(in, "k", k_);
   }

   /*
   * Load internal state from archive.
   */
   void CvHarmonicPotential::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<double>(ar, "cv0", cv0_);
      loadParameter<double>(ar, "k", k_);
   }

   /*
   * Sae internal state to archive.
   */
   void CvHarmonicPotential::save(Serializable::OArchive& ar)
   {
      ar << cv0_;
      ar << k_;
   }

   /*
   * Compute and return the bias potential value.
   */
   double CvHarmonicPotential::value(double cv)
   {
      double dcv = cv - cv0_;
      return 0.5*k_*dcv*dcv;
   }

   /*
   * Compute and return the derivative of the bias potential.
   */
   double CvHarmonicPotential::derivative(double cv)
   {
      return k_*(cv - cv0_);
   }


}
