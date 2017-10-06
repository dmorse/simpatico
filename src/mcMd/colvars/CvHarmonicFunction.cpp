/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include "CvHarmonicFunction.h"

namespace McMd
{

   using namespace Util;

   /*
   * Default constructor.
   */
   CvHarmonicFunction::CvHarmonicFunction() 
   {}

   /*
   * Default destructor.
   */
   CvHarmonicFunction::~CvHarmonicFunction()
   {}

   /*
   * Read parameters from file.
   */
   void CvHarmonicFunction::readParameters(std::istream& in)
   {
      read<double>(in, "cv0", cv0_);
      read<double>(in, "k", k_);
   }

   /*
   * Load internal state from archive.
   */
   void CvHarmonicFunction::loadParameters(Serializable::IArchive& ar)
   {
      loadParameter<double>(ar, "cv0", cv0_);
      loadParameter<double>(ar, "k", k_);
   }

   /*
   * Sae internal state to archive.
   */
   void CvHarmonicFunction::save(Serializable::OArchive& ar)
   {
      ar << cv0_;
      ar << k_;
   }

   /*
   * Compute and return the bias potential value.
   */
   double CvHarmonicFunction::value(double cv)
   {
      double dcv = cv - cv0_;
      return 0.5*k_*dcv*dcv;
   }

   /*
   * Compute and return the derivative of the bias potential.
   */
   double CvHarmonicFunction::derivative(double cv)
   {
      return k_*(cv - cv0_);
   }


}
