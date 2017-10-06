#ifndef MCMD_HARMONIC_CV_BIAS_H
#define MCMD_HARMONIC_CV_BIAS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

#include <mcMd/colvars/CvBias.h>

namespace McMd
{

   using namespace Util;

   /**
   * Harmonic function of a collective variable (CV).
   *
   * This defines a harmonic function of a collective variable, of the 
   * type often used in umbrella sampling. The potential is defined by 
   * a function  k*(cv - cv0)/2, where cv is a current value of a 
   * collective variable and where cv0 and k are parameters.
   *
   * \ingroup McMd_Colvar_Module
   */
   class HarmonicCvBias : public CvBias
   {

   public:

      /**
      * Constructor.
      */
      HarmonicCvBias();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~HarmonicCvBias();

      /**
      * Read parameters from file.
      * 
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
      * 
      * \param ar input archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      * 
      * \param ar input archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Compute and return the bias potential value.
      *
      * Returns (1/2)*k*(cv - cv0)^2
      *
      * \param cv collective variable value. 
      */
      virtual double value(double cv);

      /**
      * Compute and return the derivative of the bias potential.
      *
      * Returns value of expression k*(cv - cv0)
      *
      * \param cv collective variable value. 
      */
      virtual double derivative(double cv);

   private:

      // Central value of collective variable.
      double cv0_;

      // Spring constant
      double k_;

   };

}
#endif
