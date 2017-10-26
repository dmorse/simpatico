#ifndef MCMD_LINEAR_CV_BIAS_H
#define MCMD_LINEAR_CV_BIAS_H

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
   * Linear function of a collective variable (CV).
   *
   * This defines a simple linear function of a collective variable,
   * of the form  k*cv, where cv is a current value of a collective 
   * and k is a parameter.
   *
   * \ingroup McMd_Colvar_Module
   */
   class LinearCvBias : public CvBias
   {

   public:

      /**
      * Constructor.
      */
      LinearCvBias();

      /**
      * Destructor.
      *
      * Empty default implementation.
      */
      virtual ~LinearCvBias();

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
      * Returns k*cv
      *
      * \param cv collective variable value. 
      */
      virtual double value(double cv);

      /**
      * Compute and return the derivative of the bias potential.
      *
      * Returns value of k.
      *
      * \param cv collective variable value. 
      */
      virtual double derivative(double cv);

   private:

      // Prefactor
      double k_;

   };

}
#endif
