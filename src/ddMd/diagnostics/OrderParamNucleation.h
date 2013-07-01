#ifndef DDMD_ORDER_PARAM_NUCLEATION_H
#define DDMD_ORDER_PARAM_NUCLEATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <ddMd/diagnostics/Diagnostic.h>
#include <ddMd/simulation/Simulation.h>
#include <util/containers/DMatrix.h>               // member template

#include <util/global.h>

#include <iostream>

namespace DdMd
{

   using namespace Util;

   /**
   * \ingroup DdMd_Diagnostic_Module
   *
   * OrderParamNucleation is designed to evaluate a sort of 
   * composition profile in the specific case of nucleating a
   * lamellar ordered phase in a small region at the center of
   * a simulation cell. Composition profile is calculated in
   * terma of a quantity called cos factor which is simple the 
   * cosine component of a typical structure factor calculation. 
   * The cos factor is evaluated per atom type and for the 
   * specific wavevector corresponding to the direction 
   * perpendicular to lamellae (specified by the parameter 
   * perpIndex).
   *
   * The cos factors are accumulated in different bins based on 
   * the atom position in the direction parallel to lamellae
   * (specified by the parameter parallelIndex).
   *
   * \code
   * OrderParamNucleation{
   *    interval                      1000
   *    outputFileName    OrderParamNucleation
   *    perpIndex                        2
   *    parallelIndex                    0
   *    periodicity                      3
   *    nBin                           100
   * }
   * \endcode
   *
   * \ingroup DdMd_Diagnostic_Module
   */
   class OrderParamNucleation : public Diagnostic
   {

   public:

      /**	
      * Constructor.
      *
      * \param simulation reference to parent Simulation object
      */
      OrderParamNucleation(Simulation& simulation);

      /**	
      * Destructor.
      */
      ~OrderParamNucleation();

      /**
      * Read parameters from file.
      *
      * \param in input parameter stream
      */
      virtual void readParameters(std::istream& in);

      /**
      * Load internal state from an archive.
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
      * Clear accumulators.
      */
      virtual void clear();
   
      /**
      * Add particles to OrderParamNucleation accumulators.
      *
      * \param iStep step counter
      */
      void sample(long iStep);

      /**
      * Output results to predefined output file.
      */
      virtual void output();

   protected:

      /// Output file stream.
      std::ofstream outputFile_;

      /// Number of atom types, copied from Simulation::nAtomType().
      int  nAtomType_;

      /// Cos factors of size nAtomType*nBin
      DMatrix<double> cosFactors_;

      /// Total cos factors of size nAtomType*nBin (gathered from all processors)
      DMatrix<double> totalCosFactors_;
      
      /// Index of direction perpendicular to lamellae (direction of external field)
      int  perpIndex_;

      /// Index of direction parallel to the lamellar surfaces
      int parallelIndex_;

      /// Number of periods in the perpendicular direction
      int periodicity_;

      /// Number of bins to divide the length along parallel direction 
      int nBin_;

      /// Number of samples thus far.
      int  nSample_;

      /// Has readParam been called?
      bool isInitialized_;

   };

}
#endif
