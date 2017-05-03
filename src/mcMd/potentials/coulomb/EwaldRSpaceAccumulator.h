#ifndef MD_EWALD_RSPACE_ACCUMULATOR_h
#define MD_EWALD_RSPACE_ACCUMULATOR_h 


#include <util/space/Tensor.h>
#include <util/misc/Setable.h>

namespace McMd
{

  class PairPotential;

  using namespace Util;

  /**
  * Utility class to store r-space Coulomb energy and stress.
  *
  * \ingroup McMd_Coulomb_Module
  */
  class EwaldRSpaceAccumulator 
  {
  public:

     /*
     * Constructor.
     */
     EwaldRSpaceAccumulator()
        :rSpaceEnergy_(),
         rSpaceStress_(),
         pairPotentialPtr_(0)
     {
        rSpaceEnergy_.unset();
        rSpaceStress_.unset();
     }

     void setPairPotential(PairPotential& potential)
     {  pairPotentialPtr_ = &potential; }

     /**
     * Is the r-space energy set?
     */
     bool isSetEnergy() const
     {  return rSpaceEnergy_.isSet(); }

     /**
     * Is the r-space stress set?
     */
     bool isSetStress() const
     { return rSpaceStress_.isSet(); }

     /**
     * Return the r-space energy (compute if necessary).
     */
     double rSpaceEnergy();

     /**
     * Return the r-space stress (compute if necessary).
     */
     Tensor rSpaceStress();
     
  private:

     Setable<double> rSpaceEnergy_;

     Setable<Tensor> rSpaceStress_;

     PairPotential* pairPotentialPtr_;

  // friend:

     template <typename T>
     friend class MdEwaldPairPotentialImpl;
                  
  };
}
#endif
