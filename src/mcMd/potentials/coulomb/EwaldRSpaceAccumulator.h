#ifndef MD_EWALD_RSPACE_ACCUMULATOR_h
#define MD_EWALD_RSPACE_ACCUMULATOR_h 


#include <util/space/Tensor.h>
#include <util/misc/Setable.h>

namespace McMd
{

  using namespace Util;

  class EwaldRSpaceAccumulator 
  {
  public:

     /*
     * Constructor.
     */
     EwaldRSpaceAccumulator()
        :rSpaceEnergy_(),
         rSpaceStress_()
     {
        rSpaceEnergy_.unset();
        rSpaceStress_.unset();
     }

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
     * Return the r-space energy.
     */
     double rSpaceEnergy() const
     { return rSpaceEnergy_.value(); }

     /**
     * Return the r-space stress.
     */
     Tensor rSpaceStress() const
     { return rSpaceStress_.value(); }
     
     //Tensor rSpacePressure() const;

  private:

     Setable<double> rSpaceEnergy_;
     Setable<Tensor> rSpaceStress_;

  // friend:

     template <typename T>
     friend class MdEwaldPairPotentialImpl;
                  
  };
}
#endif
