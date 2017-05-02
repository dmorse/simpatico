/**
* Potential energy interaction functions.
*/
namespace Inter {}

namespace Inter
{

   /**
   * \defgroup Inter_Module Inter namespace
   *
   * \brief Potential energy functions.
   *
   * The interaction classes in the Inter namespace define potential 
   * energy functions for individual nonbonded or bonded pairs, 
   * 3-body angle groups, 4-body dihedral groups, etc. These
   * interaction classes are all simple non-polymorphic classes
   * that have no access to information about a parent system.
   *
   * The McMd and DdMd namespaces each contain corresponding
   * "potential" classes that provide methods to calculate energies
   * and forces for an entire system. The implementation of each
   * such "potential" class uses an associated "interaction" class 
   * for core calculations in the inner force and energy loops.
   *
   * See also:
   *
   * \ref McMd_Potential_Module "McMd Potential Module"
   *
   * \ref DdMd_Potential_Module "DdMd Potential Module"
   */

}
