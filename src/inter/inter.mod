/**
* Classes that define potential energy interaction functions.
*/
namespace Inter {}

namespace Inter
{

   /**
   * \defgroup Inter_Module Inter namespace
   *
   * \brief    Classes that define potential energy functions.
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

   /**
   * \defgroup Inter_Pair_Module Pair Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Potential energy functions for non-bonded pair interactions.
   */

   /**
   * \defgroup Inter_Bond_Module Bond Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Potential energy functions for covalent bond interactions.
   */

   /**
   * \defgroup Inter_Angle_Module Angle Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Potential energy functions for covalent angle interactions.
   */

   /**
   * \defgroup Inter_Dihedral_Module Dihedral Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Potential energy functions for covalent dihedral interactions.
   */

   /**
   * \defgroup Inter_External_Module External Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Potential energy functions for external one-body interactions.
   */

}
