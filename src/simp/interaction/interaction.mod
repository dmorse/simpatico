namespace Simp
{

   /**
   * \defgroup Simp_Interaction_Module Interaction Module
   * \ingroup Simp_Module
   *
   * \brief Potential energy functions.
   *
   * This module contains "interaction" classes that define potential
   * energy functions for individual nonbonded or bonded pair 
   * interactions, 3-body angle groups, 4-body dihedral groups, etc. 
   * These interaction classes are all simple non-polymorphic classes
   * that have no access to information about a parent system.
   *
   * Note: The McMd and DdMd namespaces each contain corresponding
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
