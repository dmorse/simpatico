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
   * Classes in the Inter namespace. These interaction classes 
   * define potential energy functions for non-bonded pairs, 
   * bonds, angles, etc. These are all simple non-polymorphic
   * classes that calculate energies and forces for small
   * groups of atoms, e.g., individual pairs, or groups of 3
   * or four atoms for angle dihedral atoms. 
   *
   * The McMd and DdMd namespaces each contain corresponding
   * "potential" classes. The potential classes are polymorphic
   * classes that provide virtual methods to calculate different
   * contributions to the energy and forces for an entire
   * system, as well as four small groups (e.g., pairs). 
   * Implementations of these methods use instances of the
   * "interaction" classes for force and energy calculations
   * in their inner loops. 
   */

   /**
   * \defgroup Inter_Pair_Module Pair Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Classes that represent non-bonded pair interactions.
   */

   /**
   * \defgroup Inter_Bond_Module Bond Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Classes that represent covalent bond interactions.
   */

   /**
   * \defgroup Inter_Angle_Module Angle Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Classes that represent covalent angle interactions.
   */

   /**
   * \defgroup Inter_Dihedral_Module Dihedral Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Classes that represent covalent dihedral interactions.
   */

   /**
   * \defgroup Inter_External_Module External Interactions
   * \ingroup  Inter_Module
   *
   * \brief    Classes that represent external one-body interactions.
   */

}
