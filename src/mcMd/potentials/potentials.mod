
namespace McMd
{

   /**
   * \defgroup McMd_Potential_Module Potential Energies
   * \ingroup  McMd_Module
   *
   * \brief Classes used by a System calculate potential energies.
   *
   * Separate potential energy classes are provided for nonbonded pairs,
   * different types of covalent potentials (pairs, angles, and 
   * dihedrals), and for one-body external potentials. Each type of 
   * potential energy * is represented by an abstract base class, which
   * declares pure virtual methods to compute energies and forces for 
   * an entire System, and a separate implementation class template. 
   * Each implementation template takes an interaction class as a 
   * template argument, and uses an instance of this interaction class 
   * to calculate energies and forces for small groups of atoms (pairs,
   * bonds, etc.). Specific instantiations of each implementation 
   * template, with an interaction class specified in by a style 
   * string in the parameter file, are instantiated by an associated 
   * Factory class.
   *
   * \sa Simp_Module
   */

   /**
   * \defgroup McMd_Pair_Module Pair Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief Nonbonded pair potentials.
   *
   * McPairPotential and MdPairPotential are abstract base classes 
   * that define slightly different interfaces use with MC and MD 
   * simulations, respectively. 
   *
   * In neutral systems (with no Coulomb interactions), these interfaces 
   * are implemented by corresponding class templates MdPairPotentialImpl 
   * and MdPairPotentialImpl, each of which takes an interaction class as 
   * a template argument.  Thus, for example, the class 
   * MdPairPotentialImpl<LJPair> implements a Lennard-Jones pair potential 
   * for use in MD simulations.
   */

   /**
   * \defgroup McMd_Bond_Module Bond Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief Covalent 2-body bond potentials.
   */

   /**
   * \defgroup McMd_Angle_Module Angle Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief Covalent 3-body angle potentials.
   */

   /**
   * \defgroup McMd_Dihedral_Module Dihedral Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief Covalent 4-body dihedral potentials.
   */

   /**
   * \defgroup McMd_Coulomb_Module Coulomb Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief Coulomb potentials, implemented using Ewald separation.
   */

   /**
   * \defgroup McMd_External_Module External Potentials
   * \ingroup McMd_Potential_Module
   *
   * \brief External one-body potentials.
   */

}
