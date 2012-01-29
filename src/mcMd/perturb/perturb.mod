
namespace McMd{

   /**  
   * \defgroup Perturb_Module Free Energy Perturbation
   * \ingroup McMd_Module
   *
   * \brief Classes that implement free energy perturbation theory and
   * replica exchange.
   *
   * This module contains the Perturbation abstract base class and its
   * concrete subclasses. Subclasses of Perturbation describes the 
   * dependence of the Boltzmann weight upon a parameter. The parameter 
   * may be either a parameter in the Hamiltonian (i.e., and interaction 
   * strength) or the inverse temperature beta = 1/(kT) that multiplies 
   * the Hamiltonian. When used in parallel mode, a Perturbation may be 
   * used to assign different values of a perturbation parameter to
   * different processors.
   *
   * This module also contains the ReplicaMove class, which implements 
   * a replica exchange or parallel tempering simulation algorithm. In 
   * replica exchange, different processors simulate systems with 
   * different values of the perturbation parameter, and occasionally 
   * exchange configurations. The ReplicaMove algorithm is always used
   * in conjunction with a concrete subclass of Perturbation, which 
   * defines how the Botzmann weight depends upon a particular tempering 
   * parameter.
   */
  
} 
