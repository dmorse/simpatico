#ifdef  MCMD_PERTURB
#ifndef MCMD_PERTURBATION_H
#define MCMD_PERTURBATION_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2012, David Morse (morse012@umn.edu)
* Distributed under the terms of the GNU General Public License.
*/

#include <util/param/ParamComposite.h>  // Base class
#include <util/containers/DMatrix.h>     // member
#include <util/containers/DArray.h>     // member


namespace McMd
{

   using namespace Util;

   /**
   * Model of parameter dependence in a free energy perturbation theory.
   *
   * A Perturbation object defines a parameter dependent statistical weight
   * of the type used in free energy perturbation theory, replica exchange,
   * and extended ensemble algorithms. In all of these applications, the
   * probability of a microstate X is given by an exponential exp(-W(X,p)),
   * where W(X,p) = beta*H depends upon a perturbation parameter p, H is the
   * system Hamiltonian and beta = 1/kT is inverse temperature. Either the
   * Hamiltonian or the inverse temperature may depend upon the perturbation
   * parameter.
   *
   * This class provides a method to modify the system parameters (beta
   * and/or the potential energy parameters) so as to correspond to a
   * given value of p, and methods to evaluate the derivative dW(X,p)/dp
   * and difference W(X, p') - W(X, p) for a given system configuration.
   *
   * \ingroup McMd_Perturb_Module
   */
   class Perturbation : public ParamComposite
   {

   public:

      /**
      * Constructor.
      *  
      * \param size number of parameter values (e.g., communicator size)
      * \param rank index of this system (e.g., communicator rank)
      */
      Perturbation(int size, int rank);

      /**
      * Destructor.
      */
      virtual ~Perturbation();

      /// \name Mutators
      //@{

      /**
      * Read perturbation parameter(s) from file.
      *
      * \param in input stream (file or std in).
      */ 
      void readParameters(std::istream& in);
 
      /**
      * Load internal state from an archive.
      *
      * \param ar input/loading archive
      */
      virtual void loadParameters(Serializable::IArchive& ar);

      /**
      * Save internal state to an archive.
      *
      * \param ar output/saving archive
      */
      virtual void save(Serializable::OArchive& ar);

      /**
      * Sets the perturbation parameter value.
      *
      * This method sets the value of the member parameter, and calls the
      * virtual void setParameter() method to modify the parameter of the
      * associated System.
      *
      * \param parameter perturbation parameter value for associated System.
      */
      void setParameter(DArray<double> parameter);
      
      //@}
      /// \name Accessors
      //@{
     
      /**
      * Get parameter i of system id.
      *
      * \param i index of perturbation parameter
      * \param id id of system in a mpi simulation
      *
      * Returns the value of the perturbation parameter p[i].
      */
      double parameter(int i, int id);
 
      /**
      * Gets the number of parameters per system.
      *
      * Returns the number of perturbation parameters for this system.
      */ 
      int getNParameters() const;

      /**
      * Get perturbation parameter p[i] of index i for this system.
      *
      * Returns the value of the perturbation parameter p[i].
      */
      virtual double parameter(int i) const = 0;

      /**
      * Get derivative of the Boltzmann weight W(X,p) with respect to the
      * perturbation parameter p[i] of index i.
      *
      * For a system in which W(X,p) = beta*H depends on a microstate X and
      * a parameter p, this method returns the value of dW(X, p)/dp[i] for the
      * current system microstate (i.e., for the current particle positions).
      */
      virtual double derivative(int i) const = 0;

      /**
      * Returns the difference W(X, p') - W(X, p).
      *
      * This method returns the difference W(X, p') - W(X, p) between 
      * value of the statistical weight W(X,p[i]) = H/kT for the current 
      * microstate, where p[i] is the value of the perturbation parameter in 
      * this system, and p'[i] = iPartnerParameter is the partner perturbation 
      * parameter value that is passed to this function.
      *
      * \param iPartnerParameter perturbation parameters for partner system.
      */
      virtual double difference(DArray<double> iPartnerParameter) const = 0;

      //@}

   protected:

      /**
      * Number of systems (e.g., communicator size)  
      */
      int size_;

      /**
      * Index for this system (e.g., communicator size)  
      */
      int rank_;

      /**
      * Number of perturbation parameters associated with a System.
      */
      int nParameters_;

      /**
      * mode 0: parameters of all replica systems are specified.
      * mode 1: parameters of first and last replica systems are specified.
      */
      int mode_;

      /**
      * Value of the perturbation parameter for the associated System. 
      * parameter is a DMatrix of dimensions 1xnParameters.
      */
      DArray<double> parameter_;
      
      /**
      * Value of the perturbation parameter for the first replica System.
      * initialParameter is a DMatrix of dimensions 1xnParameters.
      */
      DArray<double> initialParameter_;
      
      /**
      * Value of the perturbation parameter for the last replica System.
      * finalParameter is a DMatrix of dimensions 1xnParameters.
      */
      DArray<double> finalParameter_;

      /**
      * Value of the perturbation parameter for all the replica Systems.
      * parameters is a DMatrix of dimensions nProcsxnParameters.
      */
      DMatrix<double> parameters_;

      /**
      * Sets the perturbation parameter in the associated system.
      *
      * This method modifies the parameters of the associated System so as
      * to correspond to the value of the parameter member variable.
      */
      virtual void setParameter() = 0;

   };

}

#endif  // ifndef PERTURBATION_H
#endif  // ifdef  MCMD_PERTURB
