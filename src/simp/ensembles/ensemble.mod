namespace Util
{

/**
* \defgroup Ensemble_Module Statistical Ensembles
*
* Statistical mechanical ensembles for specific macro-variables.
*
* The classes EnergyEnsemble, BoundaryEnsemble, and SpeciesEnsemble 
* represent statistical ensembles for fluctuations of energy, boundary 
* volume or shape, and number of molecules of a single species, 
* respectively.  Each ensemble has a type, which specifies whether the
* variable has a constant value or if it has a Boltzman-like weight,
* and stores whatever parameters are needed to describe a Boltzmann 
* ensemble.  
*
* The type of a EnergyEnsemble can be "adiabatic" or "isothermal", and 
* a value of temperature is stored if the type is isothermal. 

* A BoundaryEnsemble can be rigid or isobaric, and stores a pressure 
* value if the type is isobaric. 
*
* A Species ensemble can be "closed" or "grand", and stores a value for 
* the chemical potential of a single species if it is grand.
*
* We envision generalizing the BoundaryEnsemble to allow for constant
* stress ensembles, and may generalize some or all of the ensembles 
* to allow for non-Boltzmann statistical ensembles. 
*
* \ingroup Util_NS_Module
*/

}

