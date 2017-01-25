#ifndef MCMD_MD_POTENTIALS_H
#define MCMD_MD_POTENTIALS_H

/*
* Simpatico - Simulation Package for Polymeric and Molecular Liquids
*
* Copyright 2010 - 2017, The Regents of the University of Minnesota
* Distributed under the terms of the GNU General Public License.
*/

/**
* \file md_potentials.h
*
* \brief Include this file to include all MD potentials at once.
*/

#ifndef INTER_NOPAIR
#include <mcMd/potentials/pair/MdPairPotential.h>
#endif
#include <mcMd/potentials/bond/BondPotential.h>
#ifdef INTER_ANGLE
#include <mcMd/potentials/angle/AnglePotential.h>
#endif
#ifdef INTER_DIHEDRAL
#include <mcMd/potentials/dihedral/DihedralPotential.h>
#endif
#ifdef INTER_EXTERNAL
#include <mcMd/potentials/external/ExternalPotential.h>
#endif
#ifdef INTER_TETHER
#include <mcMd/potentials/tether/MdTetherPotential.h>
#endif

#endif
