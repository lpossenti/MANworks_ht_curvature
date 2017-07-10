/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Simone Di Gregorio
======================================================================*/
/*! 
  @file   	Fahraeus.cpp
  @authors 	Luca Possenti <luca.possenti@polimi.it>
%		Simone Di Gregorio <simone.digre@gmail.com>
  @date   	March 2017
  @brief  	Declaration of Fahraeus effects routines.
*/

#ifndef M3D1D_FAHRAEUS_HPP_
#define M3D1D_FAHRAEUS_HPP_

#include <getfem/getfem_assembling.h> 
#include <math.h>

namespace getfem {

scalar_type
viscosity_plasma(scalar_type temperature);

scalar_type 
viscosity_vitro(scalar_type hematocrit, scalar_type R, scalar_type mu_plasma);

scalar_type 
viscosity_vivo(scalar_type hematocrit, scalar_type R, scalar_type mu_plasma);

scalar_type
fractional_Erythrocytes(scalar_type, scalar_type, scalar_type, scalar_type, scalar_type);

scalar_type
PS_compute_A(scalar_type,scalar_type, scalar_type, scalar_type);

scalar_type
PS_compute_B(scalar_type, scalar_type);

scalar_type
PS_compute_X0(scalar_type, scalar_type);

scalar_type logit(scalar_type);

bool isPositive (scalar_type);


}

#endif
