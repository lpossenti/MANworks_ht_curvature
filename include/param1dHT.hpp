/* -*- c++ -*- (enables emacs c++ mode) */
/*==============================================================================
          "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                            Politecnico di Milano
                                A.Y. 2016-2017
                  
                Copyright (C) 2017 Luca Possenti - Simone Di Gregorio
================================================================================*/
/*! 
  	@file   assembling1dHT.hpp
  	@author Luca Possenti <luca.possenti@polimi.it>
  	@author Simone Di Gregorio <simone.digre@gmail.com>
  	@date   March 2017.
	@brief  Definition of the aux class for physical parameters.
	@details 
	Assemble the dimensionless parameters of the coupled 3D/1D model:
	- Radius @f$R'(s)@f$,
	- Tissue permeability @f$\kappa_t@f$,
	- Vessel wall permeability @f$Q(s)@f$,
	- Vessel bed permeability @f$\kappa_v(s)@f$.

	being @f$s\in\Lambda@f$ the arc-lenght over the vessel network.
	\note @f$\kappa_t@f$ is assumed to be constant.

	\ingroup input
 */
 
#ifndef M3D1D_PARAM1DHT_HPP_
#define M3D1D_PARAM1DHT_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct param1dHT {

	// Dimensional physical parameters (microcirc applications)
	//! Temperature of patient [°C]
	scalar_type Temp_;
	//! Viscosity of the plasma [kg/ms]
	scalar_type mu_plasma_;
	// Dimensionless physical parameters (test-cases)
	
	// Utils
	//! Flag to choose the viscosity computation method:
	int Visco_v; //  0-> Constant choosen by user; 1-> Pries in vivo; 2-> Pries in vitro;

	ftool::md_param FILE_;

	//! Finite Element Method for hematocrit data
	getfem::mesh_fem mf_datah_;
	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname, 
			const getfem::mesh_fem & mf_datah
			) 
	{
		FILE_ = fname;
		mf_datah_ = mf_datah;
		size_type dof_datah = mf_datah_.nb_dof();
		 

		Temp_ = FILE_.real_value("Temp", "blood Temperature [°C]");
		Visco_v = FILE_.real_value("Visco_v", "Flag to choose Pries Viscosity [0] = vivo, [1] = vitro"); 

		// Compute Plasma Viscosity (Biomachines course formula)
		mu_plasma_=viscosity_plasma(Temp_);


	}
	inline scalar_type visco_plasma (void) { return mu_plasma_;  } const
	inline int visco_type (void) { return Visco_v;  } const
	
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const param1dHT & param
		)
	{ 
		out << "--- PHYSICAL PARAMS OF HEMATOCRIT PROBLEM ---" << endl;
		out << "--------------------------------" << endl;

		return out;            
	}
}; /* end of class */

} /* end of namespace */

#endif
