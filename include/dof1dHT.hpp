/* -*- c++ -*- (enables emacs c++ mode) */
/*==============================================================================
          "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                            Politecnico di Milano
                                A.Y. 2016-2017
                  
                Copyright (C) 2017 Luca Possenti - Simone Di Gregorio
================================================================================*/
/*! 
  @file   dof1dHT.hpp
  @author Luca Possenti <luca.possenti@polimi.it>
  @author Simone Di Gregorio <simone.digre@gmail.com>
  @date   March 2017.
  @brief  Definition of the aux class for the number of degrees of freedom.
 */
#ifndef M3D1D_DOF1DHT_HPP_
#define M3D1D_DOF1DHT_HPP_

namespace getfem {

//! Class to store the number of degrees of freedom of used FEMs
struct dof1dHT {

	//! Number of dof of the vessel network hematocrit
	//! It is the sum of local vessel branch dof
	size_type H_;
	//! Number of dof of the vessel network hematocrit coefficients FEM mf_coefh
	//! It is NOT the sum of local vessel branch dof
	size_type h_;
	//! Total number of dof of the hematocrit problem
	size_type hematocrit_;
	
	//! Compute the number of dof of given FEM
	void set (
			const std::vector<getfem::mesh_fem> & mf_H,
			const getfem::mesh_fem & mf_coefh
			)
	{
		H_ = 0;
		for (size_type i = 0; i<mf_H.size(); ++i) H_ += mf_H[i].nb_dof();
		h_ = mf_coefh.nb_dof();
		hematocrit_ = H_;
	}
	//! Accessor to the number of dof of mf_H
	inline size_type H (void) { return H_; } const
	//! Accessor to the number of dof of mf_coefh
	inline size_type coefh (void) { return h_; } const
	//! Accessor to the number of dof of vessel problem
	inline size_type hematocrit (void) { return hematocrit_; } const
	//! Accessor to the number of dof of coupled problem
	inline size_type tot (void) { return hematocrit_; } const

	//! Overloading of the output operator
	friend std::ostream & operator << (
			std::ostream & out, const dof1dHT & dof
			)
	{ 
		out << "--- DEGREES OF FREEDOM OF HEMATOCRIT PROBLEM --- " << endl;
		out << "  nb_dof_H      : " 		 << dof.H_  << endl;
		out << "  nb_dof_coefh  : " 		 << dof.h_  << endl;
		out << "  nb_dof_tot    : " 		 << dof.hematocrit_ << endl;
		out << "-------------------------- " << endl;
		return out;            
	}

};

} /* end of namespace */

#endif
