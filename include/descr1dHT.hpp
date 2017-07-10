/* -*- c++ -*- (enables emacs c++ mode) */
/*==============================================================================
          "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                            Politecnico di Milano
                                A.Y. 2016-2017
                  
                Copyright (C) 2017 Luca Possenti - Simone Di Gregorio
================================================================================*/
/*! 
  @file   descr1dHT.hpp
  @author Luca Possenti <luca.possenti@polimi.it>
  @author Simone Di Gregorio <simone.digre@gmail.com>
  @date   March 2017.
  @brief  Definition of the aux class for algorithm description strings.
/** @defgroup input User-defined parameters  */

#ifndef M3D1D_DESCR1DHT_HPP_
#define M3D1D_DESCR1DHT_HPP_

#include <string>
#include <iostream>

namespace getfem {

//! Class to import the descriptors of the coupled 3D/1D solver
/*!
	\ingroup input
 */
struct descr1dHT {

	//! Absolute path to the vessel mesh file
	std::string MESH_FILEH;
	//! Identifier of vessel mesh tipe
	std::string MESH_TYPEH;
	//! Identifier of vessel hematocrit FEM type
	std::string FEM_TYPEH;
	//! Identifier of vessel coefficients' FEM type
	std::string FEM_TYPEH_DATA;
	//! Identifier of vessel integration method type
	std::string IM_TYPEH;
	//! Maximum residual of hematocrit (Fixed Point Method)
	scalar_type epsH; 
	//! Under Relaxation Coefficient
	scalar_type underH;
	//! Flag to have hematocrit transport along the vessels
	bool HEMATOCRIT_TRANS; 
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Import algorithm specifications from file .param
	void import(ftool::md_param & fname) 
	{
		FILE_ = fname;
		HEMATOCRIT_TRANS = FILE_.int_value("HEMATOCRIT_TRANSPORT", "Flag to compute the transport of hematocrit");
		if(HEMATOCRIT_TRANS){
                epsH = FILE_.real_value("Residual_Hema_FPM", "Max Residuals for Hematocrit (FPM)");
		underH = FILE_.real_value("UNDER_RELAXATION_COEFFICIENT_HEMA", "Under-relaxation coefficient for hematocrit [-]");
		MESH_FILEH  = FILE_.string_value("MESH_FILEH","1D points file");
		MESH_TYPEH  = FILE_.string_value("MESH_TYPEV","1D mesh type");
		FEM_TYPEH   = FILE_.string_value("FEM_TYPEH","FEM 1D vessel - velocity");
		FEM_TYPEH_DATA = FILE_.string_value("FEM_TYPEH_DATA");
		IM_TYPEH 	= FILE_.string_value("IM_TYPEH","Name of integration method");
		if(FEM_TYPEH_DATA=="") FEM_TYPEH_DATA = "FEM_PK(1,0)";}
	}

	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const descr1dHT & descr
		)
	{ 
		cout << "---- PROBLEM DESCRIPTORS FOR HEMATOCRIT----------------" << endl;
		cout << " MESH FILE 1D problem      : " << descr.MESH_FILEH  << endl;
		cout << " IM  TYPE  1D problem      : " << descr.IM_TYPEH	<< endl;
		cout << " FEM TYPE  1D velocity     : " << descr.FEM_TYPEH   << endl;
		cout << " FEM TYPE  1D coefficients : " << descr.FEM_TYPEH_DATA << endl;
		cout << "--------------------------------------------------" << endl;

		return out;            
	}

}; /* end of class */

} /* end of namespace */

#endif
