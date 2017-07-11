/* -*- c++ -*- (enables emacs c++ mode) */
/*==============================================================================
          "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                            Politecnico di Milano
                                A.Y. 2016-2017
                  
                Copyright (C) 2017 Luca Possenti - Simone Di Gregorio
================================================================================*/
/*! 
  @file   problemHT.hpp
  @author Luca Possenti <luca.possenti@polimi.it>
  @author Simone Di Gregorio <simone.digre@gmail.com>
  @date   March 2017.
  @brief  Declaration of the main class for the problem with hematocrit transport.
 */

#ifndef M3D1D_PROBLEM3D1DHT_HPP_
#define M3D1D_PROBLEM3D1DHT_HPP_


#include <problem3d1d.hpp>
#include <assembling1dHT.hpp>
#include <dof1dHT.hpp>
#include <descr1dHT.hpp>
#include <param1dHT.hpp>
#include <mesh1dHT.hpp>
#include <Fahraeus.hpp>


namespace getfem {

//!	Main class defining the coupled 3D/1D fluid problem.
class problemHT: public problem3d1d { 

public: 

	//! Constructor
	/*! 
		It links integration methods and finite element methods to the meshes 
	*/
	problemHT(void):
		mf_coefh(meshv)
	{}
	//! Initialize the problem
	/*!
		1. Read the .param filename from standard input
		2. Import problem descriptors (file paths, GetFEM types, ...)
		3. Import mesh for tissue (3D) and vessel network (1D) + mesh for Ht in vessel network
		4. Set finite elements and integration methods
		5. Build problem parameters
		6. Build the list of tissue boundary data
		7. Build the list of vessel boundary (and junction) data
	 */
	void init (int argc, char *argv[]);
	//! Assemble the problem
	/*!
		1. Initialize problem matrices and vectors
		2. Build the monolithic matrix AM
		3. Build the monolithic rhs FM
	 */
	void assembly (void);
	void assembly_fixpoint (void);
	//! Solve the problem
	/*!
		Solve the monolithic system AM*UM=FM (direct or iterative)
	 */
	bool solve_fixpoint (void);

	//! Export results into vtk files
	/*!
		Export solutions Ut, Pt, Uv, Pv, Ht from the monolithic array UM and HT
	 */
	void export_vtk (const string & suff = "");
	//! Flag to linear or sigmoid lymphatic
	bool HEMATOCRIT_TRANSPORT(int argc, char *argv[]);

protected:
	//! Mesh for the hematocrit in network @f$\Lambda@f$ (1D)
	mesh meshh;
	//! Finite Element Method for the vessel hematocrit @f$H@f$
	//! \note Array of local FEMs on vessel branches  @f$\Lambda_i@f$ (@f$i=1,\dots,N@f$)
	vector<mesh_fem> mf_Hi;
	//! Finite Element Method for PDE coefficients defined on the vessel branches
	//! \note Array of local FEMs on vessel branches  @f$\Lambda_i@f$ (@f$i=1,\dots,N@f$)
	vector<mesh_fem> mf_coefhi;
	//! Finite Element Method for PDE coefficients defined on the network
	mesh_fem mf_coefh;

	//! Input file
	ftool::md_param PARAM;
	//! Algorithm description strings (mesh files, FEM types, solver info, ...) 
	descr1dHT descrHT;
	//! Physical parameters (dimensionless)
	param1dHT paramHT;
	//! Number of degrees of freedom
	dof1dHT  dofHT;
	//! List of BC nodes of the network
	vector< node > BCv_HT;
	//! List of junction nodes of the network
	vector< node > Jv_HT; //?
	
	////////////////////////////////////////////////////////////////////

	//! Monolithic matrix for the coupled problem
	sparse_matrix_type AM_HT;
	//! Monolithic array of unknowns for the coupled problem
	vector_type        UM_HT;		
	//! Monolithic right hand side for the coupled problem
	vector_type        FM_HT;
	//! Monolithic viscosity vector
	vector_type	   MU;

	////////////////////////////////////////////////////////////////////
	
	// Aux methods for init
	//! Import algorithm specifications
	void import_data(void);
	//! Import mesh for tissue (3D) and vessel (1D)  
	void build_mesh(void); 
	//! Set finite elements methods and integration methods 
	void set_im_and_fem(void);
	//! Build problem parameters
	void build_param(void);
	//! Build the list of vessel boundary (and junctions) data 
	void build_vessel_boundary(void);
	// Aux methods for assembly
	//! Build the monolithic matrix AM by blocks
	void assembly_mat(void);
	//! Build the monolithic rhs FM by blocks
	void assembly_rhs(void);
	//! Solve the iterative system of hematocrit problem
	vector_type iteration_solve(vector_type,vector_type);
	scalar_type calcolo_Rk(vector_type, vector_type);



}; /* end of class problem3d1d_HT */


} /* end of namespace */

#endif
