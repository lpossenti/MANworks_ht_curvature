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
  @brief  Miscelleanous assembly routines for the 1D hematocrit problem.
 */
/** @defgroup asm Assembly routines */

#ifndef M3D1D_ASSEMBLING_1D_HT_HPP_
#define M3D1D_ASSEMBLING_1D_HT_HPP_
#include <defines.hpp>
#include <node.hpp>
#include <utilities.hpp>
#include <algorithm>
#include <Fahraeus.hpp>

namespace getfem {
//! Build the advection and artificial diffusion matrices for the 1D Hematocrit transport problem,
//! @f$ D = \int_{\Lambda} \frac {\partial (\pi R'^2 u_v H)} {\partial s} w~ds
//! @f$ Df = \int_{\Lambda} D^{art} \nabla^2 H w~ds
/*!
	@param D         Computed advection matrix
	@param Df        Computed artificial diffusion matrix
	@param mim       The integration method to be used
	@param mf_h      The finite element method for the hematocrit @f$ H @f$
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for M
	@param lambdax   First cartesian component of the tangent versor  @f$ \mathbf{\lambda} @f$
	@param lambday   Second cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param lambdaz   Third cartesian component of the tangent versor @f$ \mathbf{\lambda} @f$
	@param rg        The region where to integrate
	@param U	 Velocity Field
	@param R	 Radius of the vessels
	@param diff	 Artificial diffusione coefficient

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_advection_hematocrit
	(MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_h,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_data,
	 const VEC & U, const VEC & R,
	 const VEC & lambdax, const VEC & lambday, const VEC & lambdaz,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 	
	 
	 {
generic_assembly 
	assem("l1=data$1(#2); l2=data$2(#2); l3=data$3(#2); u=data$4(#3);"
		  "t=comp(Grad(#1).Base(#1).Base(#2).Base(#3));"
		  "M$1(#1,#1)+=t(:,1,:,i,p).l1(i).u(p) +t(:,2,:,i,p).l2(i).u(p)+t(:,3,:,i,p).l3(i).u(p);"); //.u(p)


	assem.push_mi(mim);
	assem.push_mf(mf_h);
	assem.push_mf(mf_data);
	assem.push_mf(mf_u);
	assem.push_data(lambdax);
	assem.push_data(lambday);
	assem.push_data(lambdaz);
	assem.push_data(U);
	assem.push_mat(D);
	assem.assembly(rg);


} /*end asm_advection_hematocrit */

template<typename MAT, typename VEC>
void 
asm_network_artificial_diffusion
	(MAT & D,
	 const mesh_im & mim,
	 const mesh_fem & mf_h,
	 const mesh_fem & mf_data,
	 const VEC & diff,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 	
	 
	 {

 getfem::asm_stiffness_matrix_for_laplacian(D,mim,mf_h,mf_data, diff, rg);

} /*end asm_network_artificial_diffusion */

//! Build the junction matrices for the 1D Hematocrit transport problem,
//! --> conservation of mass or Pries formula for birfucation (Pries et al 2005)
/*!
	@param J         Junction matrix for fluid dynamic problem
	@param Jh        Junction matrix for hematocrit problem
	@param mim_u     The integration method to be used
	@param mf_h      The finite element method for the hematocrit @f$ H @f$
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_data_u The finite element method for the tangent versor on @f$ \Lambda @f$
	@param mf_p      The finite element method for the pressure @f$ p @f$
	@param M      	 Global matrix for linear system resolution
	@param J_data    Nodes belonging to each junction
	@param U	 Velocity Field
	@param radius	 Radius of the vessels
	@param dim	 Characteristic dimension of the problem
	@param H_old	 Hematocrit along the vessel in previous iteration

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void
asm_hematocrit_junctions
	(MAT & J, 
	 MAT & Jh,
	 const VEC & U,
	 const mesh_im & mim_u,
	 const std::vector<mesh_fem> & mf_h,
	 const mesh_fem & mf_p,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_data_u,
	 const std::vector<getfem::node> & J_data,
	 const VEC & radius,
	 const VEC & H_old,
	scalar_type dim,
	MAT & M
	 ) 
{
	GMM_ASSERT1 (mf_p.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	GMM_ASSERT1 (mf_h[0].get_qdim() == 1, 
		"invalid data mesh fem for velocity (Qdim=1 required)");
	GMM_ASSERT1 (getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK(1,0)" &&
		getfem::name_of_fem(mf_p.fem_of_element(0)) != "FEM_PK_DISCONTINUOUS(1,0)",
		"invalid data mesh fem for pressure (k>0 required)");

	sparse_matrix_type Diameters(gmm::mat_nrows(J),gmm::mat_ncols(J));
	
	for (size_type i=0; i<mf_h.size(); ++i){ /* branch loop */

		scalar_type Ri = compute_radius(mim_u, mf_data_u, radius, i);
		//cout << "Region " << i << " : radius=" << Ri << endl;

		for (size_type j=0; j<J_data.size(); ++j){

			// Identify pressure dof corresponding to junction node
			VEC psi(mf_p.nb_dof());
			asm_basis_function(psi, mim_u, mf_p, J_data[j].rg);
			size_type row = 0;
			bool found = false;
			while (!found && row<mf_p.nb_dof()){
				found = (1.0 - psi[row] < 1.0E-06);
				if (!found) row++;
			}
			GMM_ASSERT1 (row!=0 && found,  // No junction in first point
				"Error in assembling pressure basis function");
			std::vector<long signed int>::const_iterator bb = J_data[j].branches.begin();
			std::vector<long signed int>::const_iterator be = J_data[j].branches.end();
			size_type last_, first_,last_u,first_u;
			vector_type dof_enum;
			vector_type dofu_enum;
			int fine=0;
			int fine_u=0;
			for (getfem::mr_visitor mrv(mf_h[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto b : mf_h[i].ind_basic_dof_of_element(mrv.cv()))
				{dof_enum.emplace_back(b);
				fine++;}
			for (getfem::mr_visitor mrv(mf_u[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto ub : mf_u[i].ind_basic_dof_of_element(mrv.cv()))
				{dofu_enum.emplace_back(ub);
				fine_u++;}			
			first_=dof_enum[0];
			last_=dof_enum[fine-1];
			dof_enum.clear();
			first_u=dofu_enum[0];
			last_u=dofu_enum[fine_u-1];
			dofu_enum.clear();
			// Outflow branch contribution
			if (std::find(bb, be, i) != be){
				J(row, i*mf_h[i].nb_dof()+last_) -= pi*Ri*Ri*U[i*mf_u[i].nb_dof()+last_u];//col to be generalized!
				Diameters(row, i*mf_h[i].nb_dof()+last_) += 2*Ri*dim;
			}
			// Inflow branch contribution
			if (i!=0 && std::find(bb, be, -i) != be){
				J(row, i*mf_h[i].nb_dof()+first_) += pi*Ri*Ri*U[i*mf_u[i].nb_dof()+first_u];//col to be generalized!
				Diameters(row, i*mf_h[i].nb_dof()+first_) += 2*Ri*dim;
			}
		}
	}

	size_type Ncol = gmm::mat_ncols(J);
	size_type Nrow = gmm::mat_nrows(J);
	vector_type row_vec(Ncol);
	sparse_matrix_type Jq(Nrow,Ncol);
	for(size_type i=0; i<Nrow; i++)
	{
	for(size_type j=0; j<Ncol; j++)
  	gmm::mat_row(Jq,i)[j]= gmm::mat_row(J,i)[j];
	}
	

	for(size_type k=0; k<Ncol; k++)
	{
		for(size_type n=0; n<Nrow; n++)
			{
			if (Jq[n][k]>0)
				{	gmm::copy(gmm::mat_row(Jq,n), row_vec);

				int mycount=count_if(row_vec.begin(), row_vec.end(), isPositive);
					if(mycount==2){
						scalar_type value=0;
						int position=-1;
							while (value >=0)
								{
								position++;
								value=row_vec[position];	
								}
					scalar_type H_f=H_old[position];
					scalar_type FQB=Jq[n][k];
						    FQB=FQB/value*(-1);
					scalar_type D_f=Diameters[n][position];
					scalar_type D=Diameters[n][k];
					scalar_type pos_2=0;
							value=0;
							while (pos_2==k || value <=0)
								{
								pos_2++;
								value=row_vec[pos_2];
								}
					scalar_type D_2=Diameters[n][pos_2];
					scalar_type FQE=fractional_Erythrocytes(FQB, D_f, D, D_2, H_f);
					gmm::clear(row_vec);
					row_vec[position]=Jq[n][position]*FQE;
					gmm::copy(row_vec,gmm::mat_row(Jh,k));
					}
					else{
					row_vec[k]=0;
					gmm::copy(row_vec,gmm::mat_row(Jh,k));
					}
				} //end if
			} // end for of the columns
	} // end for of the elements in a row
	
} /* end of asm_junctions */

template<typename MAT, typename VEC>
void
asm_HT_bc
	(MAT & M, VEC & F,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_h,
	 const mesh_fem & mf_data,
	 const scalar_type beta,
	const std::vector<getfem::node> &  BC, 
	const VEC & radius
	 ) 
{
size_type shift=0;
for (size_type bc=0; bc < BC.size(); bc++) { 
			shift=0;
			size_type i = abs(BC[bc].branches[0]);
			if(i!=0)
			for(size_type j=0; j<i; j++)
			{
			shift = shift+mf_h[j].nb_dof();
			}
			sparse_matrix_type Mi(mf_h[i].nb_dof(),mf_h[i].nb_dof());
			gmm::copy(gmm::sub_matrix(M, 
					gmm::sub_interval(shift, mf_h[i].nb_dof()),
					gmm::sub_interval(shift, mf_h[i].nb_dof())), Mi);
			vector_type Fi(mf_h[i].nb_dof()) ;
			scalar_type Ri=compute_radius( mim, mf_data,radius, i);		
	
		if (BC[bc].label=="DIR") { // Dirichlet BC

			// Add Hin contribution to F -> F(0)=DIR

			vector_type BC_temp(mf_h[i].nb_dof(),BC[bc].value);

			getfem::assembling_Dirichlet_condition(Mi, Fi, mf_h[i], BC[bc].rg, BC_temp);

			gmm::add(Fi, 
			  gmm::sub_vector(F, 
					gmm::sub_interval(shift, mf_h[i].nb_dof())));
			gmm::clear(gmm::sub_matrix(M, 
					gmm::sub_interval(shift, mf_h[i].nb_dof()),
					gmm::sub_interval(shift, mf_h[i].nb_dof())));
			gmm::copy(Mi, 
			  gmm::sub_matrix(M, 
					gmm::sub_interval(shift, mf_h[i].nb_dof()),
					gmm::sub_interval(shift, mf_h[i].nb_dof()))); 
			gmm::clear(BC_temp);

				
		} /*end DIR condition*/
		else if (BC[bc].label=="MIX") { // Robin BC


			MAT Mi_mix(mf_h[i].nb_dof(), mf_h[i].nb_dof());
			VEC BETA(mf_data.nb_dof(), beta*pi*Ri*Ri);
			getfem::asm_mass_matrix_param(Mi_mix,
				mim, mf_h[i],mf_data, BETA, BC[bc].rg);
			
			gmm::add(gmm::scaled(Mi_mix, 1.0), 
				gmm::sub_matrix(M,
					gmm::sub_interval(shift, mf_h[i].nb_dof()),
					gmm::sub_interval(shift, mf_h[i].nb_dof())));
			gmm::clear(Mi_mix);

			// Add p0 contribution to F
			vector_type BC_temp_mix(mf_data.nb_dof(),beta*pi*Ri*Ri*BC[bc].value);

			getfem::asm_source_term(gmm::sub_vector(F, gmm::sub_interval(shift,mf_h[i].nb_dof())), 
				mim, mf_h[i], mf_data, BC_temp_mix);
			gmm::clear(BC_temp_mix);
	
		}
		else if (BC[bc].label=="OUT"){
		}
		else if (BC[bc].label=="INT") { // Internal Node
			DAL_WARNING1("internal node passed as boundary.");
		}
		else if (BC[bc].label=="JUN") { // Junction Node
			DAL_WARNING1("junction node passed as boundary.");
		}
		else {
			GMM_ASSERT1(0, "Unknown Boundary Condition"<< BC[bc].label << endl);
		}

	} /*end of for cicle*/
}/* end of asm_HT_bc*/

template<typename MAT, typename VEC>
void
asm_HT_out
	(MAT & M,
	 const mesh_im & mim,
	 const std::vector<mesh_fem> & mf_h,
	const VEC & U, const VEC & radius,
	 const std::vector<mesh_fem> & mf_u,
	 const mesh_fem & mf_data_u
	) 
{
size_type shift=0;
size_type shift_u=0;
for (size_type i=0; i < mf_h.size(); i++) { 
			if(i!=0){
			shift = shift+mf_h[i-1].nb_dof();
			shift_u=shift_u+mf_u[i-1].nb_dof();
			}

			scalar_type Ri=compute_radius( mim, mf_data_u,radius, i);		
			size_type fine_u=0, fine=0;
			vector_type dofu_enum; gmm::clear(dofu_enum);
			vector_type dof_enum; gmm::clear(dof_enum);

			for (getfem::mr_visitor mrv(mf_u[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto ub : mf_u[i].ind_basic_dof_of_element(mrv.cv()))
			{dofu_enum.emplace_back(ub);
				fine_u++;}
			for (getfem::mr_visitor mrv(mf_h[i].linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto b : mf_h[i].ind_basic_dof_of_element(mrv.cv()))
				{dof_enum.emplace_back(b);
				fine++;}

			size_type last_u=dofu_enum[fine_u-1];
			size_type last=dof_enum[fine-1];
			size_type first_u=dofu_enum[0];
			size_type first=dof_enum[0];
		if (U[i*mf_u[i].nb_dof()+last_u]>0)
	M[i*mf_h[i].nb_dof()+last][i*mf_h[i].nb_dof()+last]+=pi*Ri*Ri*U[i*mf_u[i].nb_dof()+last_u];	else
	M[i*mf_h[i].nb_dof()+first][i*mf_h[i].nb_dof()+first]-=pi*Ri*Ri*U[i*mf_u[i].nb_dof()+first_u];

	} /*end of for cicle*/


}/* end of asm_HT_out*/

//! Build the mass matrices for the 1D Poiseuille's problem updated in iterations
//! @f$ M = \int_{\Lambda} c~u~v~ds @f$ and
/*!
	@param M         Computed mass matrix
	@param mim       The integration method to be used
	@param mf_u      The finite element method for the velocity @f$ \mathbf{u} @f$
	@param mf_data   The finite element method for the tangent versor on @f$ \Lambda @f$
	@param coef      The coefficient for M
	@param rg        The region where to integrate

	@ingroup asm
 */ 
template<typename MAT, typename VEC>
void 
asm_network_poiseuilleHT
	(MAT & M,
	 const mesh_im & mim,
	 const mesh_fem & mf_u,
	 const mesh_fem & mf_data,
	 const VEC & coef,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	// Build the local mass matrix Mvvi
	getfem::asm_mass_matrix_param(M, mim, mf_u, mf_data, coef, rg);

}

} /* end of namespace */

#endif
