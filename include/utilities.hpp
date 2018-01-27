/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
  @file   utilities.hpp
  @author Domenico Notaro <domenico.not@gmail.com>
  @date   January 2016.
  @brief  Declaration of some miscelleanous auxiliary routines.
 */
#ifndef M3D1D_UTILITIES_HPP_
#define M3D1D_UTILITIES_HPP_

#include <getfem/getfem_assembling.h> 
#include <getfem/getfem_mesh.h>
#include <gmm/gmm.h>
#include <defines.hpp>

namespace getfem {

//! Aux function to compute the diameter of an element
scalar_type 
estimate_h(const mesh & mesh, const size_type i); 

//! Aux function to read an array of string split by a delim. 
//! Store the results in a pre-constructed vector
std::vector<std::string> &
split(const std::string & s, 
	  char delim, 
	  std::vector<std::string> & elems
	  ) ;

//! Aux function to read an array of string split by a delim. 
//! Return a new vector
std::vector<std::string>
split(const std::string & s, 
	  char delim
	  ) ;

//! Build the integral of FE base functions
//! @f$ \Phi = \int_{\Omega} \phi(x)~dx @f$
/*!
 @param V    The array of computed values
 @param mim  The integration method to be used
 @param mf   The finite element method
 @param rg   The region where to integrate
 */ 
template
<typename VEC>
void 
asm_basis_function
	(VEC & V,
	 const mesh_im & mim,
	 const mesh_fem & mf,
	 const mesh_region & rg = mesh_region::all_convexes()
	 ) 		
{
	GMM_ASSERT1(mf.get_qdim() == 1, 
		"invalid data mesh fem for pressure (Qdim=1 required)");
	generic_assembly 
	assem("V$1(#1)+=comp(Base(#1));");
	assem.push_mi(mim);
	assem.push_mf(mf);
	assem.push_vec(V);
	assem.assembly(rg);
}

//! Aux function to extract the radius of the ith branch, R[i] 
template
<typename VEC>
scalar_type
compute_radius
	(const mesh_im & mim,
	 const mesh_fem & mf_coef,
	 const VEC & R,
	 const size_type & rg
	 ) 
{
	/*vector_type dof_enum;
	size_type fine=0;
	for(getfem::mr_visitor mrv(mf_coef.linked_mesh().region(rg)); !mrv.finished(); ++mrv){
		for(auto b: mf_coef.ind_basic_dof_of_element(mrv.cv())){
			dof_enum.emplace_back(b);
			fine++;
		}
	}
	size_type first_=dof_enum[0];
	return R[first_];*/

	//Luca
        getfem::generic_assembly assem;
        assem.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
        assem.push_mi(mim);
        assem.push_mf(mf_coef);
        assem.push_data(R);
        std::vector<scalar_type> v(1);
        assem.push_vec(v);
        assem.assembly(rg);

        getfem::generic_assembly assem1;
        assem1.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
        assem1.push_mi(mim);
        assem1.push_mf(mf_coef);
        std::vector<scalar_type> one(gmm::vect_size(R), 1.0);
        assem1.push_data(one);
        std::vector<scalar_type> measure(1);
        assem1.push_vec(measure);
        assem1.assembly(rg);

        //cout << "Integral radius: " << v[0] << endl << "Integral measure: " << measure[0] << endl << "Radius: " << v[0]/measure[0] << endl;
	//size_type rg_size = mf_coef.linked_mesh().region(rg).size();
	//cout << "Size: " << rg_size << endl << "h: " << measure[0]/rg_size << endl;
        return v[0]/measure[0];
}

template
<typename VEC>
scalar_type
old_compute_radius
	(const mesh_im & mim,
	 const mesh_fem & mf_coef,
	 const VEC & R,
	 const size_type & rg
	 ) 
{
	VEC rbasis(mf_coef.nb_dof());
	asm_basis_function(rbasis, mim, mf_coef, rg);
	size_type firstcv = mf_coef.linked_mesh().region(rg).index().first_true();
	size_type rg_size = mf_coef.linked_mesh().region(rg).size();
	gmm::scale(rbasis, 1.0/estimate_h(mf_coef.linked_mesh(), firstcv));
	return gmm::vect_sp(rbasis, R)/rg_size;
}
//! Compute the integral of the solution
template
<typename VEC>
scalar_type 
asm_mean(const mesh_fem &mf, const mesh_im &mim, const VEC &U) 
{
	getfem::generic_assembly assem;
	assem.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
	assem.push_mi(mim);
	assem.push_mf(mf);
	assem.push_data(U);
	std::vector<scalar_type> v(1);
	assem.push_vec(v);
	assem.assembly();

        getfem::generic_assembly assem1;
        assem1.set("u=data(#1); V()+=u(i).comp(Base(#1))(i)");
        assem1.push_mi(mim);
        assem1.push_mf(mf);
        std::vector<scalar_type> one(gmm::vect_size(U), 1.0);
        assem1.push_data(one);
        std::vector<scalar_type> measure(1);
        assem1.push_vec(measure);
        assem1.assembly();
	
        return v[0]/measure[0];
}
//! Aux function to save a gmm::col_vector as a plain list of values
template 
<typename L> 
void 
write_col_vector(std::ostream &o, const L &l) 
{
    for (size_type i = 0; i < gmm::mat_nrows(l); ++i) {

	if (gmm::mat_ncols(l) != 0) o << ' ' << l(i, 0);
	for (size_type j = 1; j < gmm::mat_ncols(l); ++j) o << ", " << l(i, j); 
      o << " \n";
    }
}

} /* end of namespace */

#endif
