#ifndef DARCYPRECOND
#define DARCYPRECOND

#include "gmm_fix.hpp"

// Useful type definitions (built using the predefined types in Gmm++)
//! Special class for small (dim < 16) vectors 
using bgeot::base_small_vector;  
//! Geometrical nodes (derived from base_small_vector)
using bgeot::base_node;   	 	 
//! Double-precision FP numbers 
using bgeot::scalar_type; 	 	
//! Unsigned long integers 
using bgeot::size_type;   	 	 
//! Short integers 
using bgeot::short_type;         
//! Type for vector of integers 
typedef std::vector<size_type> vector_size_type;
//! Type for dense vector
typedef std::vector<scalar_type> vector_type;
//! Type for sparse vector (std::vector) 
typedef gmm::rsvector<scalar_type> sparse_vector_type;
//! Type for dense matrix
typedef gmm::row_matrix<vector_type> matrix_type;
//! Type for sparse matrix
typedef gmm::row_matrix<sparse_vector_type> sparse_matrix_type;   

#include <vector>
#include <getfem/getfem_generic_assembly.h>
#include <gmm/gmm_precond_diagonal.h>
#include <gmm/gmm_superlu_interface.h>

//#define USE_SAMG 1

#ifdef USE_SAMG
#include "AMG_Interface.hpp"
#endif
enum method { DIRECT=2,ITER=1 };


template <class MATRIX>
class darcy_precond
{
public:
    // TODO boundary conditions (Dirichlet, Robin) and coefficient (kappa, beta)
    darcy_precond(const MATRIX &A,
                  const getfem::mesh_fem mf_p,
                  const getfem::mesh_im mim);

    getfem::size_type nrows() const {
        return gmm::mat_nrows(A_) + gmm::mat_nrows(S_);
    }

    getfem::size_type ncols() const {
        return gmm::mat_ncols(A_) + gmm::mat_ncols(S_);
    }

    template <class L2, class L3>
    void mult(const L2 &src, L3 &dst) const
    {
        const getfem::size_type n1 = gmm::mat_ncols(A_),
                                n2 = gmm::mat_ncols(S_);

        gmm::mult(pA_, gmm::sub_vector(src, gmm::sub_interval(0, n1)),
                  gmm::sub_vector(dst, gmm::sub_interval(0, n1)));
#ifdef USE_SAMG
        std::vector<double> x(n2),b(n2);
        gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2)),b);
        gmm::clear(x);
        amg_.solve(S_, x , b , 1); 
        gmm::copy(amg_.getsol(),gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));
#else
        slu_.solve(gmm::sub_vector(dst, gmm::sub_interval(n1, n2)),
                   gmm::sub_vector(src, gmm::sub_interval(n1, n2)));
#endif
    }

private:
    const MATRIX &A_;
    gmm::diagonal_precond<MATRIX> pA_;
    MATRIX S_;

#ifdef USE_SAMG
    mutable AMG amg_;

#else
    gmm::SuperLU_factor<double> slu_;
#endif
};


namespace gmm {
    template <class MATRIX>
    struct linalg_traits<::darcy_precond<MATRIX>> {
        using this_type = ::darcy_precond<MATRIX>;
        using sub_orientation = owned_implementation;

        static size_type nrows(const this_type &m) { return m.nrows(); }
        static size_type ncols(const this_type &m) { return m.ncols(); }
    };
} // namespace gmm


template <class MATRIX>
darcy_precond<MATRIX>::darcy_precond(const MATRIX &A,
                                     const getfem::mesh_fem mf_p,
                                     const getfem::mesh_im mim)
: A_(A)
, pA_(A)
#ifdef USE_SAMG
, amg_("Schur")
#endif
{
    const getfem::size_type nb_dof_p = mf_p.nb_dof();
    const getfem::mesh &mesh = mf_p.linked_mesh();
    getfem::mesh_region inner_faces = getfem::inner_faces_of_mesh(mesh);
    getfem::mesh_region outer_faces;
    getfem::outer_faces_of_mesh(mesh, outer_faces);
    getfem::ga_workspace wp;

    std::vector<double> p(nb_dof_p);
    wp.add_fem_variable("p", mf_p, gmm::sub_interval(0, nb_dof_p), p);

    wp.add_expression("Grad_p.Grad_Test_p", mim);
    wp.add_expression("-0.5 * (Grad_p + Interpolate(Grad_p, neighbour_elt)).Normal"
                        " * (Test_p - Interpolate(Test_p, neighbour_elt))"
                      "-0.5 * (Grad_Test_p + Interpolate(Grad_Test_p, neighbour_elt)).Normal"
                        " * (p - Interpolate(p, neighbour_elt))"
                      "+2 / element_size * (p - Interpolate(p, neighbour_elt))"
                        " * (Test_p - Interpolate(Test_p, neighbour_elt))",
                      mim, inner_faces);
 // in case of dir condition
//     wp.add_expression("1/element_size*p*Test_p",
//                       mim, outer_faces);

    wp.add_expression("0.01*p*Test_p",  mim, outer_faces);
    
           // for dir=0 bc-condition
    wp.assembly(2);

    gmm::copy(wp.assembled_matrix(), S_);

#ifdef USE_SAMG
  amg_.convert_matrix(S_);
#else
    slu_.build_with(S_);
#endif
}

#endif // ifndef darcyprecond
