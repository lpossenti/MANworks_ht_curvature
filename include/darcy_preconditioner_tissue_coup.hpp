#ifndef DARCYPRECONDT
#define DARCYPRECONDT

#include "gmm_fix.hpp"
#include <vector>
#include <getfem/getfem_generic_assembly.h>
#include <gmm/gmm_precond_diagonal.h>
#include <gmm/gmm_superlu_interface.h>

//#define USE_SAMG 1
//#define USE_MP

#ifdef USE_SAMG
#include "AMG_Interface.hpp"
#endif

template <class MATRIX>
class darcy_precond_tissue_coup
{
public:
    // TODO boundary conditions (Dirichlet, Robin) and coefficient (kappa, beta)
    darcy_precond_tissue_coup(const MATRIX &A,
	          const MATRIX &Btt,
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


// Preconditionig tissue  first block
        gmm::mult(pA_, gmm::sub_vector(src, gmm::sub_interval(0, n1)),
                  gmm::sub_vector(dst, gmm::sub_interval(0, n1)));


// Preconditionig tissue second block


#ifdef USE_SAMG
       
         std::vector<double> x(n2),b(n2);
        gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2)),b);
    
        gmm::clear(x);
        amg_.solve(S_, x , b , 1);
        gmm::copy(amg_.getsol(),gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));
#else
        // one step
           slu_.solve(gmm::sub_vector(dst, gmm::sub_interval(n1, n2)),
                      gmm::sub_vector(src, gmm::sub_interval(n1, n2)));

#endif
        
         // one step
         // std::vector<double> res_temp(n2);
          //gmm::copy(gmm::sub_vector(src, gmm::sub_interval(n1, n2)),res_temp);
         // slu_ct_.solve(res_temp,gmm::sub_vector(src, gmm::sub_interval(n1, n2)));
         
         // double step
        /* std::vector<double> res_temp(n2);
         gmm::copy(res_p,res_temp);
         slu_ct_.solve(res_temp,res_p);
         gmm:add(res_temp,gmm::sub_vector(dst, gmm::sub_interval(n1, n2)));*/


    }

private:
    const MATRIX &A_;
    gmm::diagonal_precond<MATRIX> pA_;
    mutable MATRIX S_;
    const MATRIX &Ct_;
#ifdef USE_SAMG
    mutable AMG amg_;

#else
    gmm::SuperLU_factor<double> slu_;
#endif

};


namespace gmm {
    template <class MATRIX>
    struct linalg_traits<::darcy_precond_tissue_coup<MATRIX>> {
        using this_type = ::darcy_precond_tissue_coup<MATRIX>;
        using sub_orientation = owned_implementation;

        static size_type nrows(const this_type &m) { return m.nrows(); }
        static size_type ncols(const this_type &m) { return m.ncols(); }
    };
} // namespace gmm


template <class MATRIX>
darcy_precond_tissue_coup<MATRIX>::darcy_precond_tissue_coup(const MATRIX &A,
		                     const MATRIX &Btt,
                                     const getfem::mesh_fem mf_p,
                                     const getfem::mesh_im mim)
: A_(A)
, pA_(A) ,Ct_(Btt)
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

    #ifdef USE_MP
    wp.add_expression("p*Test_p",mim);
    #else
    wp.add_expression("Grad_p.Grad_Test_p", mim);
    wp.add_expression("-0.5 * (Grad_p + Interpolate(Grad_p, neighbour_elt)).Normal"
                        " * (Test_p - Interpolate(Test_p, neighbour_elt))"
                      "-0.5 * (Grad_Test_p + Interpolate(Grad_Test_p, neighbour_elt)).Normal"
                        " * (p - Interpolate(p, neighbour_elt))"
                      "+2 / element_size * (p - Interpolate(p, neighbour_elt))"
                        " * (Test_p - Interpolate(Test_p, neighbour_elt))",
                      mim, inner_faces);
    // to remove in case of mix/neumann condition
    //wp.add_expression("1/element_size*p*Test_p", mim, outer_faces); // for DIR
    wp.add_expression("0.01*p*Test_p",  mim, outer_faces);
    
     #endif
    wp.assembly(2);

     gmm::copy(wp.assembled_matrix(), S_);
    // gmm::row_matrix< std::vector<double> > M1(nb_dof_p, nb_dof_p);
    // gmm::copy(wp.assembled_matrix(), M1);
    //gmm::copy(M1, S_);   (S+C_t)
    // gmm::add(Btt,M1);
     
#ifdef USE_SAMG
    amg_.convert_matrix(S_);
#else
    slu_.build_with(S_);
#endif
}


#endif // ifndef darcyprecond
