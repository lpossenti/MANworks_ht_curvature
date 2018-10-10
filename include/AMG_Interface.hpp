// AMG_Interrface.h
#ifdef WITH_SAMG
#ifndef aamg
#define aamg


#include<iostream> 

#include "getfem/getfem_mesher.h"
#include "gmm/gmm.h"

#include "samg.h"
/* default 4 Byte integer types */
#ifndef APPL_INT
#define APPL_INT int
#endif

using bgeot::scalar_type; 

class AMG {
  private:
std::vector<scalar_type> sol_vec; /// solution vector
std::vector<int> _pt2uk; /// point to oknow vector
int _q_dof;              /// quadratic dof
int _l_dof;				/// linear dof
int _Pt, _Ut, _Pv, _Uv;
int _npts;
bool first_=true;
double * a_samg_;
APPL_INT *ja_samg_;
APPL_INT *ia_samg_;
APPL_INT nnu_,nna_;
  public:
  // ======== costructor the class ========================
  AMG(std::string name);
  void set_dof(int , int, int, int);
    // ======== destructor the class ========================
  ~AMG();  // This is the destructor: declaration
  // ======== generation af matrix
  void convert_matrix(gmm::csr_matrix<scalar_type> A_csr);
    // ======== solver of the class ========================
  void solve(gmm::csr_matrix<scalar_type> A_csr, std::vector<scalar_type> U, std::vector<scalar_type> B, int solver_type );
  // =========== set to point to uknown vector ========================
  void set_pt2uk(int * dofpt , int q_dof, int l_dof, int npts);
  
    // =========== return the solution ========================
  std::vector<scalar_type> getsol(){return sol_vec;}
};

#endif
#endif // WITH_SAMG
