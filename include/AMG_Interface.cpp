#ifdef WITH_SAMG
#include "AMG_Interface.hpp" 
#include <memory>

AMG::AMG(std::string name)
{
	std::cout<<"Build class AMG for "<< name << std::endl;
		
}


AMG::~AMG(void) {
	// delete [] 
	std::cout << "AMG deleted" << std::endl;
		
}

void AMG::solve(gmm::csr_matrix<scalar_type> A_csr, std::vector<scalar_type> U, std::vector<scalar_type> B, int solver_type)
{
	APPL_INT npnt,nsys,matrix;
	matrix=22; nsys=1;npnt= 0;   //_q_dof + _l_dof;
	APPL_INT ndiu = 1;           // dimension of (dummy) vector iu
	APPL_INT ncyc;
	APPL_INT ndip = npnt;        // dimension of (dummy) vector ip
	float told,tnew,tamg;
	//SAMG configuration
		
	APPL_INT * iu;
	if (nsys==1) iu= new APPL_INT[1];
	else {
		iu  = new APPL_INT[nnu_];
		ndiu   = nnu_;
		for(int iiu=0;iiu<nnu_;iiu++){
			//only working for nsys = 2 check it out for larger systems				 
			if(iiu < _Ut ) iu[iiu]=1;
			else if (iiu < _Pt + _Ut ) iu[iiu]=2;
			//else if (iiu < _Pt + _Ut+_Uv ) iu[iiu]=3;
			//else if (iiu < _Pt + _Ut + _Pv + _Uv ) iu[iiu]=4;
				}
		}// endl else
		//=====================================================
		/// printing format
		char ch[] = "f";int len=1;
		SAMG_SET_IOFORM(&ch[0],&len);
		
		
		//==================================================			
		// int nrc	;
		// int nptmax; // max number of points for which we switch to nrc_emergency
		// int clsolver_finest;	// 
		// SAMG_GET_NRC(&nrc); // retreive nrc=solver for coarser levels page 67 user guide
		// std::cout<<"..value of nrc======"<<nrc<<std::endl;
		// ================================
		// 1  Iterative application of currently selected smoother 
		// 2  ILU(0) preconditioned CG 
		// 3  ILUT preconditioned BiCGstab 
		// 4  Diagonally preconditioned CG 
		// 5  Diagonally preconditioned BiCGstab 
		// 6  Dense direct solver (standard)  
		// 7  Sparse direct solver (standard)  
		// 8  Least squares solver (LSQ; robust but very expensive!) 
		// 9  Dummy interface routine (→ Section 8.3.1)  
		// 10  Dense direct solver (highly tuned, → cldense_ctrl below) 
		// 11  Sparse direct solver (Intel’s pardiso) if MKL has been linked 
		// 99  Another instance of SAMG 
		// ================================



		//==================================================	



		if (solver_type==2){			
			int levelx;	SAMG_GET_LEVELX(&levelx); // retreive levelx=number of maximum coarsening levels 
			std::cout<<"..value of levelx======"<<levelx<<std::endl;
			levelx=3;SAMG_SET_LEVELX(&levelx);// change levelx=number of maximum coarsening levels 
			SAMG_GET_LEVELX(&levelx);// retreive levelx=number of maximum coarsening levels 
			std::cout<<"..check if change value of levelx======"<<levelx<<std::endl;
			int  nrc=9;int  nrc_emergency=9;int nptmax=5000;int clsolver_finest=1;
			SAMG_SET_NRC(&nrc);// change  nrc=solver for coarser levels page 67 user guide
			SAMG_SET_NRC_EMERGENCY(&nrc_emergency);
			SAMG_SET_NPTMAX(&nptmax);// change  nptmax
			SAMG_SET_CLSOLVER_FINEST(&clsolver_finest);// change  nptmax
			// SAMG_GET_NRC(&nrc);// retreive  nrc=solver for coarser levels page 67 user guide
			// std::cout<<"..check if change value of nrc======"<<nrc<<std::endl;
			ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
			} //direct_solver
	
#ifdef AMG_STAND_ALONE
		//Both ncgrad (the "nd number in ncyc) and ncgrad_default must be equal to 0 to use SAMG solver as a stand-alone solver (not as a preconditioner)

	int ncgrad_default=0;
	SAMG_SET_NCGRAD_DEFAULT(&ncgrad_default);
	ncyc      = 10050;    // V-cycle as pre-conditioner for CG; at most 50 iterations
#endif //amg_stand_alone

#ifdef AMG_ACCELERATED
	ncyc      = 11050;    // V-cycle as pre-conditioner for CG; at most 50 iterations	
#endif //amg_accelerated
	
	if (solver_type==1){
		//int levelx=2;SAMG_SET_LEVELX(&levelx);
		// // change levelx=number of maximum coarsening levels 
		int  nrc=3;int  nrc_emergency=3;
		SAMG_SET_NRC(&nrc);// change  nrc=solver for coarser levels page 67 user guide
		SAMG_SET_NRC_EMERGENCY(&nrc_emergency);
		ncyc      = 11050;   
		// int nrd=931;  SAMG_SET_NRD(&nrd);
		// int nru=931;  SAMG_SET_NRU(&nru);
		// int it_pre=2;SAMG_SET_ITER_PRE(&it_pre);
		}
	
	// APPL_INT * iscale = new APPL_INT[1];
	// this vector (iscale) indicates which uknowns require scaling if 0 no scaling
	APPL_INT * iscale = new APPL_INT[nsys]; for(int i_sys=0; i_sys<nsys; i_sys++) iscale[i_sys]=0; 
	//iscale[2]=1; 
	APPL_INT * ip     = new APPL_INT[1];
		
	// APPL_INT * ip     = new APPL_INT[npnt];
	// for (int ipt=0; ipt< npnt;  ipt++) {ip[ipt] =  _pt2uk[ipt];}
	// nsolve = napproach nxtyp internal 
		
	APPL_INT nsolve    = 2;        // results in scalar approach (current system is scalar)
	APPL_INT ifirst    = 1;        // first approximation = zero
	double  eps       = 1.0e-12;   // required (relative) residual reduction
	
	//ncyc      = 50050;
	APPL_INT n_default = 20;       // select default settings for secondary parameters

	// CURRENTLY AVAILABLE: 10-13, 15-18, 20-23, 25-28
	// NOTE: the higher the respective SECOND digit, the
	// more aggressive the coarsening (--> lower memory at
	// the expense of slower convergence)
	APPL_INT iswtch    = 4100+n_default; // complete SAMG run ....
	if (first_){
		first_=false;
		}
	else{
		iswtch    = 1100+n_default; // complete SAMG run ....
	}
	// ... memory de-allocation upon return ....
	// ... memory extension feature activated ....
	// ... residuals measured in the L2-norm .....
	// ... secondary parameter default setting # n_default
	
	// ===> Secondary parameters which have to be set if n_default=0
	//      (at the same time a demonstration of how to access secondary or hidden parameters)
		
	double a_cmplx   = 2.2;      // estimated dimensioning
	double g_cmplx   = 1.7;      // estimated dimensioning
	double w_avrge   = 2.4;      // estimated dimensioning
	double p_cmplx   = 0.0;      // estimated dimensioning (irrelevant for scalar case)
	double  chktol    = -1.0;    // input checking de-activated (we know it's ok!)
	// int  iswit=5;	SAMG_SET_ISWIT(&iswit);// iswit reuse data User at 93
	//============================================
	//     idump controls the matrix dumping of SAMG				
	// 1  Standard print output, no matrix dump. 
	// 2‐6  Write matrices to disk: level 2 up to level idmp. 
	// 7  Write matrices to disk: level 2 up to the coarsest level. 
	// 8  Write finest‐level matrix to disk (incl. right hand side etc.). 
	// 9  Write all matrices to disk. 
	APPL_INT idump     = 9;       // minimum output during setup
	//============================================
	// iout page 44 Userguide. it controls display outpu. default 2 very verbose 43
	APPL_INT iout      = -1;        // display residuals per iteration and work statistics
	// Coarsening type
	// int  ncg=54; //int  nrc_emergency=4;
	// SAMG_SET_NCG(&ncg);// change  nrc=solver for coarser levels page 67 user guide
	// output:
	APPL_INT ierr,ierrl,ncyc_done;
	double res_out,res_in;
	//end
        //video output
	// har logfile[] = "s"; int len_log = 1; 
	//	SAMG_SET_LOGFILE(&logfile,&len_log);// change levelx=number of maximum coarsening levels  
	// for(int xx=0;xx<200;xx++)std::cout<<a_samg[xx]<<" "<<a[xx]<< " -- "<<ja_samg[xx]<<" "<<ja[xx] <<std::endl;
	SAMG_CTIME(&told);
	double time_SAMG=gmm::uclock_sec();
	std::cout<<"start solving with samg " << std::endl;
//=======================================================================================
	SAMG(&nnu_,&nna_,&nsys,&ia_samg_[0],&ja_samg_[0],
			&a_samg_[0],
			&B.front(), &U.front(),
			//&b_samg[0], &u_samg[0],
			&iu[0],&ndiu,
			&ip[0],&ndip,&matrix,&iscale[0],
			&res_in,&res_out,&ncyc_done,&ierr,
			&nsolve,&ifirst,&eps,&ncyc,&iswtch,
			&a_cmplx,&g_cmplx,&p_cmplx,&w_avrge,
			&chktol,&idump,&iout);
//=======================================================================================
	std::cout << std::endl<<"*** time to solve the system using SAMG with gmm time: " << gmm::uclock_sec() - time_SAMG << " seconds\n";
	gmm::resize(sol_vec,nnu_);
	for(int i = 0 ; i < nnu_ ; i++ ){sol_vec[i]=U[i];}
	return;
}
//===============================================================
//===============================================================
//===============================================================
//===============================================================


void AMG::set_pt2uk(int * dofpt , int q_dof, int l_dof, int npts){
	
		std::cout << std::string(100, '=') << std::endl;
		std::cout<< "AMG::set_pt2uk start"  << std::endl;
		std::cout << std::string(100, '=') << std::endl;
		
		
		_q_dof=q_dof; _l_dof=l_dof; _npts=npts;
		_pt2uk.resize(_q_dof + l_dof);
		for (int idof=0; idof< q_dof+l_dof; idof++)  _pt2uk[idof] = dofpt[idof] ;
			
				
				std::cout << std::string(100, '=') << std::endl;
				std::cout<< "AMG::set_pt2uk end "  << std::endl;
				std::cout << std::string(100, '=') << std::endl;
}

void AMG::set_dof(int pt , int ut, int pv, int uv){
	_Pt=pt;_Ut=ut;_Pv=pv; _Uv=uv;
}



//==========================================================================
// Convert and store a getfem csr matrix into the AMG itnerface
//==========================================================================
void AMG::convert_matrix(gmm::csr_matrix<scalar_type> A_csr)
{
	std::cout<<" AMG::convert_matrix::Start building the matrix"<<std::endl; 
	std::cout<<"*** parameters SAMG matrix   "<<A_csr.nrows()<<std::endl;	
	std::cout<<"*** parameters SAMG matrix   "<<A_csr.ncols()<<std::endl;
	nnu_=A_csr.nc; nna_=A_csr.pr.size();
	a_samg_=new double[nna_];
	ja_samg_=new APPL_INT[nna_];
	ia_samg_=new APPL_INT[nnu_+1];
	ia_samg_[0]=1;
	unsigned int offset=0;
	for(int ia=0;ia<nnu_+1-1;ia++){
		unsigned int nonzero=A_csr.jc.at(ia+1)-A_csr.jc.at(ia);
		ia_samg_[ia+1] = nonzero + ia_samg_[ia] ;
		bool diag=true;
		for (int innz=offset; innz<offset+nonzero; innz++ ){
			if( ia == (int) A_csr.ir.at(innz))
					{	//WARNING: not working for a matrix with zero on the diagonal
						//	std::cout<<"diagonal term"<<std::endl; 
						if(fabs(A_csr.pr.at(innz))<1E-30)std::cout<<"*****diagonal zero"<<std::endl;
						a_samg_[offset]=A_csr.pr.at(innz);
						ja_samg_[offset]=A_csr.ir.at(innz)+1;
						diag=false;			
					}
				}
			if (diag) {
				std::cout<<"dsufhasdjkasdfasdf jhasdfasdfjkasdffasdfjkl"<<std::endl;
				ia_samg_[ia+1] +=1;
				a_samg_[offset]=1.0;
				ja_samg_[offset]=ia_samg_[ia];
			}
			int shift=1;
			for (int innz=offset; innz<offset+nonzero; innz++ )
				if( ia != (int) A_csr.ir.at(innz))
					{
						//	std::cout<<"non diagonal term "<< ia << " " << A_csr.ir.at(innz) <<std::endl; 
						if(fabs(A_csr.pr.at(innz))<1E-30) std::cout<<"**non diagonal zero "<<A_csr.pr.at(innz)<<std::endl;
						a_samg_[offset+shift]=A_csr.pr.at(innz);
						ja_samg_[offset+shift]=A_csr.ir.at(innz)+1;
						shift++;
					}
			//std::cout<<std::endl;
			// std::cout<<"non zero in i "<<ia<<" are "<< nonzero<<std::endl; 
			offset+=nonzero;
		}

	std::cout<<" AMG::convert_matrix:: end matrix conversion "<< std::endl;
	return;
}

#endif // if def with samg









