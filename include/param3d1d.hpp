/* -*- c++ -*- (enableMbars emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
        Course on Advanced Programming for Scientific Computing
                      Politecnico di Milano
                          A.Y. 2014-2015
                  
                Copyright (C) 2015 Domenico Notaro
======================================================================*/
/*! 
	@file   param3d1d.hpp
	@author Domenico Notaro <domenico.not@gmail.com>
	@date   January 2016.
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
 
#ifndef M3D1D_PARAM3D1D_HPP_
#define M3D1D_PARAM3D1D_HPP_

#include <mesh1d.hpp>    // import_network_radius
#include <utilities.hpp> // compute_radius

namespace getfem {

//! Class to handle the physical parameter of the coupled 3D/1D model
struct param3d1d {

	// Dimensional physical parameters (microcirc applications)
	//! Average interstitial pressure [Pa]
	scalar_type P_; 
	//! Characteristic flow speed in the capillary bed [m/s]
	scalar_type U_; 
	//! Characteristic length of the problem [m]
	scalar_type d_; 
	//! Hydraulic conductivity of the interstitium [m^2]
	scalar_type k_; 
	//! Viscosity of the blood [kg/ms]
	scalar_type mu_v_; 	
	//! Viscosity of the interstitial fluid [kg/ms]
	scalar_type mu_t_; 
	//! Hydraulic conductivity of the capillary walls [m^2 s/kg]
	scalar_type Lp_;
        //! Hydraulic conductivity of the lymphatic vessels [m s/kg] (linear case)
        scalar_type Lp_LF_;
	//! Coefficient A of the lymphatic sigmoid [s-1]
        scalar_type A_LF_;
	//! Coefficient B of the lymphatic sigmoid [s-1]
        scalar_type B_LF_;
	//! Coefficient C of the lymphatic sigmoid [Pa]
        scalar_type C_LF_;
	//! Coefficient D of the lymphatic sigmoid [Pa]
        scalar_type D_LF_;
	// Dimensionless physical parameters (test-cases)
	//! Dimensionless average radius of the vessel network
	scalar_type Rav_;
	//! Dimensionaless radii of the vessel branches
	vector_type R_;
	//! Dimensionless conductivity of the tissue
	vector_type kt_;
	//! Dimensionless conductivity of the vessel wall
	vector_type Q_;
	//! Dimensionless conductivity of the vessel bed
	vector_type kv_;
        //! Dimensionless hydraulic conductivity of the lymphatic vessels (linear case)
        vector_type Q_LF_;
	// Dimensionless parameter A of the lymphatic sigmoid 
        scalar_type QLF_A_;
	// Dimensionless parameter B of the lymphatic sigmoid 
        scalar_type QLF_B_;
	// Dimensionless parameter C of the lymphatic sigmoid
        scalar_type QLF_C_;
	// Dimensionless parameter D of the lymphatic sigmoid
        scalar_type QLF_D_;
	// Dimensionless plasma oncotic pressure
	scalar_type pi_v_;
	// Dimensionless interstitial oncotic pressure
	scalar_type pi_t_;
	// Dimensionless reflection coefficient
	scalar_type sigma_;
	// Utils
	//! File .param
	ftool::md_param FILE_;
	//! Finite Element Method for tissue data
	getfem::mesh_fem mf_datat_;
	//! Finite Element Method for vessel data
	getfem::mesh_fem mf_datav_;
	// Methods
	//! Build the arrays of dimensionless parameters
	void build(ftool::md_param & fname, 
			const getfem::mesh_fem & mf_datat,
			const getfem::mesh_fem & mf_datav
			) 
	{
		FILE_ = fname;
		mf_datat_ = mf_datat;
		mf_datav_ = mf_datav;
		size_type dof_datat = mf_datat_.nb_dof();
		size_type dof_datav = mf_datav_.nb_dof();
		 
		bool IMPORT_RADIUS = FILE_.int_value("IMPORT_RADIUS");
		bool NONDIM_PARAM  = FILE_.int_value("TEST_PARAM");
		bool EXPORT_PARAM  = FILE_.int_value("EXPORT_PARAM");
		bool LINEAR_LYMPHATIC_DRAINAGE = FILE_.int_value("LINEAR_LYMPHATIC_DRAINAGE");

		// Check
		if (IMPORT_RADIUS)
			GMM_ASSERT1(NONDIM_PARAM == 0,
				"try to import non constant (dimensionless) radius:" 
				"please insert dimensional parameters");
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless radius R'... "   << endl;
		#endif
		if (!IMPORT_RADIUS) { 	/* case R' = const */

			if (NONDIM_PARAM) // we assume it is already non-dimensional
				Rav_ = FILE_.real_value("RADIUS", "Vessel average radius");
			else // to be non-dimensionalized
				Rav_ = FILE_.real_value("RADIUS", "Vessel average radius")/FILE_.real_value("d");
			R_.assign(dof_datav, Rav_);

		} else { 				/* case R' = R'(s) */
			std::string RFILE = FILE_.string_value("RFILE"); 
			cout << "  Importing radius values from file " << RFILE << " ..." << endl;
			std::ifstream ist(RFILE);
			if (!ist) cerr << "impossible to read from file " << RFILE << endl;
			import_network_radius(R_, ist, mf_datav_);
		}
		#ifdef M3D1D_VERBOSE_
		cout << "  Assembling dimensionless permeabilities kt, Q, kv ... "   << endl;
		#endif
		if (NONDIM_PARAM) {
			// Import dimensionless params from FILE_
			scalar_type ktval = FILE_.real_value("Kt"); 
			scalar_type Qval  = FILE_.real_value("Q"); 
			scalar_type kvval = FILE_.real_value("Kv");
			pi_t_ = FILE_.real_value("pi_t_adim");
			pi_v_ = FILE_.real_value("pi_v_adim");
			sigma_ = FILE_.real_value("sigma");
			if(!LINEAR_LYMPHATIC_DRAINAGE)
			{
                        QLF_A_ = FILE_.real_value("QLF_A", "Dimensionless parameter A of lymphatic drainage");
                        QLF_B_ = FILE_.real_value("QLF_B", "Dimensionless parameter B of lymphatic drainage");
                        QLF_C_ = FILE_.real_value("QLF_C", "Dimensionless parameter C of lymphatic drainage");
                        QLF_D_ = FILE_.real_value("QLF_D", "Dimensionless parameter D of lymphatic drainage");
 			Q_LF_.assign(dof_datat, 0);
			}
			else
			{
                        scalar_type QLFval  = FILE_.real_value("Q_LF");
                        Q_LF_.assign(dof_datat, QLFval);
			}
			// Fill the data arrays
			  kt_.assign(dof_datat,  ktval);
			  kv_.assign(dof_datav,  kvval);
			   Q_.assign(dof_datav,   Qval);

		} 
		else {
			// Import dimensional params from FILE_
			P_  = FILE_.real_value("P", "average interstitial pressure [Pa]"); 
			U_  = FILE_.real_value("U", "characteristic flow speed in the capillary bed [m/s]"); 
			d_  = FILE_.real_value("d", "characteristic length of the problem [m]"); 
			k_  = FILE_.real_value("k", "permeability of the interstitium [m^2]"); 
			mu_v_ = FILE_.real_value("mu_v", "blood viscosity [kg/ms]"); 
			mu_t_ = FILE_.real_value("mu_t", "interstitial fluid viscosity [kg/ms]"); 
			pi_t_ = FILE_.real_value("Pi_t", "interstitial oncotic pressure [Pa]"); 
			pi_v_ = FILE_.real_value("Pi_v", "fluid oncotic pressure [Pa]"); 
			sigma_ = FILE_.real_value("sigma", "reflection coefficient [-]"); 
			Lp_   = FILE_.real_value("Lp", "permeability of the vessel walls [m^2 s/kg]"); 
			if(LINEAR_LYMPHATIC_DRAINAGE)
			{
                        Lp_LF_ = FILE_.real_value("Lp_LF", "permeability of lymphatic vessels [(m s) / kg] ");
                        Q_LF_.assign(dof_datat, Lp_LF_*P_*d_/U_);
			}
			// Compute the dimenless params
			kt_.assign(dof_datat, k_/mu_t_*P_/U_/d_);
			pi_t_ = pi_t_/P_;
			pi_v_ = pi_v_/P_; 

			for (auto r : R_){ // C++11-only!
				kv_.emplace_back(pi/8.0/mu_v_*P_*d_/U_*r*r*r*r);
				Q_.emplace_back(2*pi*Lp_*P_/U_*r);
			}

			// Fixed Point Method for Lymphatic System (Sigmoid)
			if(!LINEAR_LYMPHATIC_DRAINAGE)
			{
                        A_LF_ = FILE_.real_value("A_LF", "First Coefficient (A) of the lymphatic flow [s-1] ");
                        B_LF_ = FILE_.real_value("B_LF", "Second Coefficient (B) of the lymphatic flow [s-1]");
                        C_LF_ = FILE_.real_value("C_LF", "Third Coefficient (C) of the lymphatic flow [Pa]");
                        D_LF_ = FILE_.real_value("D_LF", "Fourth Coefficient (D) of the lymphatic flow [Pa]");
			QLF_A_ = A_LF_/U_*d_;	
			QLF_B_ = B_LF_/U_*d_;
			QLF_C_ = C_LF_/P_;
			QLF_D_ = D_LF_/P_;
 			Q_LF_.assign(dof_datat, 0);
			}
		}
		// Check values
		GMM_ASSERT1(kt_[0] != 0, "wrong tissue conductivity (kt>0 required)"); 
		GMM_ASSERT1(kv_[0] != 0, "wrong vessel bed conductivity (kv>0 required)");
		if (Q_[0] == 0) cout << "Warning: uncoupled problem (Q=0)" << endl;
		
		if (EXPORT_PARAM){
			std::string ODIR = FILE_.string_value("OutputDir","OutputDirectory");
			getfem::vtk_export exp(ODIR+"radius.vtk");
			exp.exporting(mf_datav_);
			exp.write_mesh();
			exp.write_point_data(mf_datav_, R_, "R");
			getfem::vtk_export expQ(ODIR+"conductivity.vtk");
			expQ.exporting(mf_datav_);
			expQ.write_mesh();
			expQ.write_point_data(mf_datav_, Q_, "Q");
		}

	}


	//! Get the radius at a given dof
	inline scalar_type R  (size_type i) { return R_[i];  } const
	//! Get the tissue permeability at a given dof
	inline scalar_type kt (size_type i) { return kt_[i]; } const
	//! Get the vessel bed permeability at a given dof
	inline scalar_type kv (size_type i) { return kv_[i]; } const
	//! Get the vessel wall permeability at a given dof
	inline scalar_type Q  (size_type i) { return Q_[i];  } const
	//! Get the radius at a given mesh_region
	scalar_type R  (const getfem::mesh_im & mim, const size_type rg) { 
		return compute_radius(mim, mf_datav_, R_, rg);  
	}
	//! Get the vessel bed permeability at a given mesh_region
	scalar_type kv  (const getfem::mesh_im & mim, const size_type rg) { 
		return compute_radius(mim, mf_datav_, kv_, rg);  
	}
	//! Get the radius
	vector_type & R (void) { return R_; }
	//! Get the vessel wall permeabilities
	vector_type & Q (void) { return Q_; }
	//! Get the vessel bed permeabilities
	vector_type & kv (void) { return kv_; }
        //! Get the lymphatic vessels permeability
        inline scalar_type Q_LF (size_type i) { return Q_LF_[i]; } const
	//! Get the coefficient of lymphatic vessel
	inline scalar_type QLF_a  (void) { return QLF_A_; } const
	inline scalar_type QLF_b  (void) { return QLF_B_;  } const
	inline scalar_type QLF_c  (void) { return QLF_C_;  } const
	inline scalar_type QLF_d  (void) { return QLF_D_;  } const
	//! Get the interstitial oncotic pressure
	inline scalar_type pi_t  (void) { return pi_t_; } const
	//! Get the plasma oncotic pressure
	inline scalar_type pi_v  (void) { return pi_v_; } const
	//! Get the reflection coefficient
	inline scalar_type sigma  (void) { return sigma_; } const
	//! Overloading of the output operator
	friend std::ostream & operator << (
		std::ostream & out, const param3d1d & param
		)
	{ 
		out << "--- PHYSICAL PARAMS ------" << endl;
		out << "  R'     : "                << param.R_[0] << endl; 
		out << "  kappat : "                << param.kt_[0] << endl; 
		out << "  Q      : "                << param.Q_[0] << endl; 
                out << "  kappav : "                << param.kv_[0] << endl;
		out << "--------------------------" << endl;

		return out;            
	}
}; /* end of class */

} /* end of namespace */

#endif
