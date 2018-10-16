/* -*- c++ -*- (enables emacs c++ mode) */
/*==============================================================================
          "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                            Politecnico di Milano
                                A.Y. 2016-2017
                  
                Copyright (C) 2017 Luca Possenti - Simone Di Gregorio
================================================================================*/
/*! 
  @file   problemHT.cpp
  @author Luca Possenti <luca.possenti@polimi.it>
  @author Simone Di Gregorio <simone.digre@gmail.com>
  @date   March 2017.
  @brief  Definition of the main class for the hematocrit transport problem.
 */

#include <problemHT.hpp>
#include <cmath>



namespace getfem {

/////////// Initialize the problem ///////////////////////////////////// 
void 
problemHT::init(int argc, char *argv[])
{
	/*//1. Read the .param filename from standard input
	PARAM.read_command_line(argc, argv);
	//2. Import data (algorithm specifications, boundary conditions, ...)
	import_data();*/ // DONE IN problemHT::HEMATOCRIT_TRANSPORT
	//3. Import mesh vessel network (1D)
	build_mesh();
	//4. Set finite elements and integration methods
	set_im_and_fem();
	//5. Build problem parameters
	build_param();
	//6. Build the list of vessel boundary (and junction) data
	build_vessel_boundary();
}
bool
problemHT::HEMATOCRIT_TRANSPORT(int argc, char *argv[])
{	//1. Read the .param filename from standard input
	PARAM.read_command_line(argc, argv);
	//2. Import data (algorithm specifications, boundary conditions, ...)
	import_data();
	return descrHT.HEMATOCRIT_TRANS;
}
void
problemHT::import_data(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing descriptors for hematocrit problems ..." << endl;
	#endif
	descrHT.import(PARAM);
	#ifdef M3D1D_VERBOSE_
	cout << descrHT;
	#endif
}


void
problemHT::build_mesh(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Importing the 1D mesh for the vessel in hematocrit problem... "   << endl;
	#endif
	std::ifstream ifs(descrHT.MESH_FILEH);
	GMM_ASSERT1(ifs.good(), "impossible to read from file " << descrHT.MESH_FILEH);
 	vector_type Uv( dof.Uv()); gmm::clear(Uv);
	gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())),  Uv);
	import_pts_file_HT(ifs, meshh, BCv_HT, nb_vertices, Uv, descr.MESH_TYPEV, mimv, mf_Uvi);
	nb_branches = nb_vertices.size();
	ifs.close();
}

void
problemHT::set_im_and_fem(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEMs for hematocrit problems ..." << endl;
	#endif
	bgeot::pgeometric_trans pgt_h = bgeot::geometric_trans_descriptor(descrHT.MESH_TYPEH);
	pfem pf_H = fem_descriptor(descrHT.FEM_TYPEH);
	pfem pf_coefh = fem_descriptor(descrHT.FEM_TYPEH_DATA);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting IMs and FEMs for vessel branches (hematocrit)..." << endl;
	#endif
	mf_Hi.reserve(nb_branches);
	mf_coefhi.reserve(nb_branches);
	for(size_type i=0; i<nb_branches; ++i){

		mesh_fem mf_tmp(meshv);
		mf_tmp.set_finite_element(meshv.region(i).index(), pf_coefh);
		mf_coefhi.emplace_back(mf_tmp);
		mf_tmp.clear();

		mf_tmp.set_finite_element(meshv.region(i).index(), pf_H);
		mf_Hi.emplace_back(mf_tmp);
		mf_tmp.clear();
	}
	mf_coefh.set_finite_element(meshv.convex_index(), pf_coefh);

	#ifdef M3D1D_VERBOSE_
	cout << "Setting FEM dimensions for hematocrit problems ..." << endl;
	#endif
	dofHT.set(mf_Hi, mf_coefv);

	#ifdef M3D1D_VERBOSE_
	cout << std::scientific << dofHT;
	#endif
}

void
problemHT::build_param(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building parameters for hematocrit problems ..." << endl;
	#endif
	paramHT.build(PARAM, mf_coefv);
	#ifdef M3D1D_VERBOSE_
	cout << paramHT;
	#endif
}

void 
problemHT::build_vessel_boundary(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Building hematocrit boundary ..." << endl;
	#endif
try {

	dal::bit_vector junctions; // global idx of junctions vertices in meshv
	dal::bit_vector extrema;   // global idx of extreme vertices in meshv

	Jv_HT.clear();
	nb_extrema=0; 
	nb_junctions=0;
	
	size_type fer = nb_branches; // first empty region
	GMM_ASSERT1(meshh.has_region(fer)==0, 
		"Overload in meshv region assembling!");
	
	// List all the convexes
	dal::bit_vector nn = meshh.convex_index();
	bgeot::size_type cv;
	for (cv << nn; cv != bgeot::size_type(-1); cv << nn) {
		
		bgeot::pconvex_structure cvs = meshh.structure_of_convex(cv);
		if (cvs->nb_points()>2) 
			cerr << "Error: convex #" << cv << "has more than 2 vertices!" << endl;
		if (cvs->nb_faces()>2)  
			cerr << "Error: convex #" << cv << "has more than 2 faces!" << endl;

		// Build regions for BCs and junctions
		// Global idx of mesh vertices
		size_type i0 = meshh.ind_points_of_convex(cv)[cvs->ind_points_of_face(1)[0]];
		size_type i1 = meshh.ind_points_of_convex(cv)[cvs->ind_points_of_face(0)[0]];
		// Identify vertex type
		if (meshh.convex_to_point(i0).size()==1){ /* inflow extremum */
			// Update information
			extrema.add(i0);
			nb_extrema++;
			// Build a new region made by a single face
			GMM_ASSERT1(meshh.has_region(fer)==0, 
				"Overload in meshv region assembling!");
			meshh.region(fer).add(cv, 1);
			// Store the current index and then update it
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_HT.size())) {
				found = (i0 == BCv_HT[bc].idx);
				if (!found) bc++;
			}
			GMM_ASSERT1(found=true, "Miss a boundary node in BCv list!");
			BCv_HT[bc].rg = fer; 
			fer++;
			// Store the containing branch index
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshh.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			BCv_HT[bc].branches.emplace_back(branch); 
		}
		else if (meshh.convex_to_point(i0).size()==2){ /* trivial inflow junction */
			// DO NOTHING
		}
		else if (meshh.convex_to_point(i0).size()>=2){ /* non-trivial inflow junction */
			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i0);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i0);
				nb_junctions++;
				GMM_ASSERT1(meshh.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshh.region(fer).add(cv, 1); // single-face region
				// Create a new junction node
				Jv_HT.emplace_back("JUN", 0, i0, fer);
				fer++;
			}
			// Search for index of containing branch (\mathcal{P}^{in}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshh.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");
			// Add the inflow branch (to the right junction node)
			size_type jj = 0;
			bool found = false;
			while (!found && jj < nb_junctions){
				found = (i0 == Jv_HT[jj].idx);
				if (!found) jj++;
			}
			//cout << "Branch -" << branch << " added to junction " << jj << endl;
			Jv_HT[jj].value += param.R(mimv, branch);
			Jv_HT[jj].branches.emplace_back(-branch);
			GMM_ASSERT1(branch>0, 
				"Error in network labeling: -0 makes no sense");
		}
		
		if (meshh.convex_to_point(i1).size()==1){ 
			size_type bc = 0; 
			bool found = false;
			while (!found && (bc<BCv_HT.size())) {
				found = (i1 == BCv_HT[bc].idx);
				if (!found) bc++;
			}
			if (found){ /* outlow extremum */
				extrema.add(i1); 
				nb_extrema++; 
				// Build a new region made by a single face
				GMM_ASSERT1(meshh.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshh.region(fer).add(cv, 0);
				// Store the current index and then update it
				BCv_HT[bc].value *= 1.0;
				BCv_HT[bc].rg = fer; 
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshh.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_HT[bc].branches.emplace_back(branch); 
			}
			else { /* interior -> Mixed point */
				// "MIX" label via post-processing
				// Build a new region made by a single face
				GMM_ASSERT1(meshh.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				meshh.region(fer).add(cv, 0);
				BCv_HT.emplace_back("MIX", 0.0, i1, fer);
				fer++;
				// Store the containing branch index
				size_type branch = 0; 
				bool contained = false;
				while (!contained && branch<nb_branches ) {
					contained = meshh.region(branch).is_in(cv);
					if (!contained) branch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				BCv_HT.back().branches.emplace_back(branch); 
			}
		}
		else if (meshh.convex_to_point(i1).size()==2){ /* trivial outflow junction */

			// Search for index of first containing branch (\mathcal{P}^{out}_j)
			size_type firstbranch = 0; 
			bool contained = false;
			while (!contained && firstbranch<nb_branches ) {
				contained = meshh.region(firstbranch).is_in(cv);
				if (!contained) firstbranch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i1!");

			// Check if i1 is a trivial junction (or a INT point)
			size_type cv1 = meshh.convex_to_point(i1)[0];
			size_type cv2 = meshh.convex_to_point(i1)[1];
			bool is_junc = (meshh.region(firstbranch).is_in(cv1) < 1 ||
							meshh.region(firstbranch).is_in(cv2) < 1 );
							
			if (is_junc){
				cout << "Found a trivial junction at i1 = " << i1 << endl;
				// Check if jucntion has been already stored, 
				// if not add to the junction list (J) and build a new region
				dal::bit_vector tmp; tmp.add(i1);
				if(!junctions.contains(tmp)){
					// Store the junction vertex
					junctions.add(i1);
					nb_junctions++;
					GMM_ASSERT1(meshh.has_region(fer)==0, 
						"Overload in meshv region assembling!");
					// Build a new region with idx "first empty region"
					meshh.region(fer).add(cv, 0);
					// Create a new junction node
					Jv_HT.emplace_back("JUN", 0, i1, fer);
					fer++;
				// Search for index of second containing branch (\mathcal{P}^{out}_j)
				size_type secondbranch = 0; 
				size_type secondcv = (( cv1 == cv) ? cv2 : cv1);
				size_type firstcv = (( cv1 != cv) ? cv2 : cv1);
				contained = false;
				while (!contained && secondbranch<nb_branches ) {
					if (secondbranch!=firstbranch)
					contained = meshh.region(secondbranch).is_in(secondcv);
					if (!contained) secondbranch++;
				}
				GMM_ASSERT1(contained=true, "No branch region contains node i1!");
				// Add the two branches
				scalar_type in;
				in=0;
				if (meshh.ind_points_of_convex(firstcv)[0]==i1) in=-1;
				else if (meshh.ind_points_of_convex(firstcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in firstbranch convex index");
				Jv_HT.back().branches.emplace_back(in*firstbranch);

				in=0;
				if (meshh.ind_points_of_convex(secondcv)[0]==i1) in=-1;
				else if (meshh.ind_points_of_convex(secondcv)[1]==i1) in=+1;
				GMM_ASSERT1(in!=0, "There's something wrong in secondbranch convex index");
				Jv_HT.back().branches.emplace_back(in*secondbranch);
				Jv_HT.back().value += param.R(mimv, firstbranch);
				Jv_HT.back().value += param.R(mimv, secondbranch);
				}
			}
		}
		else if (meshh.convex_to_point(i1).size()>=2){ /* non-trivial outflow junction */

			// Search for index of containing branch (\mathcal{P}^{out}_j)
			size_type branch = 0; 
			bool contained = false;
			while (!contained && branch<nb_branches ) {
				contained = meshh.region(branch).is_in(cv);
				if (!contained) branch++;
			}
			GMM_ASSERT1(contained=true, "No branch region contains node i0!");

			// Check if jucntion has been already stored, 
			// if not add to the junction list (J) and build a new region
			dal::bit_vector tmp; tmp.add(i1);
			if(!junctions.contains(tmp)){
				// Store the junction vertex
				junctions.add(i1);
				nb_junctions++;
				GMM_ASSERT1(meshh.has_region(fer)==0, 
					"Overload in meshv region assembling!");
				// Build a new region with idx "first empty region"
				meshh.region(fer).add(cv, 0);
				// Create a new junction node
				Jv_HT.emplace_back("JUN", 0, i1, fer);
				// Add the outflow branch
				Jv_HT.back().branches.emplace_back(+branch);
				Jv_HT.back().value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << i1 << endl;
				fer++;
			}
			else {
				// Add the outflow branch (to the right junction node)
				size_type jj = 0;
				bool found = false;
				while (!found && jj < nb_junctions){
					found = (i1 == Jv_HT[jj].idx);
					if (!found) jj++;
				}
				Jv_HT[jj].branches.emplace_back(+branch);
				Jv_HT[jj].value += param.R(mimv, branch);
				//cout << "Branch " << branch << " added to junction " << jj << endl;
			}
		}

	} /* end of convexes loop */
	
	// Ckeck network assembly
	#ifdef M3D1D_VERBOSE_
	cout << "--- NETWORK ASSEMBLY ------------------ "   << endl;
        /*cout << "  Branches:   " << nb_branches << endl
		 << "  Vertices:   " << nn.size()+1 << endl;
	cout << "  Extrema:    " << extrema << endl;	  
	for (size_type i=0; i<BCv_HT.size(); ++i)
		cout << "    -  label=" << BCv_HT[i].label 
			 << ", value=" << BCv_HT[i].value << ", ind=" << BCv_HT[i].idx 
			 << ", rg=" << BCv_HT[i].rg << ", branches=" << BCv_HT[i].branches << endl; 
	cout << "  Junctions: " << junctions << endl;
	for (size_type i=0; i<Jv_HT.size(); ++i)
                cout << "    -  label=" << Jv_HT[i].label
			 << ", value=" << Jv_HT[i].value << ", ind=" << Jv_HT[i].idx 
                         << ", rg=" << Jv_HT[i].rg << ", branches=" << Jv_HT[i].branches << endl;
        cout << "---------------------------------------- "   << endl;*/
	#endif

} 
GMM_STANDARD_CATCH_ERROR; // catches standard errors

} /* end of build_vessel_boundary */

//////// Assemble the problem ////////////////////////////////////////// 
void
problemHT::assembly(void)
{	
	//1 Build the monolithic matrix AM
	assembly_mat();
	//2 Build the monolithic rhs FM
	assembly_rhs();
}

void
problemHT::assembly_fixpoint(void)
{
assembly();
}

void 
problemHT::assembly_mat(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Allocating AM_HT, UM_HT, FM_HT ..." << endl;
	#endif
	gmm::resize(AM_HT, dofHT.H(), dofHT.H()); gmm::clear(AM_HT);
	gmm::resize(FM_HT, dofHT.H()); gmm::clear(FM_HT);
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling the monolithic matrix AM_HT ..." << endl;
	#endif
	// Junction compatibility matrix for the hematocrit problem
	sparse_matrix_type Jh(dofHT.H(), dofHT.H());
	sparse_matrix_type Jvv(dof.Pv(), dofHT.H());

	
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Bh and Jh ..." << endl;
	#endif

	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling the tangent versor (Hematocrit)..." << endl;
	#endif


	vector_type Uv( dof.Uv()); gmm::clear(Uv);
	gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())) ,  Uv);

	scalar_type Theta = PARAM.real_value("THETA", "Theta Number");
	vector_type element_size(dofHT.H());
	scalar_type max_size=0;
	scalar_type max_product=0;//Luca
	scalar_type temp;
	
	//Diffusivity Simone
	/*for(dal::bv_visitor k(meshv.convex_index()); !k.finished();++k){
		temp=meshv.convex_area_estimate(k,2);
		if(temp>max_size) max_size=temp;
		}
	scalar_type max_U;
	scalar_type max_U_positive=*max_element(Uv.begin(), Uv.end());
	scalar_type max_U_negative=*min_element(Uv.begin(), Uv.end());
		if(max_U_positive > fabs(max_U_negative))
			max_U=max_U_positive;
		else
			max_U=fabs(max_U_negative);
	
	scalar_type Diffusivity=max_size*max_U/2*Theta;*/

	//Diffusivity Luca
	#ifdef M3D1D_VERBOSE_
	cout << endl << "Assembling artificial diffusivity" << endl;
	#endif
	size_type shift_U = 0;
        size_type shift_H = 0;
        for(size_type i=0; i<nb_branches; ++i){
		#ifdef M3D1D_VERBOSE_
            	cout << "Branch " << i << endl;
		#endif
		//Estimate maximum h for the i-th branch
            	for(dal::bv_visitor k(meshv.region(i).index()); !k.finished();++k){
			temp=meshv.convex_area_estimate(k,2);
			if(temp>max_size) max_size=temp;
		}
		#ifdef M3D1D_VERBOSE_
            		cout << "Maximum h: " << max_size << endl;
		#endif		
		
		//Estimate maximum u for the i-th branch
            	if(i>0) shift_U += mf_Uvi[i-1].nb_dof();
            	//Obtain the vector of velocity in branch i
            	vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
            	gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift_U, mf_Uvi[i].nb_dof())) ,  Uvi);
            	//maximum u
            	scalar_type max_U;
            	scalar_type max_U_positive=*max_element(Uvi.begin(), Uvi.end());
           	scalar_type max_U_negative=*min_element(Uvi.begin(), Uvi.end());
            	if(max_U_positive > fabs(max_U_negative))
                     max_U=max_U_positive;
            	else
                     max_U=fabs(max_U_negative);
		#ifdef M3D1D_VERBOSE_
            		cout << "Maximum velocity: " << max_U << endl;
		#endif
            	temp = max_U * max_size;
            	if(temp> max_product) 
			max_product=temp;
        }
        scalar_type Diffusivity =  max_product * Theta/2;



	#ifdef M3D1D_VERBOSE_
	//cout << "Max Element Size	  : " << max_size << endl;
	//cout << "Max Velocity    	  : " << max_U << endl;
	cout << "Max (Velocity*h)    	  : " << max_product << endl;
	cout << "Artificial Diffusivity   : " << Diffusivity << endl;
	if(Diffusivity)
		{
			scalar_type Pe_h = max_product/2/Diffusivity;
			cout << "max(Pe_h)		  : " << Pe_h << endl;
		}
	#endif


	
	// Local matrices
	//size_type shift_U = 0;//luca
	//size_type shift_H = 0;
	shift_U = 0;
	shift_H = 0;//luca -dif
	for(size_type i=0; i<nb_branches; ++i){
		if(i>0) shift_U += mf_Uvi[i-1].nb_dof();
		if(i>0) shift_H += mf_Hi[i-1].nb_dof();
		//Obtain the vector of velocity in branch i
		sparse_matrix_type Bhi(mf_Hi[i].nb_dof(),mf_Hi[i].nb_dof());clear(Bhi);
		sparse_matrix_type Dhi(mf_Hi[i].nb_dof(),mf_Hi[i].nb_dof());clear(Dhi);

 		vector_type Uvi( mf_Uvi[i].nb_dof()); gmm::clear(Uvi);
		gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt()+shift_U, mf_Uvi[i].nb_dof())) ,  Uvi);

		//Obtain the radius of branch i
		scalar_type Ri = param.R(mimv, i);
//scalar_type Ri = param.Ri(i);
		vector_type R_veci(mf_coefvi[i].nb_dof(),Ri);
		//Obtain the flow in the branch i
		gmm::scale(Uvi,pi*Ri*Ri);
		// Allocate temp local tangent versor
			#ifdef M3D1D_VERBOSE_
		cout << "Assembling Advection Matrix for branch n° " << i << endl;
			#endif
		// Build Bhi
	
		asm_advection_hematocrit(Bhi, mimv, mf_Hi[i], mf_Uvi[i],
								mf_coefvi[i], Uvi, R_veci,
									param.lambdax(i), param.lambday(i), param.lambdaz(i), meshv.region(i));
		
		// cout << "---> --> --> --> versore tangente ematocrito ramo.."<< i << ": ..." << param.lambday(i) << endl;

			#ifdef M3D1D_VERBOSE_
		cout << "Assembling Artificial Viscosity Matrix for branch n° " << i << endl;
			#endif
		// Build Dhi

		vector_type diff(mf_coefvi[i].nb_dof(),Diffusivity);
		gmm::scale(diff,pi*Ri*Ri);
		asm_network_artificial_diffusion (Dhi, mimv, mf_Hi[i], mf_coefvi[i], diff, meshv.region(i));
		// Copy Bhi and Dhi
		gmm::scale(Dhi,1.0);
		gmm::scale(Bhi,-1.0);
		gmm::add(Bhi,
				gmm::sub_matrix(AM_HT, 
				gmm::sub_interval(shift_H, mf_Hi[i].nb_dof()),
				gmm::sub_interval(shift_H, mf_Hi[i].nb_dof())));
		


		gmm::add(Dhi,
				gmm::sub_matrix(AM_HT, 
				gmm::sub_interval(shift_H, mf_Hi[i].nb_dof()),
				gmm::sub_interval(shift_H, mf_Hi[i].nb_dof())));
		gmm::clear(Dhi);
		gmm::clear(Bhi); 

	} /* end of branches loop */

		asm_HT_out(AM_HT, mimv, mf_Hi, Uv, param.R(), mf_Uvi, mf_coefv);

			#ifdef M3D1D_VERBOSE_
		cout << "Assembling Hematocrit Junctions..."<< endl;
			#endif
		scalar_type dim = PARAM.real_value("d", "characteristic length of the problem [m]");
		dim=dim*1E6; // unit of measure in Pries formula is micrometers

		asm_hematocrit_junctions(Jvv, Jh,Uv, mimv,mf_Hi, mf_Pv, mf_Uvi,mf_coefv, Jv_HT, param.R(),UM_HT,dim, AM_HT);
 
		// Copy Jh
		gmm::add(Jh,AM_HT);



}

void 
problemHT::assembly_rhs(void)
{
	#ifdef M3D1D_VERBOSE_
	cout << "Assembling rhs of FM_HT ... " << endl;
	#endif

	#ifdef M3D1D_VERBOSE_
	cout << "  Initializing RHS for FM_HT ..." << endl;
	#endif

	// Coefficients for BCs
	scalar_type bcoef  = PARAM.real_value("BETA_H", "Coefficient for mixed BC of Ht");

	vector_type beta(dofHT.coefh(), bcoef);

	#ifdef M3D1D_VERBOSE_
	cout << "  Building hematocrit boundary term ..." << endl;
	#endif

 	vector_type Uv( dof.Uv()); gmm::clear(Uv);
	gmm::copy(gmm::sub_vector(UM, gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())),  Uv);

	asm_HT_bc (AM_HT, FM_HT, mimv, mf_Hi, mf_coefv, bcoef, BCv_HT, param.R());

	gmm::clear(beta);

}

////////// Solve the problem ///////////////////////////////////////////    

vector_type 
problemHT::iteration_solve(vector_type U_O,vector_type F_N){
	
	#ifdef M3D1D_VERBOSE_
	cout << "Solving the hematocrit system ... " << endl;
	#endif
	
	scalar_type alfa=descrHT.underH;
	gmm::csc_matrix<scalar_type> A_HT;
	gmm::clean(AM_HT, 1E-12);
	gmm::copy(AM_HT, A_HT);
	scalar_type cond;
	vector_type U_new;
	gmm::resize(U_new, dofHT.H()); gmm::clear(U_new);

	//Solving with SuperLU method
	gmm::SuperLU_solve(A_HT, U_new, F_N, cond);
	//cout << "Old Ht is " << gmm::sub_vector(U_O, gmm::sub_interval(dof.Ut(), dof.Pt())) << endl;

	//UNDER-RELAXATION
	if(alfa!=1){
	gmm::scale(U_new,alfa);
	gmm::scale(U_O,(1-alfa));
	gmm::add(U_O,U_new);
	gmm::scale(U_O,1/(1-alfa));}
	return U_new;


}

scalar_type
problemHT::calcolo_Rk(vector_type U_N, vector_type U_O){

// The residual is computed as ||V(k)-V(k-1)||/||V(k-1)|| with ||V(k)|| Eucledian norm

	vector_type Sol_old(dofHT.H()); // H(k-1) 
	vector_type Sol_new(dofHT.H());
	scalar_type R_H, N_H;
	scalar_type residual;
	//obtain H(k-1)
	gmm::copy(U_O,Sol_old);
	//obtain H(k)
	gmm::copy(U_N,Sol_new);

	gmm::scale(Sol_old, -1.0);

	gmm::add(Sol_old, Sol_new);

	R_H=gmm::vect_norm2(Sol_new);

	N_H=gmm::vect_norm2(Sol_old);
	
	residual=R_H/N_H;


	return residual; // gives if residual is bigger than the max value given in input 


}
bool
problemHT::solve_fixpoint(void)
{
/*solver 
1- Declaration of variables
2- Save the constant matrices (that don't change during the iterative process)
3- Get the intial guess for hematocrit (the fluid dynamic solution is already done in the main)
4- Iterative Process
			a-compute the viscosity in each vessel
			b-modify the mass matrix for fluid dynamic problem
			c- add the lymphatic contribution
			d-1 find the new solution as AM *U(k+1) = F(k)
			d-2 under-relaxation process U(k+1)= alfa*U(k+1) + (1-alfa)U(k)
			e-1 find the new solution for hematocrit as AM_HT *H(k+1) = F(k)
			e-2 under-relaxation process H(k+1)= alfa*H(k+1) + (1-alfa)H(k)
			f-compute TFR
			g-compute lymphatic total flow rate
			h-compute total FR going in or out the interstitial domain

			i- check residuals Rk: the iterative procedure ends when either maximum number of iteration is reached or the residuals of method are less than a pre-established value. 
			To this purpose THREE RESIDUALS ARE DEFINED:
				I-solution residual=||ut(k)-ut(k-1)||/||ut(k-1)||+||pt(k)-pt(k-1)||/||pt(k-1)||+||uv(k)-uv(k-1)||/||uv(k-1)||+||pv(k)-pv(k-1)||/||pv(k-1)|| < epsilon_f
				I-solution residual=||H(k)-H(k-1)||/||H(k-1)|| < epsilon_h
					with || V || Euclidean norm
				II-conservation mass residual \sum_i (Dtt*Ut(k)+Btt*Pt(k)-Btv*Pv(k)+Btv*DeltaPi(k+1)+Mlf*FI(k+1))/TFR < epsilon_m
			l- Update the value of U(k-1) with U(k)
			m- Update the value of H(k-1) with H(k)

						-Every N iteration the solution is saved


*/
// 1 - declaration of variables
	vector_type U_new; 
	vector_type U_old;
	gmm::resize(U_new, dof.tot()); gmm::clear(U_new);
	gmm::resize(U_old, dof.tot()); gmm::clear(U_old);
	vector_type F_new;
	gmm::resize(F_new, dof.tot()); gmm::clear(F_new);
	bool print_res=descr.print_residual;
	scalar_type epsSol=descr.epsSol;
	scalar_type resSol=epsSol*100;
	scalar_type epsCM=descr.epsCM;
	scalar_type resCM=epsCM*100;
	scalar_type max_iteration=descr.Max_it;
	int iteration_save=descr.Save_it;
	int iteration=0;
	bool RK=1;
	clock_t t;
	clock_t time_G;
	vector_type F_LF; gmm::resize(F_LF, dof.Pt());
	vector_type Uphi(dof.Pv()); 
	sparse_matrix_type Bvt(dof.Pv(), dof.Pt());
	sparse_matrix_type Bvv(dof.Pv(), dof.Pv());
	sparse_matrix_type Btv(dof.Pt(), dof.Pv());
	vector_type Pt(dof.Pt()); 
	vector_type Pv(dof.Pv()); 
	scalar_type Pi_t=param.pi_t();
	scalar_type Pi_v=param.pi_v();
	scalar_type sigma=param.sigma();
	scalar_type picoef;
        vector_type DeltaPi, Ones(dof.Pv(),1.0);
	gmm::resize(DeltaPi, dof.Pv()); gmm::clear(DeltaPi);
        vector_type auxOSt(dof.Pt());
        vector_type auxOSv(dof.Pv());
        vector_type auxCM(dof.Pt()); gmm::clear(auxCM);
	// Gnuplot gp;
	vector_type RES_SOL(max_iteration), RES_CM(max_iteration);

	//Hematocrit variables
	vector_type H_old(dofHT.H());
	vector_type H_new(dofHT.H());
	scalar_type dim = PARAM.real_value("d", "characteristic length of the problem [m]");
	dim=dim*1E6;
	scalar_type mu_plasma=paramHT.visco_plasma();
	int visco_v=paramHT.visco_type();
	scalar_type mu_start = PARAM.real_value("mu_v", "blood viscosity [kg/ms]"); 
	vector_type RES_H(max_iteration);
	scalar_type epsH=descrHT.epsH;
	scalar_type resH=epsH*100;
// 2 - Saving the constant matrices
	#ifdef M3D1D_VERBOSE_
	cout << "Saving the constant matrices ... " << endl;
	#endif
	//Extracting matrices Bvt, Bvv
	gmm::copy(gmm::sub_matrix(AM, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv()	, dof.Pv()),
			gmm::sub_interval(dof.Ut(), dof.Pt())),
				Bvt); 
	gmm::copy(gmm::sub_matrix(AM, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv()), 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
				Bvv);

	//Extracting matrix Btv
	gmm::copy(gmm::sub_matrix(AM, 
			gmm::sub_interval(dof.Ut(), dof.Pt()),
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())),
				Btv);
	gmm::scale(Btv,-1.0);
	//Extracting Oncotic term				
	picoef=sigma*(Pi_v-Pi_t);
	gmm::copy(Ones, DeltaPi);
        gmm::scale(DeltaPi,picoef);
       	gmm::mult(Btv,DeltaPi,auxOSt);
        gmm::mult(Bvv,DeltaPi,auxOSv);

	//Extracting Mvv_kv
	#ifdef M3D1D_VERBOSE_
	cout << "  Assembling Mvv0 in FixPoint Hematocrit..." << endl;
	#endif
	// Local matrices
	size_type shift = 0;
	sparse_matrix_type Mvv_mu(dof.Uv(), dof.Uv());
	for(size_type i=0; i<nb_branches; ++i){

		if(i>0) shift += mf_Uvi[i-1].nb_dof();
//scalar_type Ri = param.Ri(i);
		scalar_type Ri = param.R(mimv, i);
// cout << " ----------raggio.." <<Ri << endl;
//scalar_type kvi = param.kv(i);
		 scalar_type kvi = param.kv(mimv, i);
// cout << " ----------kvi.." <<kvi << endl;

		/* Va questo nel caso curvo?
		// Coefficient \pi^2*Ri'^4/\kappa_v
		vector_type ci(mf_coefvi[i].nb_dof(), pi*pi*Ri*Ri*Ri*Ri/kvi);
		*/ //o questo?
		vector_type ci(mf_coefvi[i].nb_dof()); //gmm::clear(ci);
		for(size_type j=0; j<mf_coefvi[i].nb_dof(); ++j){
			ci[j]=pi*pi*Ri*Ri*Ri*Ri/kvi*(1.0+param.Curv(i,j)*param.Curv(i,j)*Ri*Ri);
// 			cout<<" ci"<<ci[j];
// 			cout<<" Ri"<<Ri;
// 			cout<<" kvi"<<kvi;
// 			cout<<" curv"<<param.Curv(i,j)<<endl;

		}
// 		cout<<"\n\n\n ci="<<ci[0]<<"    u="<<3.5/ci[0]*pi*Ri*Ri<<"\n\n\n";
		// Allocate temp local matrices
		sparse_matrix_type Mvvi(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		sparse_matrix_type Dvvi(dof.Pv(), mf_Uvi[i].nb_dof());
		// Build Mvvi and Dvvi
		asm_network_poiseuille(Mvvi, Dvvi, 
			mimv, mf_Uvi[i], mf_Pv, mf_coefvi[i],
			ci, param.lambdax(i), param.lambday(i), param.lambdaz(i), meshv.region(i));
		gmm::scale(Dvvi, pi*Ri*Ri);
		// Copy Mvvi and Dvvi
		gmm::add(Mvvi, 
			gmm::sub_matrix(Mvv_mu, 
				gmm::sub_interval(shift, mf_Uvi[i].nb_dof()), 
				gmm::sub_interval(shift, mf_Uvi[i].nb_dof()))); 
		gmm::clear(Dvvi);
		gmm::clear(Mvvi);
		
	} /* end of branches loop */
		sparse_matrix_type Mvv_bc(dof.Uv(),dof.Uv());
		sparse_matrix_type Mvv(dof.Uv(),dof.Uv());
		gmm::copy(gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())),Mvv);
		gmm::scale(Mvv_mu,-1.0);
		gmm::add(Mvv,Mvv_mu,Mvv_bc);
		gmm::clear(Mvv_mu);
		gmm::clear(Mvv);

	// Opening file to save number of iteration and residual
	std::ofstream SaveResidual;
	SaveResidual.open(descr.OUTPUT+"Residuals.txt");
	SaveResidual << "Iteration" << "\t" << "Solution Residual" << "\t" << "Mass Conservation Residual" << "\t" << "Hematocrit Residual" << endl;
	gmm::copy(UM,U_old);

	time_G=clock();

	scalar_type H_start = PARAM.real_value("H_START", "hematocrit start");
	gmm::resize(UM_HT,dofHT.H()); gmm::clear(UM_HT);
	vector_type ones_H(dofHT.H(),1.00);
	gmm::add(ones_H, UM_HT);
	gmm::scale(UM_HT,H_start);	

		assembly();

// 3 - Get the initial guess H0
	#ifdef M3D1D_VERBOSE_
	cout << "Solving the hematocrit system ... " << endl;
	#endif

	gmm::csc_matrix<scalar_type> A_HT;
	gmm::clean(AM_HT, 1E-12);
	gmm::copy(AM_HT, A_HT);
	scalar_type cond;
	//Solving with SuperLU method
	gmm::SuperLU_solve(A_HT, UM_HT, FM_HT, cond);

	#ifdef M3D1D_VERBOSE_
	cout << "Solved the initial guess for hematocrit" << endl;
	#endif

	gmm::copy(UM_HT,H_old);
//4- Iterative Process
while(RK && iteration < max_iteration)
	{	

// a-compute the viscosity in each vessel
	#ifdef M3D1D_VERBOSE_
	cout << "Computing Viscosity - Iteration "<< iteration << "..." << endl;
	#endif
		shift = 0;
		size_type shift_h=0;
		scalar_type shift_coef=0;
		gmm::clear(MU); gmm::resize(MU, mf_coefv.nb_dof());
	for(size_type i=0; i<nb_branches; ++i){

		vector_type Hi(mf_Hi[i].nb_dof());
		vector_type H_const(mf_coefvi[i].nb_dof());
		scalar_type Ri = param.R(mimv, i);
		vector_type mui; gmm::clear(mui);


		if(i>0) shift_h += mf_Hi[i-1].nb_dof();
		if(i>0) shift_coef += mf_coefvi[i-1].nb_dof();
		if(i>0) shift += mf_Uvi[i-1].nb_dof();


		gmm::copy(gmm::sub_vector(H_old, 
			gmm::sub_interval(shift_h, mf_Hi[i].nb_dof())), Hi);

// cout << "-------- inizia interpolation" << endl;
		getfem::interpolation(mf_Hi[i], mf_coefvi[i], Hi, H_const, 0);


		switch(visco_v)
				{
				 case(0):{ //cout << "-------- case 0 " << endl;
					for (auto h : H_const)
						{
							if(h==0)
							{mui.emplace_back(mu_plasma);
				// cout << "-------- case 0  ---- h==0 " << endl;
							}
							else
							{mui.emplace_back(viscosity_vivo(h, Ri*dim, mu_plasma));
				// cout << "-------- case 0  ---- h!=0 " << endl;
				// cout << "-------- mui  ---- h!=0 "<< mui << endl; 
							}
						}
					}break;
				 case(1):{ //cout << "-------- case 1 " << endl;

					for (auto h : H_const)
						{
							if(h==0)
							mui.emplace_back(mu_plasma);
							else
							mui.emplace_back(viscosity_vitro(h, Ri*dim, mu_plasma));
						}
						}break;
				default:
					cerr << "Invalid value for Visco_v " << visco_v << endl;
				}
// cout << "-------- end switch " << endl;
			size_type pos=0;
			for (getfem::mr_visitor mrv(mf_coefv.linked_mesh().region(i)); !mrv.finished(); ++mrv)
			for (auto muu : mf_coefv.ind_basic_dof_of_element(mrv.cv()))
				{
				MU[muu] = mui[pos];
				pos++;}

//b-modify the mass matrix for fluid dynamic problem
	#ifdef M3D1D_VERBOSE_
	cout << "Modify Mvvk - Iteration "<< iteration << "..." << endl;
	#endif
		scalar_type kvi = param.kv(mimv, i);
		// Coefficient \pi^2*Ri'^4/\kappa_v
		vector_type ci(mf_coefvi[i].nb_dof()); //gmm::clear(ci);
		for(size_type j=0; j<mf_coefvi[i].nb_dof();j++)
			{ci[j]=pi*pi*Ri*Ri*Ri*Ri/kvi*(1.0+param.Curv(i,j)*param.Curv(i,j)*Ri*Ri)/mu_start*mui[j];
// cout << "-------- ci  "<<ci[j]<< " ";
// 			cout<<" Ri"<<Ri;
// 			cout<<" kvi"<<kvi;
// 			cout<<" curv"<<param.Curv(i,j)<<endl;
// 	cout << "mu_start" << mu_start << endl;
// cout << "mui[j]" <<mui[j] << endl;
			}
// 		cout<<"\n\n\n ci="<<ci[0]<<"    u="<<3.5/ci[0]*pi*Ri*Ri<<"\n\n\n";
		// Allocate temp local matrices
// cout << "-------- entra Mvv_mui "<< endl;
		sparse_matrix_type Mvv_mui(mf_Uvi[i].nb_dof(), mf_Uvi[i].nb_dof());
		// Build Mvv_mui
// cout << "-------- entra netw_pois "<< endl;		
		asm_network_poiseuilleHT(Mvv_mui, mimv, mf_Uvi[i], mf_coefvi[i], ci, meshv.region(i));
		// Copy Mvv_mui in Mvv_mu
		gmm::add(Mvv_mui, 
			gmm::sub_matrix(Mvv_mu, 
				gmm::sub_interval(shift, mf_Uvi[i].nb_dof()), 
				gmm::sub_interval(shift, mf_Uvi[i].nb_dof()))); 
		gmm::clear(Mvv_mui);
	} /* end of branches loop */
	//Update Mvv and AM matrix
		gmm::add(Mvv_mu,Mvv_bc,Mvv);
		gmm::clear(gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
		gmm::copy(Mvv,
				gmm::sub_matrix(AM, 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv()), 
				gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));
		gmm::clear(Mvv);
		gmm::clear(Mvv_mu);

//c- add the lymphatic contribution
	#ifdef M3D1D_VERBOSE_
	cout << "Adding Lymphatic Contribution - Iteration "<< iteration << "..." << endl;
	#endif
	gmm::copy(FM,F_new);
	if(!LINEAR_LYMPH()){
	//Adding lymphatic contribution
	F_new=problem3d1d::modify_vector_LF(U_old,F_new);
        //gmm::copy(F_new,FM);
	}

//d-1 find the new solution as AM *U(k+1) = F(k)
//d-2 under-relaxation process U(k+1)= alfa*U(k+1) + (1-alfa)U(k)
	#ifdef M3D1D_VERBOSE_
	cout << "Solving the fluid dynamic problem - Iteration "<< iteration << "..." << endl;
	#endif

	t=clock();

	U_new=problem3d1d::iteration_solve(U_old,F_new);
	gmm::copy(U_new,UM);
	
// 	gmm::clear(UM);
// 	problem3d1d::solve_samg();
// 	gmm::copy(UM,U_new);
	

	t=clock()-t;
//e-1 find the new solution for hematocrit as AM_HT *H(k+1) = F(k)
//e-2 under-relaxation process H(k+1)= alfa*H(k+1) + (1-alfa)H(k)
	#ifdef M3D1D_VERBOSE_
	cout << "Solving the hematocrit problem - Iteration "<< iteration << "..." << endl;
	#endif
		assembly();
		H_new=iteration_solve(H_old, FM_HT);

//f-compute TFR
//g-compute lymphatic total flow rate
//h-compute total FR going in or out the interstitial domain
	#ifdef M3D1D_VERBOSE_
	cout << "Computing Flow Rate - Iteration "<< iteration << "..." << endl;
	#endif
		// Extracting solutions Pt, Pv 
		gmm::copy(gmm::sub_vector(U_new, 
			gmm::sub_interval(dof.Ut(), dof.Pt())), Pt);
		gmm::copy(gmm::sub_vector(U_new, 
			gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(), dof.Pv())), Pv);
		// Lymphatic Contribution
		if(!LINEAR_LYMPH())
		F_LF=problem3d1d::compute_lymphatics(U_new);
		else
		{       
		vector_type Pl(dof.Pt(),PARAM.real_value("PL"));
		vector_type Pl_aux(dof.Pt());
		sparse_matrix_type Mlf (dof.Pt(), dof.Pt());
		scalar_type lf_coef=param.Q_LF(0);//scalar then uniform untill now
		asm_tissue_lymph_sink(Mlf, mimt, mf_Pt);
		gmm::scale(Mlf,lf_coef);
		gmm::mult(Mlf,Pl,Pl_aux);
		gmm::scale(Pl_aux,-1);
		gmm::mult(Mlf,Pt,F_LF);
		gmm::add(F_LF,Pl_aux,F_LF);
		}
		// Computing Bvv*Pv - Bvt*Pt
		gmm::mult(Bvt, Pt, Uphi);
		gmm::mult_add(Bvv, Pv, Uphi);
        	//oncotic term
		picoef=sigma*(Pi_v-Pi_t);
		gmm::copy(Ones, DeltaPi);
		gmm::scale(DeltaPi,-1.0);
        	gmm::scale(DeltaPi,picoef);
		gmm::mult_add(Bvv, DeltaPi, Uphi);
		//Computing TFR
		TFR = std::accumulate(Uphi.begin(), Uphi.end(), 0.0);
        	//computing flowrate of lymphatic system
		FRlymph = std::accumulate(F_LF.begin(), F_LF.end(), 0.0);
        	//computing flowrate from the cube
        	FRCube = TFR - FRlymph;

		if(RK && iteration < max_iteration && print_res && (iteration % iteration_save) == 0)
				{
				problem3d1d::export_vtk();
				export_vtk();
				cout << "Solution at iteration " << iteration+1 << " saved" << endl;
				cout << "TFR                 = " << TFR << endl;
				cout << "Lymphatic Flow Rate = " << FRlymph << endl;
				cout << "Flow Rate of cube   = " << FRCube << endl;
				}
//i- check residuals Rk
	#ifdef M3D1D_VERBOSE_
	cout << "Checking Residuals - Iteration "<< iteration << "..." << endl;
	#endif

	//Solution residual
			resSol=problem3d1d::calcolo_Rk(U_new, U_old);
	//Hematocrit residual
			resH=calcolo_Rk(H_new,H_old);
	//Conservation of mass residual
			gmm::mult(gmm::sub_matrix(AM, 
						gmm::sub_interval(dof.Ut(), dof.Pt()),
						gmm::sub_interval(0, dof.tot())),
						U_new,
							auxCM);
			gmm::add(auxOSt,auxCM);

			if(!LINEAR_LYMPH())
			gmm::add(F_LF,auxCM);

			scalar_type resCM=0;	
			if(TFR!=0)	
			resCM=std::accumulate(auxCM.begin(), auxCM.end(), 0.0)/TFR;

	RK=resSol>epsSol || fabs(resCM) > epsCM || resH > epsH; // all the residual must reach convergence to exit the "while"

	iteration++;

	//Saving residual values in an output file
	SaveResidual << iteration << "\t" << resSol << "\t" << resCM << "\t" << resH << endl;

			if(print_res)  {
			cout << "\nStep n°:" << iteration << "\nSolution Residual = " << resSol << "\nMass Residual = " << fabs(resCM) << "\nHematocrit Residual "<< resH << endl;
			cout << "\t\t\t\tTime: " <<  ((float)t)/CLOCKS_PER_SEC << " s "<< endl;
					}
			cout << "********************************************************" << endl;

//l- Update the value of U(k-1) with U(k)
//m- Update the value of H(k-1) with H(k)
	#ifdef M3D1D_VERBOSE_
	cout << "Updating Solution - Iteration "<< iteration << "..." << endl;
	#endif

	gmm::copy(U_new,U_old);
	gmm::copy(H_new,H_old);
	gmm::copy(H_old,UM_HT);

	//plotting residuals
	RES_SOL[iteration-1]=fabs(resSol);
	RES_CM[iteration-1]=fabs(resCM);
	RES_H[iteration-1]=fabs(resH);
// 	gp << "set logscale y; set xlabel 'iteration';set ylabel 'residual'; plot '-' w lines title 'Solution Residual', '-' w lines title 'Mass Conservation Residual','-' w lines title 'Hematocrit Residual'\n";
// 
// 	gp.send1d(RES_SOL);
// 	gp.send1d(RES_CM);
// 	gp.send1d(RES_H);
// 	gp.flush();

	//De-allocate memory
	gmm::clear(F_LF);
	gmm::clear(U_new); gmm::clear(auxCM); gmm::clear(F_new); gmm::clear(H_new);
	} //Exit the while
	
	gmm::copy(U_old,UM);

	time_G=clock()-time_G;
	cout<< "Iterative Process Time = " << ((float)time_G)/CLOCKS_PER_SEC << " s"<< endl;
	SaveResidual.close();
	if(RK)
		cout << "The method has NOT reached convergence for minimum residual" << endl;

// 	Gnuplot gp1;
// 	gp1 << "plot '-' w lines title 'Hematocrit'\n";
// 	gp1.send1d(UM_HT);
// 	gp1.flush();
/*	
	//De-allocate memory
	gmm::clear(auxOSt);
	gmm::clear(auxOSv);

// 	Gnuplot gp2;
// 	gp2 << "plot '-' w lines title 'Vessel velocity'\n";
// 	gp2.send1d(gmm::sub_vector(UM,gmm::sub_interval(dof.Ut()+dof.Pt(), dof.Uv())));   
// 
// 	gp2.flush();
	//De-allocate memory
	gmm::clear(auxOSt);
	gmm::clear(auxOSv);


	Gnuplot gp3;
	gp3 << "plot '-' w lines title 'Vessel pressure'\n";
	gp3.send1d(gmm::sub_vector(UM,gmm::sub_interval(dof.Ut()+dof.Pt()+dof.Uv(),dof.Pv())));   

	gp3.flush();
	//De-allocate memory
	gmm::clear(auxOSt);
	gmm::clear(auxOSv);*/

return true;
}

////////// Export results into vtk files ///////////////////////////////
void 
problemHT::export_vtk(const string & suff) //ODIFCA 
{
  if (PARAM.int_value("VTK_EXPORT"))
  {
	#ifdef M3D1D_VERBOSE_
	cout << "Exporting the solution (vtk format) to " << descr.OUTPUT << " ..." << endl;
	#endif
	 
	#ifdef M3D1D_VERBOSE_
	cout << "  Exporting Ht ..." << endl;
	#endif
	size_type start = 0;
	size_type length = 0;
	size_type k, last_dof;
	last_dof=nb_branches*mf_Hi[nb_branches-1].nb_dof();
	for(size_type i=0; i<nb_branches; ++i){
		if(i>0) start += mf_Hi[i-1].nb_dof();
		length = mf_Hi[i].nb_dof();
		vtk_export exp_Ht(descr.OUTPUT+"Ht"+suff+std::to_string(i)+".vtk");
			exp_Ht.exporting(mf_Hi[i]);
			exp_Ht.write_mesh();
			exp_Ht.write_point_data(mf_Hi[i], 
			gmm::sub_vector(UM_HT, gmm::sub_interval(start, length)), "Ht"); 
			}

	getfem::vtk_export expMU(descr.OUTPUT+"MU.vtk");
	expMU.exporting(mf_coefv);
	expMU.write_mesh();
	expMU.write_point_data(mf_coefv, MU, "mu");

	#ifdef M3D1D_VERBOSE_
	cout << "... export done, visualize the data file with (for example) Paraview " << endl; 
	#endif

  }
} /* end of export_vtk */



} /* end of namespace */
