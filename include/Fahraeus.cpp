/* -*- c++ -*- (enables emacs c++ mode) */
/*======================================================================
    "Mixed Finite Element Methods for Coupled 3D/1D Fluid Problems"
                      Politecnico di Milano
                          A.Y. 2016-2017
                  
                Copyright (C) 2016 Simone Di Gregorio
======================================================================*/
/*! 
  @file   	Fahraeus.cpp
  @authors 	Luca Possenti <luca.possenti@polimi.it>
%		Simone Di Gregorio <simone.digre@gmail.com>
  @date   	March 2017
  @brief  	Definition of Fahraeus effects routines.
 */
#include <Fahraeus.hpp>

namespace getfem {

// Compute the viscosity of plasma (see Biomachines course)
scalar_type
viscosity_plasma(scalar_type temperature)
{
	scalar_type viscosity;

	viscosity = 1.808 / ( 1 + 0.0337 * temperature + 0.00022 * temperature * temperature) * 1E-3 * 1.8;

	return viscosity;


}

// Compute the viscosity of blood in vitro (Pries et al 1992)
scalar_type 
viscosity_vitro(scalar_type hematocrit, scalar_type R, scalar_type mu_plasma)
{

		scalar_type visco_ast_45,e,viscosity,relative_viscosity;
		scalar_type D; //diameter of vessel
		scalar_type C;

		D = R*2;
		e=2.7182818284;

		visco_ast_45 = 220*pow(e,-1.3*D)+3.2-2.44*pow(e,-0.06*pow(D,0.645));

		C = (0.8+pow(e,-0.075*D))*(-1+1/(1+10*pow(D/10,12)))+1/(1+10*pow(D/10,12));

		relative_viscosity =  1 + ( visco_ast_45 - 1 ) * (pow((1-hematocrit),C)-1) / (pow((1-0.45),C)-1);

		viscosity = relative_viscosity*mu_plasma;

		return viscosity;

}

// Compute the viscosity of blood in vivo (Pries et al 1994)
scalar_type 
viscosity_vivo(scalar_type hematocrit, scalar_type R, scalar_type mu_plasma)
{
		
		scalar_type visco_ast_45,e,viscosity,relative_viscosity;
		scalar_type D; //diameter of vessel
		scalar_type C;

		D = R*2;
		e = 2.7182818284;
		/*cout << "diametro" << D << endl;
		cout << "mu_plasma" << mu_plasma << endl;
		cout << "hematocrit" << hematocrit << endl;*/

		visco_ast_45 = 6*pow(e,-0.085*D)+3.2-2.44*pow(e,-0.06*pow(D,0.645));

		//cout << "viscoast" << visco_ast_45 << endl;

		C = (0.8+pow(e,-0.075*D))*(-1+1/(1+10*pow(D/10,12)))+1/(1+10*pow(D/10,12));
		//cout << "C" << C << endl;
		relative_viscosity =  (1 + ( visco_ast_45 - 1 ) * (pow((1-hematocrit),C)-1) / (pow((1-0.45),C)-1)*pow(D/(D-1.1),2))*pow(D/(D-1.1),2);

		//cout << "relative_viscosity" << relative_viscosity << endl;
		viscosity = relative_viscosity*mu_plasma;
		//cout << viscosity << endl;

		//Check for admissible value of Viscosity (Remember that viscosity is expressed in Pa*s)
		if(viscosity > 0.05) //That is 50 cP (Blood viscosity in a vessel with radius 4 \mum and H=0.45: 9.33 cP; 1cP = 10^-3 Pa*s) 
			cout << "------------------------------------" << endl << "WARNING! Viscosity is higher than 50cP in at least one vessel. Please check your parameters" << endl <<  "------------------------------------"<< endl;
		return viscosity;

}
//Compute the fractional RBC flow in a bifurcation
scalar_type
fractional_Erythrocytes(scalar_type FQB, scalar_type D_f, scalar_type D, scalar_type D_2, scalar_type H_f)
{
/*cout << "FQB " <<FQB << endl;
cout << "D " <<D << endl;
cout << "D-f " <<D_f << endl;
cout << "D2 " <<D_2 << endl;
cout << "H " <<H_f << endl;*/
scalar_type FQE;
scalar_type A=PS_compute_A(D,D_2, H_f, D_f);
//cout << "A" << A << endl;
scalar_type B=PS_compute_B(H_f, D_f);
//cout << "B" << B << endl;
scalar_type X0=PS_compute_X0(H_f, D_f);
//cout << "x0" << X0 << endl;
scalar_type aux,logitFQE;

if(FQB <= X0)
FQE=0;
else if (FQB < (1 - X0) && FQB > X0){

aux=(FQB-X0)/(1-2*X0);
logitFQE=A+B*logit(aux);
FQE=pow(2.71828,logitFQE)/(1+pow(2.71828,logitFQE));

}
else if (FQB >= (1-X0))
FQE=1;
else
GMM_ASSERT1(0, "Error in SP computation " <<  endl);
//cout << "FQE" << FQE << endl;
return FQE;
}

scalar_type PS_compute_A(scalar_type D,scalar_type D_2, scalar_type H_f, scalar_type D_f)
{

scalar_type A=-13.29*((pow(D/D_2,2)-1)/(pow(D/D_2,2)+1))*(1-H_f)/D_f;

return A;
}
scalar_type PS_compute_B(scalar_type H_f, scalar_type D_f)
{

scalar_type B=1+6.98*(1-H_f)/D_f;

return B;
}
scalar_type PS_compute_X0(scalar_type H_f, scalar_type D_f)
{

scalar_type X0=0.964*(1-H_f)/D_f;

return X0;
}

scalar_type logit(scalar_type x)
{
scalar_type result=log(x/(1-x));

return result;
}

bool isPositive (scalar_type i) { return i>0;};


} /* end of namespace */
