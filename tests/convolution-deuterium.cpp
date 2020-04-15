#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "specdeut.h"
#include "constants.h"
#include "F2.h"

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

using namespace std;

// Step size for integration
const double dTheta = 0.05;

// Initialize spectral function class
specdeut * kaptari_sf = new specdeut();
F2 * F2p = new F2(0);
F2 * F2d = new F2(1);

// Spec functions for n/p in d:
map<double,double> mom_dist = kaptari_sf->getMomDist();
// Iterators to loop over:
map<double,double>::iterator itP;


// Functions with parameters to be minimized
double offshell_const( double virt, double xB , double off_a );
double offshell_lin( double virt, double xB , double off_a , double off_b );
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c );

// Function to create theory point
void calc_theo( double x, double Q2, double &theo_d_const , double &theo_d_lin , 
		const double *pars_const , const double *pars_lin ,double mom_cut );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [OutputTextFile]\n";
		return -1;
	}
	

	ofstream outfile;
	outfile.open(argv[1]);


	// Resulting constant offshell_const fit from SF
	const double pars_cstOff_SF[5] 			= {0.641326,2.3434,0.436479, -1.25864,		1};
	const double pars_linOff_xshift_SF[6] 		= {0.511169,2.13985,0.388139,-0.269376,-4.33521,1};
	const double pars_no_cstOff_SF[5] 			= {0.641326,2.3434,0.436479, 0,		1};
	const double pars_no_linOff_xshift_SF[6] 		= {0.511169,2.13985,0.388139,0,0,	1};

	for( double x = 0.2 ; x <= 0.96 ; x+=0.001 ){
		double Q2 = 14.*x;
		double h2_const = 0;
		double h2_lin = 0;
		cout << x << " ";

		// Deuterium convolution with offshell and no cut
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution no offshell and no cut
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution with and without offshell, 50MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.05);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.05);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution with and without offshell, 100MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.1);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.1);
		cout << h2_const << " " << h2_lin << " ";


		// Deuterium convolution with and without offshell, 200MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.2);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.2);
		cout << h2_const << " " << h2_lin << " ";
		
		// Deuterium convolution with and without offshell, 400MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.4);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.4);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution with and without offshell, 500MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.5);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.5);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution with and without offshell, 500MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.6);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.6);
		cout << h2_const << " " << h2_lin << " ";

		// Deuterium convolution with and without offshell, 500MeV
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_cstOff_SF,		pars_linOff_xshift_SF,		0.75);
		cout << h2_const << " " << h2_lin << " ";
		calc_theo(x,	Q2,	h2_const,	h2_lin,	pars_no_cstOff_SF,	pars_no_linOff_xshift_SF,	0.75);
		cout << h2_const << " " << h2_lin << " ";

		cout << F2p->Eval(x,Q2) << " " << 
			f2nf2p(x,pars_no_cstOff_SF[0],pars_no_cstOff_SF[1],pars_no_cstOff_SF[2]) << " " << 
			f2nf2p(x,pars_no_linOff_xshift_SF[0],pars_no_linOff_xshift_SF[1],pars_no_linOff_xshift_SF[2]) 
		<< "\n";
	}
	outfile.close();
	
	return 0;
}

double offshell_const( double virt, double xB , double off_a  ){
	return 1 + (off_a)*virt*virt;
}
double offshell_lin( double virt, double xB , double off_a , double off_b ){
	return 1 + (off_a + off_b*(xB-0.5))*virt*virt;
	//return 1 + (off_b*(xB-off_a))*virt*virt;
}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;
}

void calc_theo( double x, double Q2, double &theo_d_const , double &theo_d_lin , 
		const double *pars_const , const double *pars_lin ,double mom_cut ){
	int A 		= 2;
	double N,Z;
	theo_d_const = 0;
	theo_d_lin = 0;
	for( itP = mom_dist.begin() ; itP != mom_dist.end() ; itP++){
		double p_m = itP->first; // GeV
		double E_i = mD - sqrt(mN*mN + p_m*p_m);
		if( mom_cut != 0){ if( p_m < mom_cut ) continue; }

		// Grab momentum distribution value:
		double sp = itP->second; // GeV^-3

		for( double theta = 0 ; theta <= M_PI ; theta += dTheta){
			// Define y and nu from (E_i,p_i):
			double y = (E_i + p_m*cos(theta) )/mP;
			double nu = (E_i*E_i - p_m*p_m - mP*mP)/(mP*mP); //unitless

			double jacobian = sin(theta); // p^2 dp already inside momentum function normalization
			double phi_int = 2*M_PI;
			//double flux_fact = 1. + (p_m*cos(theta))/sqrt(p_m*p_m + mP*mP);
			double flux_fact = 1. + (p_m*cos(theta))/E_i;
			double mass_fact = mD/(mP*A);

			// Flag problematic kinematics for the delta conservation (?):
			if( y < x ) continue;
			if( y > A) continue;
			if( nu > 0. ) continue;

			Z = 1; N = 1;
			theo_d_const 
				+= jacobian * flux_fact * phi_int * mass_fact * dTheta
				* ( Z*sp*offshell_const(nu,x/y,pars_const[3])
				   + N*sp*offshell_const(nu,x/y,pars_const[3])*f2nf2p(x/y,pars_const[0],pars_const[1],pars_const[2]) )
				* F2p->Eval(x/y,Q2)
				* pars_const[4];
				
			theo_d_lin 
				+= jacobian * flux_fact * phi_int * mass_fact * dTheta
				* ( Z*sp*offshell_lin(nu,x/y,pars_lin[3],pars_lin[4])
				   + N*sp*offshell_lin(nu,x/y,pars_lin[3],pars_lin[4])*f2nf2p(x/y,pars_lin[0],pars_lin[1],pars_lin[2]) )
				* F2p->Eval(x/y,Q2)
				* pars_lin[5];

		}// end loop over theta
	}// end loop over p miss
	theo_d_const *= (1.) / ( F2p->Eval(x,Q2) * (1+f2nf2p(x,pars_const[0],pars_const[1],pars_const[2]) ) );
	theo_d_lin   *= (1.) / ( F2p->Eval(x,Q2) * (1+f2nf2p(x,pars_lin[0],pars_lin[1],pars_lin[2]) ) );
}
