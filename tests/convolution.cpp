#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>
//#include <Eigen/Dense>

#include "spec.h"
#include "constants.h"
#include "F2.h"

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

using namespace std;

// Step size for integration
const double dTheta = 0.05;

// Initialize spectral function class
spec * kaptari_sf = new spec();
F2 * F2p = new F2(0);
F2 * F2d = new F2(1);

// Spec functions for n/p in He-3
map<double,map<double,double>> test_n = kaptari_sf->getFullN();
map<double,map<double,double>> test_p = kaptari_sf->getFullP();
// 	iterators to loop over 
map<double,map<double,double>>::iterator itK_n;
map<double,map<double,double>>::iterator itK_p;
map<double,double>::iterator itE_n;
map<double,double>::iterator itE_p;
vector<double> alphas = kaptari_sf->getAlphas();
vector<double> nus = kaptari_sf->getNus();
vector<double> rho_p = kaptari_sf->getRhoProtons();
vector<double> rho_n = kaptari_sf->getRhoNeutrons();

// Functions with parameters to be minimized
double offshell( double virt, double xB , double off_a );
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c );

// Function to create theory point
void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut );
void calc_theo_rho( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [OutputTextFile]\n";
		return -1;
	}
	

	ofstream outfile;
	outfile.open(argv[1]);


	// Resulting constant offshell fit from SF
	const double pars_cstOff_SF[6] 			= {0.641326,2.3434,0.436479, -1.25864,1.01506,0.989354};
	const double pars_cstOff_no_SF[6] 		= {0.641326,2.3434,0.436479,0,	      1.01506,0.989354};
	// Resulting constant offshell fit from rho
	const double pars_cstOff_rho[6] 		= {0.670572,2.37343,0.413435,-0.866635,1.00928,0.990707};
	const double pars_cstOff_no_rho[6] 		= {0.670572,2.37343,0.413435,0,        1.00928,0.990707};


	for( double x = 0.2 ; x <= 0.96 ; x+=0.01 ){
		double Q2 = 14.*x;
		double he3 = 0; double h3 = 0;
		cout << "\tx = " << x << "\n";
		outfile << x << " " << Q2 << " ";

		// Resulting constant offshell fit from SF and rho
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.);
		outfile << he3 << " " << h3 << " ";
		calc_theo_rho(x,	Q2,	he3,	h3,	pars_cstOff_rho,	0.);
		outfile << he3 << " " << h3 << " ";

		// Resulting constant offshell fit from SF and rho with offshell FORCED to 0
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,		0.);
		outfile << he3 << " " << h3 << " ";
		calc_theo_rho(x,	Q2,	he3,	h3,	pars_cstOff_no_rho,		0.);
		outfile << he3 << " " << h3 << " ";
		


		// Now only allow for nucleons > 50 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.05);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.05);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 100 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.1);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.1);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 200 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.2);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.2);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 400 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.4);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.4);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 500 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.5);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.5);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 750 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		0.75);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	0.75);	// without offshell
		outfile << he3 << " " << h3 << " ";
		// Now only allow for nucleons > 1000 MeV/c in convolution with the spectral function
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_SF,		1.0);	// with offshell
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,		Q2,	he3,	h3,	pars_cstOff_no_SF,	1.0);	// without offshell
		outfile << he3 << " " << h3 << "\n";



		//calc_theo(x,Q2,he3,h3,pars,0.);			// Full EMC ratio
		//outfile << x << " " << Q2 << " " << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.);		// Only motion effects to EMC ratio
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.05);		// EMC contribution due to nucleons > 50 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.05);		// Only motion effects for > 50 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.100);		// EMC contribution due to nucleons > 100 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.100);	// Only motion effects for > 100 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.200);		// EMC contribution due to nucleons > 200 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.200);	// Only motion effects for > 200 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.400);		// EMC contribution due to nucleons > 400 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.400);	// Only motion effects for > 400 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.600);		// EMC contribution due to nucleons > 600 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.600);	// Only motion effects for > 600 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,0.800);		// EMC contribution due to nucleons > 800 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,0.800);	// Only motion effects for > 800 MeV/c
		//outfile << he3 << " " << h3 << " ";

		//calc_theo(x,Q2,he3,h3,pars,1.000);		// EMC contribution due to nucleons > 1000 MeV/c
		//outfile << he3 << " " << h3 << " ";
		//calc_theo(x,Q2,he3,h3,pars_noOff,1.000);	// Only motion effects for > 1000 MeV/c
		//outfile << he3 << " " << h3 << "\n";
		
	}
	outfile.close();
	
	return 0;
}

double offshell( double virt, double xB , double off_a  ){
	return 1 + (off_a)*virt*virt;
}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	//return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;	// -- simpler parameterization so that x=1 is just 1 parameter
}


void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut ){
	int A 		= 3;
	double N,Z;
	theo_he3 = 0; theo_h3 = 0;
	double counts = 0; double xalpha = 0;
	for(	itK_n = test_n.begin(), itK_p = test_p.begin(); 
			itK_n != test_n.end() && itK_p!= test_p.end(); itK_n++, itK_p++){

		double p_m = itK_n->first/1000.;
		if( p_m < mom_cut ) continue;
		for( itE_n = itK_n->second.begin(), itE_p = itK_p->second.begin(); 
				itE_n != itK_n->second.end() && itE_p != itK_p->second.end(); itE_n++, itE_p++){

			double E_m = itE_n->first/1000.;
			for( double theta = 0 ; theta <= M_PI ; theta += dTheta){

				// Convert from Kaptari E to initial E:
				double E_i = mHe3 - sqrt( pow( mHe3 - mP + E_m , 2 ) + p_m*p_m  );
				// Define y and nu from (E_i,p_i):
				double y = (E_i + p_m*cos(theta) )/mP;
				double nu = (E_i*E_i - p_m*p_m - mP*mP)/(mP*mP); //unitless

				// Grab spectral function values:
				double sp_n = itE_n->second;
				double sp_p = itE_p->second;

				double jacobian = sin(theta); // p^2 dp dE already inside spectral function definition
				double phi_int = 2*M_PI;
				//double flux_fact = 1. + (p_m*cos(theta))/sqrt(p_m*p_m + mP*mP);
				double flux_fact = 1. + (p_m*cos(theta))/E_i;
				double mass_fact = mHe3/(mP*3.);

				// Flag problematic kinematics for the delta conservation (?):
				if( y < x ) continue;
				if( y > A) continue;
				if( nu > 0. ) continue;

				xalpha += x/y;
				counts++;

				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are for the SF
				theo_he3 += jacobian * flux_fact * phi_int * mass_fact * dTheta 
					* ( Z*sp_p*offshell(nu,x/y,pars[3]) 
						+ N*sp_n*offshell(nu,x/y,pars[3])*f2nf2p(x/y,pars[0],pars[1],pars[2]) )
					* F2p->Eval(x/y, Q2 )
					* pars[4];
				Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n for only the SF
				theo_h3 += jacobian * flux_fact * phi_int * mass_fact * dTheta 
					* ( Z*sp_n*offshell(nu,x/y,pars[3]) 
						+ N*sp_p*offshell(nu,x/y,pars[3])*f2nf2p(x/y,pars[0],pars[1],pars[2]) )
					* F2p->Eval(x/y, Q2 )
					* pars[5];

			}// end loop over theta
		}// end loop over E miss
	}// end loop over p miss
	//cout << x << " " << xalpha/counts << "\n";

	theo_he3 *= (1./A) / F2d->Eval(x,Q2);
	theo_h3 *= (1./A) / F2d->Eval(x,Q2);
}


void calc_theo_rho( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut ){
	int A 		= 3;
	double N,Z;
	theo_he3 = 0; theo_h3 = 0;
	double baryon = 0, momentum = 0;
	double xalpha = 0;
	double counts = 0;
	for( int i = 0 ; i < rho_p.size() ; i++ ){

			double alpha = alphas.at(i);
			double nu = -nus.at(i)/(mP*mP);

			// Grab spectral function values:
			double sp_n = rho_n.at(i);
			double sp_p = rho_p.at(i);


			//cout << alpha << " " << nu*mP*mP << " " << sp_p << "\n";
			baryon += (2.*sp_p + 1.*sp_n);
			momentum += (2*sp_p + 1*sp_n)*alpha;

			// Flag problematic kinematics for the delta conservation:
			if( alpha < x ) continue;
			else if( alpha > A ) continue;
			if( nu > 0. ) continue;

			xalpha += x/alpha;
			counts++;

			Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are for the SF
			theo_he3 += ( Z*sp_p*offshell(nu,x/alpha,pars[3]) 
					+ N*sp_n*offshell(nu,x/alpha,pars[3])*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
				* F2p->Eval(x/alpha, Q2 )
				* pars[4];
			Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n for only the SF
			theo_h3 += ( Z*sp_n*offshell(nu,x/alpha,pars[3]) 
					+ N*sp_p*offshell(nu,x/alpha,pars[3])*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
				* F2p->Eval(x/alpha, Q2 )
				* pars[5];

	}// end loop over i
	//cout << x << " " << xalpha/counts << "\n";
	//cout << baryon/3. << " " << momentum/3. << "\n";
	theo_he3 *= (1./A) / F2d->Eval(x,Q2);
	theo_h3 *= (1./A) / F2d->Eval(x,Q2);
}
