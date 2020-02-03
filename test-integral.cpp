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

// Functions with parameters to be minimized
double offshell( double virt, double xB , double off_a );
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c );

// Function to create theory point
void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [OutputTextFile]\n";
		return -1;
	}
	
	// He-3 + H-3 with offshell and np starting from init
	//const double pars[5] = {-1.49198,1.35858,0.694808,1.29278,-2.56522};
	// He-3 only with offshell and np starting from init
	//const double pars[5] = {-1.28977,0.791046,0.777979,0.905046,1.50513};
	// H-3 ony with offshell and np starting from init
	//const double pars[5] = {-1.07099,0.481792,0.876975,0.704521,2.95809};
	// He-3 only with fixed offshell starting from init				-- DONE
	//const double pars[5] = {-2.97006,1.81705,1.58944,0.862649,0};
	// H-3 ony with fixed offshell starting from init
	//const double pars[5] = {-1.29765,0.933887,0.78149,1.01652,0};
	// He-3 and H-3 with fixed offshell and startiing from init
	//const double pars[5] = {-1.32074,0.970307,0.770459,1.03452,0};
	
	// He-3 + H-3 with offshell and np starting with new param and norm err from init
	//const double pars[6] = {0.439694,2.6614,0.378713,0.415886,0.970452,1.0135};
	// H-3 with offshell and np starting with new param and norm err from init
	//const double pars[6] = {0.394644,1.87921,0.315853,2.59768,1,0.999951};
	// He-3 with offshell and np starting with new param and norm err from init
	//const double pars[6] = {0.34658,2.3002,0.246982,1.82858,1.00002,1};
	
	// He-3 + H-3 fit with jacobian fix
	//const double pars[6] 		= {0.568697,2.20877,0.415621,-0.118305,1.00649,0.995528};
	//const double pars_noOff[6] 	= {0.568697,2.20877,0.415621,0.0,1.00649,0.995528};
	
	// He-3 + H-3 fit with new framework fix
	const double pars[6]		= {0.641236,2.34294,0.436482,-1.25748,1.015,0.989314};
	const double pars_noOff[6]	= {0.641236,2.34294,0.436482,0,1.015,0.989314};

	ofstream outfile;
	outfile.open(argv[1]);

	for( double x = 0.2 ; x <= 0.96 ; x+=0.2 ){
		double Q2 = 14.*x;
		double he3 = 0; double h3 = 0;
		cout << "\tx = " << x << "\n";

		calc_theo(x,Q2,he3,h3,pars,0.);			// Full EMC ratio
		outfile << x << " " << Q2 << " " << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.);		// Only motion effects to EMC ratio
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.05);		// EMC contribution due to nucleons > 50 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.05);		// Only motion effects for > 50 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.100);		// EMC contribution due to nucleons > 100 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.100);	// Only motion effects for > 100 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.200);		// EMC contribution due to nucleons > 200 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.200);	// Only motion effects for > 200 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.400);		// EMC contribution due to nucleons > 400 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.400);	// Only motion effects for > 400 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.600);		// EMC contribution due to nucleons > 600 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.600);	// Only motion effects for > 600 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,0.800);		// EMC contribution due to nucleons > 800 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,0.800);	// Only motion effects for > 800 MeV/c
		outfile << he3 << " " << h3 << " ";

		calc_theo(x,Q2,he3,h3,pars,1.000);		// EMC contribution due to nucleons > 1000 MeV/c
		outfile << he3 << " " << h3 << " ";
		calc_theo(x,Q2,he3,h3,pars_noOff,1.000);	// Only motion effects for > 1000 MeV/c
		outfile << he3 << " " << h3 << "\n";
		
	}
	outfile.close();
	
	return 0;
}

double offshell( double virt, double xB , double off_a ){
	return 1 + off_a*virt*virt;
}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	//return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;	// -- simpler parameterization so that x=1 is just 1 parameter
}


void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars , double mom_cut ){
	int A 		= 3;
	double N,Z;
	theo_he3 = 0; theo_h3 = 0;
	double meanxalpha = 0;
	double count = 0;
	double baryon = 0, moment = 0;
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

				meanxalpha += x/y;
				count++;

				//double jacobian = (1./alpha) *sin(theta)
				//	* (mHe3 - mP + E_m)
				//	/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
				double jacobian = sin(theta); // p^2 dp dE already inside spectral function definition
				double phi_int = 2*M_PI;
				//double flux_fact = 1. + (p_m*cos(theta))/sqrt(p_m*p_m + mP*mP);
				double flux_fact = 1. + (p_m*cos(theta))/E_i;
				double mass_fact = mHe3/(mP*3.);

				// Baryon and momentum sum rule are from y=(0,A), so no limit on y has to be < x
				baryon += 1./A*(Z*sp_p + N*sp_n)*mass_fact*flux_fact*phi_int*jacobian*dTheta		;
				moment += 1./A*(Z*sp_p + N*sp_n)*mass_fact*flux_fact*phi_int*jacobian*y*dTheta		;
				// Flag problematic kinematics for the delta conservation (?):
				if( y < x ) continue;
				if( y > A) continue;
				if( nu > 0. ) continue;



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
	theo_he3 *= (1./A) / F2d->Eval(x,Q2);
	theo_h3 *= (1./A) / F2d->Eval(x,Q2);
	cout << baryon << " " << moment << "\n";

	//cout << x << " " << meanxalpha/count << "\n";
}


