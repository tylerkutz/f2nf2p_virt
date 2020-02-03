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
const double dTheta = 0.01;

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
void calc_theo( double x, double Q2, double &baryon, double &momentum , const double *pars );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [OutputTextFile]\n";
		return -1;
	}
	
	// He-3 + H-3 fit with jacobian fix
	const double pars[6] 		= {0.568697,2.20877,0.415621,-0.118305,1.00649,0.995528};
	const double pars_noOff[6] 	= {0.568697,2.20877,0.415621,0.0,1.00649,0.995528};

	ofstream outfile;
	outfile.open(argv[1]);

	for( double x = 0.2 ; x <= 0.96 ; x+=0.01 ){
		double Q2 = 14.*x;
		double baryon = 0, momentum = 0;
		cout << "\tx = " << x << "\n";
		
		calc_theo( x , Q2 , baryon , momentum , pars );
		cout << x << " " << Q2 << " " << baryon << " " << momentum << "\n";
		outfile << x << " " << Q2 << " " << baryon << " " << momentum << "\n";

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


void calc_theo( double x, double Q2, double &baryon, double &momentum , const double *pars ){
	int A 		= 3;
	double N,Z;
	for(	itK_n = test_n.begin(), itK_p = test_p.begin(); 
			itK_n != test_n.end() && itK_p!= test_p.end(); itK_n++, itK_p++){

		double p_m = itK_n->first/1000.;
		for( itE_n = itK_n->second.begin(), itE_p = itK_p->second.begin(); 
				itE_n != itK_n->second.end() && itE_p != itK_p->second.end(); itE_n++, itE_p++){

			double E_m = itE_n->first/1000.;
			for( double theta = 0 ; theta <= M_PI ; theta += dTheta){

				// Convert from Kaptari E to initial E:
				double E_i = mHe3 - sqrt( pow( mHe3 - mP + E_m , 2 ) + p_m*p_m  );
				// Define alpha and nu from (E_i,p_i):
				double alpha = (E_i - p_m*cos(theta) )/mP;
				double nu = (E_i*E_i - p_m*p_m - mP*mP)/(mP*mP); //unitless

				// Grab spectral function values:
				double sp_n = itE_n->second;
				double sp_p = itE_p->second;

				// Flag problematic kinematics:
				//if( x/alpha < 0 || x/alpha > 1 ){
				//	continue;
				//}
				//if( alpha < 0 ) continue;
				//if( x / alpha > 1 ) continue;

				//double jacobian = (1./alpha) *sin(theta)
				//	* (mHe3 - mP + E_m)
				//	/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
				double jacobian = sin(theta); // p^2 dp dE already inside spectral function definition
				double phi_int = 2*M_PI;
				double flux_fact = 1. + (p_m*cos(theta))/sqrt(p_m*p_m + mP*mP);


				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are for the SF
				baryon += jacobian * flux_fact * phi_int * dTheta
					* ( Z*sp_p + N*sp_n );
				momentum += jacobian * flux_fact * phi_int * dTheta
					* ( Z*sp_p + N*sp_n )
					* alpha;

			}// end loop over theta
		}// end loop over E miss
	}// end loop over p miss
	baryon *= (1./3);
	momentum *= (1./3);
}


