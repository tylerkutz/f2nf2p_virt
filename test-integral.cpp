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
void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars);


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
	const double pars[6] = {0.34658,2.3002,0.246982,1.82858,1.00002,1};

	ofstream outfile;
	outfile.open(argv[1]);
	
	for( double x = 0.2 ; x < 0.9 ; x+=0.01 ){
		double Q2 = 14.*x;
		double he3=0;
		double h3=0;
		calc_theo(x,Q2,he3,h3,pars);
		
		outfile << x << " " << Q2 << " " << he3 << " " << h3 << "\n";
	}
	outfile.close();
	
	return 0;
}

double offshell( double virt, double xB , double off_a ){
	return 1 + off_a*virt*virt;
}
//double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c, double f2nf2p_d ){
//	return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
//}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	//return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;	// -- simpler parameterization so that x=1 is just 1 parameter
}


void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars){
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
				if( x/alpha < 0 || x/alpha > 1 ){
					continue;
				}

				double jacobian = (1./alpha) *sin(theta)
					* (mHe3 - mP + E_m)
					/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
				double phi_int = 2*M_PI;

				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are
				theo_he3 += jacobian * phi_int * dTheta 
					* ( Z*sp_p + N*sp_n*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* offshell(nu,x/alpha,pars[3])
					* F2p->Eval(x/alpha, Q2 )
					* pars[4];	
				Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n
				theo_h3 += jacobian * phi_int * dTheta 
					* ( Z*sp_n + N*sp_p*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* offshell(nu,x/alpha,pars[3])
					* F2p->Eval(x/alpha, Q2 )
					* pars[5];	

			}// end loop over theta
		}// end loop over E miss
	}// end loop over p miss
	theo_he3 *= (1./A) / F2d->Eval(x,Q2);
	theo_h3 *= (1./A) / F2d->Eval(x,Q2);
}


