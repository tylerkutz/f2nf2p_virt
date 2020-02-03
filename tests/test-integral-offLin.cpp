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
double offshell( double virt, double xB , double off_a0 , double off_a1);
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c );

// Function to create theory point
void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars, double mom_cut );
void calc_sum_rules( double x, double Q2, double &baryon_he3, double &baryon_h3 );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [OutputTextFile]\n";
		return -1;
	}
	
	// He-3 and H-3 with offshell linear parameterization
	const double pars[7] = { 0.540931, 2.70716, 0.381848, -2.29884, 4.21216, 0.977637, 1.00757 };
	// H-3 with offshell linear parameterization
	//const double pars[7] = { 0.409517, 1.94394, 0.321044, 1.98499, 0.722875, 1, 0.999986 };
	// He-3 with offshell linear parameterization
	//const double pars[7] = { 0.257715, 1.86088, 0.220412, 3.67788, -2.2937, 1, 1 };

	ofstream outfile;
	outfile.open(argv[1]);

	for( double x = 0.2 ; x <= 0.95 ; x+=0.01 ){
		double Q2 = 14.*x;
		//double baryon_he3 = 0;
		//double baryon_h3 = 0;
		//calc_sum_rules(x,Q2,baryon_he3,baryon_h3);
		//cout << x << " " << baryon_he3/3.*100 << " " << baryon_h3/3.*100 << "\n";
		double he3_0 = 0, he3_200 = 0, he3_400 = 0, he3_600 = 0, he3_800 = 0;
		double h3_0 = 0,  h3_200 = 0,  h3_400 = 0,  h3_600 = 0, h3_800 = 0;
		calc_theo(x,Q2,he3_0,h3_0,pars,0.);
		//cout << x << " " << Q2 << " " << he3_0 << " " << h3_0 << "\n";
		outfile << x << " " << Q2 << " " << he3_0 << " " << h3_0 << "\n";
		//calc_theo(x,Q2,he3_200,h3_200,pars,0.05);
		//calc_theo(x,Q2,he3_400,h3_400,pars,0.075);
		//calc_theo(x,Q2,he3_600,h3_600,pars,0.100);
		//calc_theo(x,Q2,he3_800,h3_800,pars,0.200);
		//cout << x << " " << Q2 << " " << he3_0 << " " << h3_0 << " " << he3_200 << " " << h3_200 << " " << he3_400 << " " << h3_400 << " " << he3_600 << " " << h3_600 << " " << he3_800 << " " << h3_800 << "\n";
		//outfile << x << " " << Q2 << " " << he3_0 << " " << h3_0 << " " << he3_200 << " " << h3_200 << " " << he3_400 << " " << h3_400 << " " << he3_600 << " " << h3_600 << " " << he3_800 << " " << h3_800 << "\n";
	}
	outfile.close();
	
	return 0;
}

double offshell( double virt, double xB , double off_a0 , double off_a1){
	return 1 + (off_a0 + off_a1*xB)*virt*virt;
}
//double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c, double f2nf2p_d ){
//	return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
//}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	//return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;	// -- simpler parameterization so that x=1 is just 1 parameter
}

void calc_sum_rules( double x, double Q2, double &baryon_he3, double &baryon_h3 ){
	int A 		= 3;
	double N,Z;
	double he3 = 0;
	double h3 = 0;
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

				//double jacobian = (1./alpha) *sin(theta)
				//	* (mHe3 - mP + E_m)
				//	/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
				double jacobian = (1.)*sin(theta);
				double phi_int = 2*M_PI;

				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are
				baryon_he3 += jacobian * phi_int * dTheta 
					* ( Z*sp_p + N*sp_n );
				Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n
				baryon_h3 += jacobian * phi_int * dTheta 
					* ( Z*sp_n + N*sp_p );
			}// end loop over theta
		}// end loop over E miss
	}// end loop over p miss

}

void calc_theo( double x, double Q2, double &theo_he3, double &theo_h3 , const double *pars, double mom_cut ){
	int A 		= 3;
	double N,Z;
	double nu_min = 1e4; double nu_max = -1e4;
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
				if( nu > nu_max) nu_max = nu;
				if( nu < nu_min) nu_min = nu;

				double jacobian = (1./alpha) *sin(theta)
					* (mHe3 - mP + E_m)
					/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
				double phi_int = 2*M_PI;

				cout << sp_p*jacobian*phi_int*dTheta << " " << alpha << " " << nu << "\n";

				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are
				theo_he3 += jacobian * phi_int * dTheta 
					* ( Z*sp_p + N*sp_n*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* offshell(nu,x/alpha,pars[3],pars[4])
					* F2p->Eval(x/alpha, Q2 )
					* pars[5];	
				Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n
				theo_h3 += jacobian * phi_int * dTheta 
					* ( Z*sp_n + N*sp_p*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* offshell(nu,x/alpha,pars[3],pars[4])
					* F2p->Eval(x/alpha, Q2 )
					* pars[6];	

			}// end loop over theta
		}// end loop over E miss
	}// end loop over p miss
	theo_he3 *= (1./A) / F2d->Eval(x,Q2);
	theo_h3 *= (1./A) / F2d->Eval(x,Q2);
}


