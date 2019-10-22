#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "spec.h"
#include "constants.h"

using namespace std;

// Initialize spectral function class
spec * kaptari_sf = new spec();

double offshell( double virt, double xB );
const double offshell_a = 1;
double f2nf2p( double xB );
const double f2nf2p_a = -1.21721713;
const double f2nf2p_b = 0.8622478;
const double f2nf2p_c = 0.82047886;
const double f2nf2p_d = 0.96399233;
double f2pf2d( double xB , double Q2 );

// Main
int main(int argc, char ** argv){

	if (argc<3){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./gen_eep [nEvents] [outputFile]\n";
		return -1;
	}

	int A = 3;
	int Z = 2;
	int N = A-Z;
	
	// Loop over xB
	for( double x = 0.2; x < 0.95 ; x+=0.01 ){
	
		// Spec functions for n/p in He-3
		std::map<double,std::map<double,double>> test_n = kaptari_sf->getFullN();
		std::map<double,std::map<double,double>> test_p = kaptari_sf->getFullP();
		// 	iteraotrs to loop over 
		std::map<double,std::map<double,double>>::iterator itK_n;
		std::map<double,std::map<double,double>>::iterator itK_p;
		std::map<double,double>::iterator itE_n;
		std::map<double,double>::iterator itE_p;
		double integral_p = 0;
		double integral_n = 0;
		for( itK_n = test_n.begin(), itK_p = test_p.begin(); itK_n != test_n.end() && itK_p!= test_p.end(); itK_n++, itK_p++){
			double k = itK_n->first/1000.;
			for( itE_n = itK_n->second.begin(), itE_p = itK_p->second.begin(); itE_n != itK_n->second.end() && itE_p != itK_p->second.end(); itE_n++, itE_p++){
				double E = itE_n->first/1000.;
				double v = (k*k + E*E - mN*mN);
				double sp_n = itE_n->second;
				double sp_p = itE_p->second;

				integral_p += Z * sp_p * offshell(v,x);
				integral_n += N * sp_n * offshell(v,x) * f2nf2p(x);
			}

		}
		
		cout << (integral_p + integral_n) * f2pf2d(x, 14.*x ) * (2./A) << "\n";

	}

	return 0;
}

double offshell( double virt, double xB ){
	return 1 ;//+ offshell_a*virt*virt;
}
double f2nf2p( double xB ){
	return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
}
double f2pf2d( double xB , double Q2 ){
	double C1  [2] = { 1.417, 0.948}; 
	double C2  [2] = {-0.108,-0.115}; 
	double C3  [2] = { 1.486, 1.861}; 
	double C4  [2] = {-5.979,-4.733}; 
	double C5  [2] = { 3.524, 2.348}; 
	double C6  [2] = {-0.011,-0.065}; 
	double C7  [2] = {-0.619,-0.224}; 
	double C8  [2] = { 1.385, 1.085}; 
	double C9  [2] = { 0.270, 0.213}; 
	double C10 [2] = {-2.179,-1.687}; 
	double C11 [2] = { 4.722, 3.409}; 
	double C12 [2] = {-4.363,-3.255}; 
	double beta [2] = { 1.	, pow( 1. - exp( -min(20.,7.7*(1./xB + mP*mP/Q2 - 1. ) ) ) , -1 ) };
	double A = 1.22*exp(3.2*xB);
	double lambda2[2] = C6 + C7*xB + C8*xB*xB;

	/*double lambda2[2] = C6 + C7*x + C8*x*x

	lambda1 = C9 + C10*x + C11*x*x + C12*x*x*x

	F2_thr = C1*np.power(1.-x,3) + C2*np.power(1.-x,4) + C3*np.power(1.-x,5) + C4*np.power(1.-x,6) + C5*np.power(1.-x,7)
	
	if( Q2 > A ): lambda2 = [0.,0.]
	mult0 = 1. + lambda1[0]*np.log(Q2/A) + lambda2[0]*np.power( np.log(Q2/A) , 2 )
	mult1 = 1. + lambda1[1]*np.log(Q2/A) + lambda2[1]*np.power( np.log(Q2/A) , 2 )
	return (beta[0]*F2_thr[0]*mult0)/(beta[1]*F2_thr[1]*mult1)/2.*/
	return 1.;

}

