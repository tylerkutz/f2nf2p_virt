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

//using Eigen::Vector2d;
using namespace std;

// Initialize spectral function class
spec * kaptari_sf = new spec();
F2 * F2p = new F2(0);
F2 * F2d = new F2(1);

double offshell( double virt, double xB );
const double offshell_a = 1;
double f2nf2p( double xB );
const double f2nf2p_a = -1.21721713;
const double f2nf2p_b = 0.8622478;
const double f2nf2p_c = 0.82047886;
const double f2nf2p_d = 0.96399233;
const double f2p_a = 3.07;
const double f2p_b = 0.75;
const double f2p_c = 0.14;
const double f2p_d = -0.19;
const double f2p_e = -2.93;
const double f2p_f = -0.05;
const double f2p_g = 3.65;
//double f2pf2d( double xB , double Q2 );
//double f2p( double xB, double Q2 );

// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./minimizer [OutputTextFile]\n";
		return -1;
	}

	int A = 3;
	int Z = 1;
	int N = A-Z;
	double Q2 = 10;
	double dTheta = 0.05;
	
	TFile * outFile = new TFile("out.root","RECREATE");
	TTree * outTree = new TTree("kin","kin");
	double xB, q2, pmiss, Emiss, Estruck, Alpha, spec_p, p_F2, np_F2, Theta;
	outTree->Branch("xB",		&xB);
	outTree->Branch("Q2",		&q2);
	outTree->Branch("p_m",		&pmiss);
	outTree->Branch("E_m",		&Emiss);
	outTree->Branch("theta",	&Theta);
	outTree->Branch("E_s",		&Estruck);
	outTree->Branch("alpha",	&Alpha);
	outTree->Branch("sp_p",		&spec_p);
	outTree->Branch("F2p",		&p_F2);
	outTree->Branch("F2nF2p",	&np_F2);
	TH2D * hEmPm = new TH2D("hEmPm","hEmPm",1000,0,5,1000,0,5);
	TH1D * hCos  = new TH1D("hCos","hCos",500,-1,1);
	
	ofstream outfile;
	outfile.open(argv[1]);
	// Loop over xB
	for( double x = 0.0; x <= 1. ; x+=0.005 ){
		// Spec functions for n/p in He-3
		map<double,map<double,double>> test_n = kaptari_sf->getFullN();
		map<double,map<double,double>> test_p = kaptari_sf->getFullP();
		// 	iteraotrs to loop over 
		map<double,map<double,double>>::iterator itK_n;
		map<double,map<double,double>>::iterator itK_p;
		map<double,double>::iterator itE_n;
		map<double,double>::iterator itE_p;

		double integral = 0;
		for( itK_n = test_n.begin(), itK_p = test_p.begin(); itK_n != test_n.end() && itK_p!= test_p.end(); itK_n++, itK_p++){
			double p_m = itK_n->first/1000.;
			for( itE_n = itK_n->second.begin(), itE_p = itK_p->second.begin(); itE_n != itK_n->second.end() && itE_p != itK_p->second.end(); itE_n++, itE_p++){
				double E_m = itE_n->first/1000.;
				for( double theta = 0 ; theta <= M_PI ; theta += dTheta){
					
					// Convert from Kaptari E to initial E:
					double E_i = mHe3 - sqrt( pow( mHe3 - mP + E_m , 2 ) + p_m*p_m  );
					// Define alpha and nu from (E_i,p_i):
					double alpha = (E_i - p_m*cos(theta) )/mP;
					double nu = (E_i*E_i - p_m*p_m - mP*mP);

					// Grab spectral function values:
					double sp_n = itE_n->second;
					double sp_p = itE_p->second;
	

					// Flag problematic kinematics:
					if( x/alpha < 0 || x/alpha > 1 ){
						outTree->Fill();
						continue;
					}
					
					double jacobian = (1./alpha) *sin(theta)
								* (mHe3 - mP + E_m)
								/ (mP * sqrt(pow(mHe3-mP+E_m,2) + p_m*p_m) );
					double phi_int = 2*M_PI;
					integral += jacobian * phi_int * dTheta 
							* ( Z*sp_p + N*sp_n*f2nf2p(x/alpha) )
							* offshell(nu,x/alpha)
							* F2p->Eval(x/alpha, Q2 );	

				}
			}
		}
		outfile << x << " " << integral << " " <<  (1./A)*integral / F2d->Eval(x,Q2) << "\n";
	}
	outfile.close();
	outFile->cd();
	outTree->Write();
	outFile->Close();

	return 0;
}

double offshell( double virt, double xB ){
	return 1;
}
double f2nf2p( double xB ){
	return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
}

//double f2p( double xB, double Q2 ){
//	return (f2p_a*pow(x,f2p_b) + f2p_c*pow(x,f2p_d)*(1.+f2p_e*sqrt(x))*(log(Q2)+f2p_f*pow(log(Q2),2))) * pow(1.-x,f2p_g);
//}






/*
double f2pf2d( double xB , double Q2 ){
	// Using Whitlow parameterization for the moment
	Vector2d C1( 1.417, 0.948); 
	Vector2d C2(-0.108,-0.115); 
	Vector2d C3( 1.486, 1.861); 
	Vector2d C4(-5.979,-4.733); 
	Vector2d C5( 3.524, 2.348); 
	Vector2d C6(-0.011,-0.065); 
	Vector2d C7(-0.619,-0.224); 
	Vector2d C8( 1.385, 1.085); 
	Vector2d C9( 0.270, 0.213); 
	Vector2d C10(-2.179,-1.687); 
	Vector2d C11( 4.722, 3.409); 
	Vector2d C12(-4.363,-3.255); 

	Vector2d beta(1., pow( 1. - exp( -min(20.,7.7*(1./xB + mP*mP/Q2 - 1. ) ) ) , -1 ) );
	double A = 1.22*exp(3.2*xB);

	Vector2d lambda1, lambda2, F2_thr;
	lambda2 = C6 + C7*xB + C8*xB*xB;
	lambda1 = C9 + C10*xB + C11*xB*xB + C12*xB*xB*xB;
	F2_thr = C1*pow(1.-xB,3) + C2*pow(1.-xB,4) + C3*pow(1.-xB,5) + C4*pow(1.-xB,6) + C5*pow(1.-xB,7);

	if( Q2 > A ) lambda2(0.,0.);
	double mult0 = 1. + lambda1[0]*log(Q2/A) + lambda2[0]*pow( log(Q2/A) , 2 );
	double mult1 = 1. + lambda1[1]*log(Q2/A) + lambda2[1]*pow( log(Q2/A) , 2 );

	return (beta[0]*F2_thr[0]*mult0)/(beta[1]*F2_thr[1]*mult1)/2.;
}
*/
