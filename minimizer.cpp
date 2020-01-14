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

#include "Minuit2/FCNBase.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TVirtualFitter.h"
//#include "TFitterMinuit.h"

using namespace std;

// Data reader
void readData();
// Data storage
double data_x	 	[22]= {0.};
double data_he3 	[22]= {0.};
double data_he3_er 	[22]= {0.};
double data_h3 		[22]= {0.};
double data_h3_er 	[22]= {0.};

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
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c, double f2nf2p_d );

// Function to minimize
double Chi2( const double *pars );


// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./minimizer [OutputTextFile]\n";
		return -1;
	}

	// Starting parameters
	const double np_a = -1.21721713;
	const double np_b = 0.8622478;
	const double np_c = 0.82047886;
	const double np_d = 0.96399233;
	const double of_a = 0.;

	readData();

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","MIGRAD");
	min->SetMaxFunctionCalls(100000);
	//min->SetTolerance(0.5);
	min->SetTolerance(5);
	min->SetPrintLevel(1);

	ROOT::Math::Functor f(&Chi2,5);
	min->SetFunction(f);
	min->SetVariable(0,	"np_a",	np_a, 	0.1	);
	min->SetVariable(1,	"np_b",	np_b, 	0.1	);
	min->SetVariable(2,	"np_c",	np_c, 	0.1	);
	min->SetVariable(3,	"np_d",	np_d, 	0.1	);
	min->SetVariable(4,	"of_a",	of_a, 	0.1	);
	min->Minimize();

	//double covar[] = {0.};
	//min->GetCovMatrix(covar);

	// Print covar result:
	ofstream outfile;
	outfile.open(argv[1]);
	cerr 		<< "\n************ Fit covar results *************\n";
	outfile 	<< "\n************ Fit covar results *************\n";
	cerr 		<< "\t\t [ i*ndim + j ] \n";
	outfile 	<< "\t\t [ i*ndim + j ] \n";
	//for( auto i : covar ){
	//	auto val = i;
	//	cerr 	<< val << "\n";
	//	outfile << val << "\n";
	//}
	for( int i = 0 ; i < 5 ; i++ ){
		for( int j = 0 ; j < 5 ; j++){
			cerr << min->CovMatrix(i,j) << "\n";
			outfile << min->CovMatrix(i,j) << "\n";
		}
	}
	cerr 		<< "**********************************************" << endl;
	outfile 	<< "**********************************************" << endl;
	outfile.close();
	delete min;

	/*
	int second = sizeof(covar[0]);
	cout << second << std::endl;
	int manual = first/second;
	cout << manual << std::endl;
	for( int i = 0 ; i < size ; i ++ ){
		double val = covar[i];
		cout << val << "\n";
	}
	//for( int i = 0 ; i < size ; i++ ){
	//	cout << covar[i] << "\n";
	//}
	*/
	return 0;
}

double offshell( double virt, double xB , double off_a ){
	return 1 + off_a*virt*virt;
}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c, double f2nf2p_d ){
	return f2nf2p_a + f2nf2p_b*xB + f2nf2p_c*exp( f2nf2p_d*(1.-xB) );
}

void readData(){
	ifstream f;

	f.open("../exp_data/marathon_prelim_emc.txt");
	double x, he3, he3e, h3, h3e, he3h3, he3h3e;
	int i = 0;
	string dummyLine;
	getline(f, dummyLine);
	while(!f.eof()){
		f >> x;
		f >> he3;
		f >> he3e;
		f >> h3;
		f >> h3e;
		f >> he3h3;
		f >> he3h3e;
		data_x[i] = x;
		data_he3[i] = he3;	data_he3_er[i] = he3e;
		data_h3[i] = h3;	data_h3_er[i] = h3e;
		i++;
	}
	return;
}


double Chi2( const double *pars ){
	// par 0-3 	= f2n/f2p parameters
	// par 4 	= virt parameter
	
	// Want to return chi2 which is calculated given the parameters to be minimized
	double chi2 = 0.;
	cerr << "************ Fitting in progress *************\n";
	// Perform loop over exp data
	for( int i = 0 ; i < 22 ; i++ ){
		if( i%5==0) cerr << "\ton point " << i << "\n";
		double x 	= data_x[i];
		double Q2 	= 14.*x;
		int A 		= 3;
		double N,Z;

		// For each value in x, need to evaluate integral for both He-3, H-3:
		double theo_he3 = 0;
		double theo_h3	= 0;
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
						* ( Z*sp_p + N*sp_n*f2nf2p(x/alpha,pars[0],pars[1],pars[2],pars[3]) )
						* offshell(nu,x/alpha,pars[4])
						* F2p->Eval(x/alpha, Q2 );	
					Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n
					theo_h3 += jacobian * phi_int * dTheta 
						* ( Z*sp_n + N*sp_p*f2nf2p(x/alpha,pars[0],pars[1],pars[2],pars[3]) )
						* offshell(nu,x/alpha,pars[4])
						* F2p->Eval(x/alpha, Q2 );	

				}// end loop over theta
			}// end loop over E miss
		}// end loop over p miss
		theo_he3 *= (1./A) / F2d->Eval(x,Q2);
		theo_h3 *= (1./A) / F2d->Eval(x,Q2);

		// Calculate chi2:
		chi2 += 	pow(	(data_he3[i] - theo_he3) /data_he3_er[i]	, 2 );
		chi2 += 	pow(	(data_h3[i] - theo_h3)   /data_h3_er[i]		, 2 );

	}// end loop over data
	
	cerr << "------------Finished calculations!------------\n";
	cerr << "\tCurrent parameters:\n";
	cerr << "\t\tnp_a: " << pars[0] << "\n";
	cerr << "\t\tnp_b: " << pars[1] << "\n";
	cerr << "\t\tnp_c: " << pars[2] << "\n";
	cerr << "\t\tnp_d: " << pars[3] << "\n";
	cerr << "\t\tof_a: " << pars[4] << "\n";
	cerr << "\tCurrent chi2:\n";
	cerr << "\t\t" << chi2 << "\n";
	cerr << "**********************************************\n\n";

	return chi2;
}
