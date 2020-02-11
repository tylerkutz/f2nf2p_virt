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
#include "TMatrixD.h"

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

vector<double> alphas = kaptari_sf->getAlphas();
vector<double> nus = kaptari_sf->getNus();
vector<double> rho_p = kaptari_sf->getRhoProtons();
vector<double> rho_n = kaptari_sf->getRhoNeutrons();

// Functions with parameters to be minimized
double offshell( double virt, double xB , double off_a );
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c );

// Function to minimize
double Chi2_Rho( const double *pars );

int opt = -1;
// Main
int main(int argc, char ** argv){

	if( argc < 3 ){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./minimizer [Fit Opt] [OutputTextFile]\n";
		cerr << "\t[Fit Opt]:\n\t\t0 = He-3 + H-3\n";
		cerr << "\t\t1 = H-3\n";
		cerr << "\t\t2 = He-3\n\n";
		return -1;
	}

	
	// New parameters from new param of F2n/F2p
	// 	OG
	const double np_a 	= 0.568697;
	const double np_b 	= 2.20877;
	const double np_c 	= 0.41562;
	const double of_p 	= 0.;
	const double of_n 	= 0.;
	const double Nhe3 	= 1.;
	const double Nh3 	= 1.;

	readData();

	ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2","MIGRAD");
	min->SetMaxFunctionCalls(100000);
	min->SetTolerance(5);
	min->SetPrintLevel(1);

	ROOT::Math::Functor f(&Chi2_Rho,7);
	min->SetFunction(f);
	min->SetVariable(0,	"np_a",		np_a, 		0.1	);
	min->SetVariable(1,	"np_b",		np_b, 		0.1	);
	min->SetVariable(2,	"np_c",		np_c, 		0.1	);
	min->SetVariable(3,	"of_p",		of_p, 		0.1	);
	min->SetVariable(4,	"of_n",		of_n, 		0.1	);
	min->SetVariable(5,	"N_he3",	Nhe3,   	0.1	);
	min->SetVariable(6,	"N_h3",		Nh3,    	0.1	);
	opt = atoi(argv[1]);
	if( opt == 1 ){ cout << "Fixing He-3 norm, i.e. only doing H-3 fit\n"; min->FixVariable(5); }
	if( opt == 2 ){ cout << "Fixing H-3 norm, i.e. only doing He-3 fit\n"; min->FixVariable(6); }
	
	min->Minimize();

	cerr << "Printing covariance matrix:\n";
	TMatrixD covMatrix(min->NDim(), min->NDim() ); 
	min->GetCovMatrix( covMatrix.GetMatrixArray() ); 
	covMatrix.Print();
	cerr << "\n";

	// Print covar result:
	ofstream outfile;
	outfile.open(argv[2]);
	cerr 		<< "\n************ Fit covar results *************\n";
	outfile 	<< "\n************ Fit covar results *************\n";
	cerr 		<< "\t\t [ i*ndim + j ] \n";
	outfile 	<< "\t\t [ i*ndim + j ] \n";
	for( int i = 0 ; i < 7 ; i++ ){
		for( int j = 0 ; j < 7 ; j++){
			cerr << min->CovMatrix(i,j) << " " << min->Correlation(i,j) << " " << i << " " << j << "\n";
			outfile << min->CovMatrix(i,j) << " " << min->Correlation(i,j) << " " << i << " " << j << "\n";
		}
	}
	cerr 		<< "**********************************************" << endl;
	outfile 	<< "**********************************************" << endl;
	outfile.close();
	delete min;

	return 0;
}

double offshell( double virt, double xB , double off_a ){
	return 1 + (off_a)*virt*virt;
}
double f2nf2p( double xB , double f2nf2p_a, double f2nf2p_b, double f2nf2p_c ){
	return f2nf2p_a * pow( 1.-xB , f2nf2p_b ) + f2nf2p_c;
}

void readData(){
	ifstream f;

	f.open("../../exp_data/marathon_prelim_emc.txt");
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


double Chi2_Rho( const double *pars ){
	// par 0-2 	= f2n/f2p parameters
	// par 3-4 	= virt parameter
	// par 5-6	= norm uncertainty
	
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
		for( int i = 0 ; i < rho_p.size() ; i++ ){

				double alpha = alphas.at(i);
				double nu = -nus.at(i)/(mP*mP);

				// Grab spectral function values:
				double sp_n = rho_n.at(i);
				double sp_p = rho_p.at(i);

				// Flag problematic kinematics for the delta conservation:
				if( alpha < x ) continue;
				else if( alpha > A ) continue;
				if( nu > 0. ) continue;

				Z = 2; N = A-Z; // He3 - as it's He3 SF, keep p and n as they are for the SF
				theo_he3 += ( Z*sp_p*offshell(nu,x/alpha,pars[3]) 
						+ N*sp_n*offshell(nu,x/alpha,pars[4])*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* F2p->Eval(x/alpha, Q2 );
				Z = 1; N = A-Z; // H3 - as it's He3 SF, swap p and n for only the SF
				theo_h3 += ( Z*sp_n*offshell(nu,x/alpha,pars[3]) 
						+ N*sp_p*offshell(nu,x/alpha,pars[4])*f2nf2p(x/alpha,pars[0],pars[1],pars[2]) )
					* F2p->Eval(x/alpha, Q2 );

		}// end loop over i
		theo_he3 *= (1./A) / F2d->Eval(x,Q2);
		theo_h3 *= (1./A) / F2d->Eval(x,Q2);

		// Calculate chi2:
		if( opt == 0 || opt == 2 ){
			chi2 += 	pow(	(data_he3[i] - pars[5]*theo_he3) /data_he3_er[i]	, 2 );
			chi2 +=		pow(	(pars[5] - 1.)/0.05	, 2 );	// He-3 normalization 
		}
		if( opt == 0 || opt == 1 ){
			chi2 += 	pow(	(data_h3[i]  - pars[6]*theo_h3)  /data_h3_er[i]	, 2 );
			chi2 += 	pow(	(pars[6] - 1.)/0.05	, 2 );  // H-3 normalization
		}


	}// end loop over data
	
	cerr << "------------Finished calculations!------------\n";
	cerr << "\tCurrent parameters:\n";
	cerr << "\t\tnp_a: " 	<< pars[0] << "\n";
	cerr << "\t\tnp_b: " 	<< pars[1] << "\n";
	cerr << "\t\tnp_c: " 	<< pars[2] << "\n";
	cerr << "\t\tof_p: " 	<< pars[3] << "\n";
	cerr << "\t\tof_n: " 	<< pars[4] << "\n";
	cerr << "\t\tN_he3: "	<< pars[5] << "\n";
	cerr << "\t\tN_h3: " 	<< pars[6] << "\n";
	cerr << "\tCurrent chi2:\n";
	cerr << "\t\t" << chi2 << "\n";
	cerr << "**********************************************\n\n";

	return chi2;
}
