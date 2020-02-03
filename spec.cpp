#include "spec.h"

using namespace std;

// Constructor
spec::spec(){
	fill_arrays();
	fill_rho();
}
// Destructor
spec::~spec(){}

// Private member function
void spec::fill_rho(){
	std::ifstream f_corr("../../include/jacksonLC_A3.dat");
	double alpha, dalpha, nu, dnu, rho_p, rho_n;
	if (f_corr.is_open() == true){
		while (!f_corr.eof()) { 
			f_corr >> alpha >> dalpha >> nu >> dnu >> rho_p >> rho_n;
			full_rho_proton[alpha][nu] = rho_p;
			full_rho_neutron[alpha][nu] = rho_n;
		}
	}
	else{ cout << "\tFailed to find rho...\n"; }
}

// Private member function
void spec::fill_arrays(){
	std::ifstream f_corr("../../include/kaptarisf_3par.dat");

	char line[256]; 
	for( int i = 0 ; i < 1 ; i ++) f_corr.getline(line, 256);  // skip the first line
	double k, E, SF_p, SF_n, dk, dE;
	int ctr_E = 0;
	if (f_corr.is_open() == true){ 
		while (!f_corr.eof()) { 
			f_corr >> k >> E >> SF_p >> SF_n >> dk >> dE;
			SF_p = SF_p / (4.*M_PI);
			SF_n = SF_n / (4.*M_PI);
			ks.push_back(k);
			Es.push_back(E);
			kSteps.push_back(dk);
			ESteps.push_back(dE);
			full_SF_proton[k][E] = SF_p;
			full_SF_neutron[k][E] = SF_n;
			if( ctr_E == 0 ){ // 2-body breakup channel in He-3
				Em_2body = E;
				n2_proton[k] = SF_p;
				n2_neutron[k] = SF_n; // should be 0
			}
			else{ // continuum breakup channel in He-3
				if( ctr_E == 1 ) Em_thres = E;
				contin_proton[k][E] = SF_p;
				contin_neutron[k][E] = SF_n;
			}
			ctr_E++;
		} 
	} 
	else{ cout << "\tFailed to find wavefunction...\n"; }  
}

// Public getter member functions
std::map<double,std::map<double,double>> spec::getFullN(void){
	return full_SF_neutron;
}
std::map<double,std::map<double,double>> spec::getFullP(void){
	return full_SF_proton;
}
std::map<double,std::map<double,double>> spec::getFullRhoN(void){
	return full_rho_neutron;
}
std::map<double,std::map<double,double>> spec::getFullRhoP(void){
	return full_rho_proton;
}
std::map<double,std::map<double,double>> spec::getContinuumN(void){
	return contin_neutron;
}
std::map<double,std::map<double,double>> spec::getContinuumP(void){
	return contin_proton;
}
std::map<double,double> spec::get2bodyN(void){
	return n2_neutron;
}
std::map<double,double> spec::get2bodyP(void){
	return n2_proton;
}
double spec::getEm_2body(void){
	return Em_2body;
}
double spec::getEm_thres(void){
	return Em_thres;
}
std::vector<double> spec::getKs(void){
	return ks;
}
std::vector<double> spec::getEs(void){
	return Es;
}
std::vector<double> spec::getKsteps(void){
	return kSteps;
}
std::vector<double> spec::getEsteps(void){
	return ESteps;
}
