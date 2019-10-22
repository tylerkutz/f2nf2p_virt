#ifndef __SPEC_H__
#define __SPEC_H__

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>

class spec
{
	public:
		spec();
		~spec();
		// Getter functions
		std::vector<double> getKs(void);
		std::vector<double> getEs(void);
		std::vector<double> getKsteps(void);
		std::vector<double> getEsteps(void);
		std::map<double,std::map<double,double>> getFullN(void);
		std::map<double,std::map<double,double>> getFullP(void);
		std::map<double,std::map<double,double>> getContinuumN(void);
		std::map<double,std::map<double,double>> getContinuumP(void);
		std::map<double,double> get2bodyN(void);
		std::map<double,double> get2bodyP(void);
		double getEm_2body(void);
		double getEm_thres(void);
	private:
		// Private storage
		std::map<double,std::map<double,double>> full_SF_proton;
		std::map<double,std::map<double,double>> full_SF_neutron;
		std::map<double,std::map<double,double>> contin_proton;
		std::map<double,std::map<double,double>> contin_neutron;
		double Em_2body;
		double Em_thres;
		std::map<double,double> n2_proton;
		std::map<double,double> n2_neutron;
		std::vector<double> ks;
		std::vector<double> Es;
		std::vector<double> kSteps;
		std::vector<double> ESteps;
		void fill_arrays();

};

#endif
