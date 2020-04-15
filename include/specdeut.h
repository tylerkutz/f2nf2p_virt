#ifndef __SPECDEUT_H__
#define __SPECDEUT_H__

#include "TSpline.h"
#include "constants.h"
#include "TH1.h"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <map>


class specdeut
{
	public:
		specdeut();
		~specdeut();
		// Getter function for momentum distribution
		double getDensity(double k);

		std::map<double,double> getMomDist(void);

	private:
		// Private storage for the momentum distribution
		std::map<double,double> full_mom_dist;

		double density[201];
		double kPts[201];
		TSpline3 * splDensity;
		void fill_arrays();

};

#endif
