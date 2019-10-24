#ifndef __F2_H__
#define __F2_H__

#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <cmath>
#include "constants.h"

class F2
{
	public:
		F2(int);
		~F2();
		void makeProton(void);
		void makeDeuterium(void);
		double Eval( double, double );
		
	private:
		double m02;
		double mp2;
		double mr2;
		double Q02;
		double lam02;
		double ap1;
		double ap2;
		double ap3;
		double bp1;
		double bp2;
		double bp3;
		double cp1;
		double cp2;
		double cp3;
		double ar1;
		double ar2;
		double ar3;
		double br1;
		double br2;
		double br3;
		double cr1;
		double cr2;
		double cr3;
		
		double makeT( double );

		double F2p( double, double );
		double F2r( double, double );

		double ap( double );
		double bp( double );
		double cp( double );
		double ar( double );
		double br( double );
		double cr( double );

		double xp( double x, double Q2 );
		double xr( double x, double Q2 );

};

#endif
