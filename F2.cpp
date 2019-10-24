#include "F2.h"

using namespace std;

// Constructor
F2::F2( int choice ){
	if( choice == 0 ) makeProton();
	else if( choice == 1 ) makeDeuterium();
}
// Destructor
F2::~F2(){}

void F2::makeProton(void){
	m02 = 0.5063;
	mp2 = 34.75;
	mr2 = 0.03190;
	Q02 = 1.374;
	lam02 = 0.06527;
	ap1 = -0.11895;
	ap2 = -0.4783;
	ap3 = 1.353;
	bp1 = 1.0833;
	bp2 = 2.656;
	bp3 = 1.771;
	cp1 = 0.3638;
	cp2 = 0.1211;
	cp3 = 1.166;
	ar1 = 0.3425;
	ar2 = 1.0603;
	ar3 = 0.5164;
	br1 = -10.408;
	br2 = 14.857;
	br3 = 0.07739;
	cr1 = 1.3633;
	cr2 = 2.256;
	cr3 = 2.209;
}
void F2::makeDeuterium(void){
	m02 = 0.426;
	mp2 = 0.00007713;
	mr2 = 0.2293;
	Q02 = 2.65;
	lam02 = 0.06527;
	ap1 = -0.4287;
	ap2 = 0.2891;
	ap3 = 0.3931;
	bp1 = -27.212;
	bp2 = 30.687;
	bp3 = 0.04577;
	cp1 = 0.00073;
	cp2 = 0.9741;
	cp3 = 0.8722;
	ar1 = 0.2986;
	ar2 = 3.615;
	ar3 = 1.1455;
	br1 = 1.987;
	br2 = 7.150;
	br3 = 0.9350;
	cr1 = 1.0316;
	cr2 = 26.36;
	cr3 = 3.024;
}

double F2::Eval( double x, double Q2 ){
	return Q2 / ( Q2 + m02 ) * ( F2p(x,Q2) + F2r(x,Q2 ) );
}

double F2::makeT( double Q2){
	return log( log((Q2 + Q02)/lam02) / log(Q02/lam02) );
}

double F2::F2p( double x, double Q2 ){
	double t = makeT( Q2 );
	return cp(t) * pow( xp(x,Q2) , ap(t) ) * pow(1 - x , bp(t) );
}
double F2::ap( double t ){
	return ap1 + (ap1-ap2)*(1./(1.+pow(t,ap3)) - 1.);
}
double F2::bp( double t ){
	return bp1 + bp2*pow(t,bp3);
}
double F2::cp( double t ){
	return cp1 + (cp1-cp2)*(1./(1.+pow(t,cp3)) - 1.);
}
double F2::xp( double x, double Q2 ){
	double W2 = mP*mP + Q2*( 1./x - 1);
	return 1./(1. + (W2 - mP*mP)/(Q2 + mp2) );
}


double F2::F2r( double x, double Q2 ){
	double t = makeT( Q2 );
	return cr(t) * pow( xr(x,Q2) , ar(t) ) * pow(1 - x , br(t) );
}
double F2::ar( double t ){
	return ar1 + ar2*pow(t,ar3);
}
double F2::br( double t ){
	return br1 + br2*pow(t,br3);
}
double F2::cr( double t ){
	return cr1 + cr2*pow(t,cr3);
}
double F2::xr( double x, double Q2 ){
	double W2 = mP*mP + Q2*( 1./x - 1);
	return 1./(1. + (W2 - mP*mP)/(Q2 + mr2) );
}
