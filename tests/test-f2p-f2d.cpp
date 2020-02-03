#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "spec.h"
#include "constants.h"
#include "F2.h"

using namespace std;

// Initialize spectral function class
spec * kaptari_sf = new spec();
F2 * F2p = new F2(0);
F2 * F2d = new F2(1);

// Main
int main(int argc, char ** argv){

	if (argc<2){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./minimizer [compareF2p.txt] [compareF2d.txt]\n";
		return -1;
	}

	int A = 3;
	int Z = 2;
	int N = A-Z;
	
	ifstream infile; infile.open("../exp_data/f2p.txt");
	ofstream outfile; outfile.open(argv[1]);
	double xBin, x, Q2, F2, stat, pid, model, mis, rad;
	if (infile.is_open() == true){ 
		while (!infile.eof()) { 
			infile >> xBin >> x >> Q2 >> F2 >> stat >> pid >> model >> mis >> rad;
			outfile << F2*(1-(stat+pid+model+mis+rad)/100.) << " " << F2 << " " << F2*(1+(stat+pid+model+mis+rad)/100.) << " " << F2p->Eval(x,Q2) << " " << fabs(F2p->Eval(x,Q2) - F2) << " " << fabs(F2p->Eval(x,Q2) - F2)/F2 << "\n";
		}
	}

	ifstream infile2; infile2.open("../exp_data/f2d.txt");
	ofstream outfile2; outfile2.open(argv[2]);
	if (infile2.is_open() == true){ 
		while (!infile2.eof()) { 
			infile2 >> xBin >> x >> Q2 >> F2 >> stat >> pid >> model >> mis >> rad;
			outfile2 << F2*(1-(stat+pid+model+mis+rad)/100.) << " " << F2 << " " << F2*(1+(stat+pid+model+mis+rad)/100.) << " " << F2d->Eval(x,Q2) << " " << fabs(F2d->Eval(x,Q2) - F2) << " " << fabs(F2d->Eval(x,Q2) - F2)/F2 << "\n";
		}
	}
	

	outfile.close();
	outfile2.close();

	return 0;
}

