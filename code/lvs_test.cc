#include "odes.hh"

int main() {
	double p[4]={1.5,1.,3.,1.};
	lvs_dop54 lvs(p);
	lvs.fixed_solve(10.,500,true);
}
