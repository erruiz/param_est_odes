#include "odes.hh"

int main() {
	double p[4]={1.5,1.,3.,1.};
	lv_dop54 lv(p);
	lv.fixed_solve(10.,40,true);
}
