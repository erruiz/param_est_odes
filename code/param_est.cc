#include <fstream>

#include "odes.hh"

// The file path for the reference solution.
std::string fname="../report/data/lv_test.out";

// The maximum number of epochs to run.
const int max_epochs=10000;

// The learning rate to use for gradient descent.
const double lr=1e-4;

// The true parameters of the system.
const double p[4]={1.5,1.,3.,1.};

int main() {
	printf("# The true parameters:\n");
	printf("# p = [%g %g %g %g]\n",p[0],p[1],p[2],p[3]);

	// Define the perturbed parameters. These will get updated after each epoch.
	double p_hat[4]={1.2,0.8,2.8,0.8};

	// Optimize the parameters by generating solutions with the perturbed
	// parameters and comparing these with the reference solution. The
	// parameters are updated according to the gradient descent algorithm.
	double loss;
	int epoch=0;
	do {
		lvs_dop54 lvs(p_hat);
		lvs.solve(10.,40,fname,lr,true);
		loss=lvs.loss;
		epoch++;
	} while (loss>1e-9 and epoch<max_epochs);
	
	printf("# The estimated parameters:\n");
	printf("# p_hat = [%g %g %g %g]\n",p_hat[0],p_hat[1],p_hat[2],p_hat[3]);
}
