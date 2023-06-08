#include "dopri54.hh"

/** The class constructor initializes the constants and dynamically allocates
 * memory for the solution and all intermediate calculations. */
dopri54::dopri54(int dof_,double *p_,double t0): dof(dof_), t(t0), fcount(0),
	dq(new double[dof]), k1(new double[dof]), k2(new double[dof]), 
	k3(new double[dof]), k4(new double[dof]), k5(new double[dof]),
	k6(new double[dof]), k7(new double[dof]), q(new double[dof]), 
	q_hat(new double[dof]), p(p_) {}

/** The class destructor frees the dynamically allocated memory. */
dopri54::~dopri54() {
	delete [] q_hat;
	delete [] q;
	delete [] k7;
	delete [] k6;
	delete [] k5;
	delete [] k4;
	delete [] k3;
	delete [] k2;
	delete [] k1;
	delete [] dq;
}

/** Prints the current solution to the terminal. */
void dopri54::print(double t_,double *in) {
	printf("%g",t_);
	for (int i=0;i<dof;i++) printf(" %g",in[i]);
	puts("");
}

/** Performs an integration step following the Dormand-Prince 5(4) scheme. */
void dopri54::step(double dt) {
	// Calculate k2.
	for (int i=0;i<dof;i++) dq[i]=q[i]+(dt/5.)*k1[i];
	f(t+dt/5.,dq,k2,p);

	// Calculate k3.
	for (int i=0;i<dof;i++) dq[i]=q[i]+(dt/40.)*(3*k1[i]+9*k2[i]);
	f(t+(3/10.)*dt,dq,k3,p);

	// Calculate k4.
	for (int i=0;i<dof;i++) {
		dq[i]=q[i]+dt*(44*k1[i]/45.-56*k2[i]/15.+32*k3[i]/9.);
	}
	f(t+(4/5.)*dt,dq,k4,p);

	// Calculate k5.
	for (int i=0;i<dof;i++) {
		dq[i]=q[i]+dt*(19372*k1[i]/6561.-25360*k2[i]/2187.+64448*k3[i]/6561.
			-212*k4[i]/729.);
	}
	f(t+(8/9.)*dt,dq,k5,p);

	// Calculate k6.
	for (int i=0;i<dof;i++) {
		dq[i]=q[i]+dt*(9017*k1[i]/3168.-355*k2[i]/33.+46732*k3[i]/5247.+49*k4[i]/176.
			-5103*k5[i]/18656.);
	}
	f(t+dt,dq,k6,p);

	// Calculate k7.
	for (int i=0;i<dof;i++) {
		dq[i]=q[i]+dt*(35*k1[i]/384.+500*k3[i]/1113.+125*k4[i]/192.-2187*k5[i]/6784.
			+11*k6[i]/84.);
	}
	f(t+dt,dq,k7,p);

	// Update the function count.
	fcount+=6;

	// The complete solutions. Here, q_hat is used for error estimation only.
	for (int i=0;i<dof;i++) {
		q[i]+=dt*(35*k1[i]/384.+500*k3[i]/1113.+125*k4[i]/192.-2187*k5[i]/6784.
			+11*k6[i]/84.);
		q_hat[i]+=dt*(5179*k1[i]/57600.+7571*k3[i]/16695.+393*k4[i]/640.
			-92097*k5[i]/339200.+187*k6[i]/2100.+k7[i]/40.);
	}

	// Update the current time.
	t+=dt;

	// Reuse k7 as k1 of the next step. To do this, we will just swap k1<->k7.
	double *tmp=k1;k1=k7;k7=tmp;
}

/** Integrates an ODE system with a fixed time step following the Dormand-Prince
 * 5(4) scheme. */
void dopri54::fixed_solve(double t_end,int n_steps,bool verbose) {
	// Set the initial conditions.
	init();
	
	// Calculate the step size.
	double dt=(t_end-t)/n_steps;
	
	// Print the initial condition to the terminal, if specified by the user.
	if (verbose) print(t,q);

	// Calculate k1 and update the function count.
	f(t,q,k1,p);
	fcount++;
	
	for (int i=0;i<n_steps;i++) {
		// Perform an integration step and print the current solution to the
		// terminal, if specified by the user.
		step(dt);
		if (verbose) print(t,q);
	}
}

/** Updates the error estimate of the solution with respect to some reference
 * solution. The reference solution is typically obtained either analytically or
 * via some highly accurate integration scheme. */
void dopri54::error(double *q_,double *qh_) {
	double s=0.,del;
	for (int i=0;i<dof;i++) {
		del=q_[i]-qh_[i];
		s+=del*del;
	}
	abs_err=sqrt(s/dof);
}

/** Prints performance information about the solver. */
void dopri54::work_precision() {
	printf("%.15g %d\n",abs_err,fcount);
}
