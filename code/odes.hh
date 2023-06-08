#ifndef ODES_HH
#define ODES_HH

#include <fstream>
#include <string>

#include "dopri54.hh"

#include "blas.h"

class lotka_volterra {
public:
	/** Specifies the right-hand side of the Lotka-Volterra system. */
	virtual void lv_f(double t_,double *in,double *out,double *p) {
		out[0]=p[0]*in[0]-p[1]*in[0]*in[1];
		out[1]=-p[2]*in[1]+p[3]*in[0]*in[1];
	}
	/** Sets the initial condition for the Lotka-Volterra system. */
	void lv_init(double *q_,double *qh_) {
		for (int i=0;i<2;i++) q_[i]=qh_[i]=1.;
	}
};

class lv_dop54 : public dopri54, public lotka_volterra {
public:
	lv_dop54(double *p_,double t0=0.): dopri54(2,p_,t0) {}
	virtual void init() {lv_init(q,q_hat);}
	virtual void f(double t_,double *in,double *out,double *p_) {
		lv_f(t_,in,out,p_);
	}
};

class lv_sensitivity {
public:
	virtual void lvs_f(double t_,double *in,double *out,double *p) {
		// Equations that govern the state variables.
		out[0]=p[0]*in[0]-p[1]*in[0]*in[1];
		out[1]=-p[2]*in[1]+p[3]*in[0]*in[1];
		// Forward sensitivity equations.
		out[2]=(p[0]-p[1]*in[1])*in[2]-p[1]*in[0]*in[3]+in[0];
		out[3]=p[3]*in[1]*in[2]+(-p[2]+p[3]*in[0])*in[3];
		out[4]=(p[0]-p[1]*in[1])*in[4]-p[1]*in[0]*in[5]-in[0]*in[1];
		out[5]=p[3]*in[1]*in[4]+(-p[2]+p[3]*in[0])*in[5];
		out[6]=(p[0]-p[1]*in[1])*in[6]-p[1]*in[0]*in[7];
		out[7]=p[3]*in[1]*in[6]+(-p[2]+p[3]*in[0])*in[7]-in[1];
		out[8]=(p[0]-p[1]*in[1])*in[8]-p[1]*in[0]*in[9];
		out[9]=p[3]*in[1]*in[8]+(-p[2]+p[3]*in[0])*in[9]+in[0]*in[1];
	}

	/** Sets the initial conditions of the system, including the sensitivity
	 * equations. */
	void lvs_init(double *q_,double *qh_) {
		for (int i=0;i<10;i++) q_[i]=qh_[i]=i<2?1.:0.;
	}
};

class lvs_dop54 : public dopri54, public lv_sensitivity {
public:
	/** The number of state variables. */
	const int m;
	/** The number of parameters. */
	const int l;
	/** The loss. */
	double loss;
	/** The gradient of the loss with respect to the solution. */
	double *dldq;
	/** The gradient of the solution with respect to the parameters. */
	double *dqdp;
	/** The gradient of the loss with respect to the parameters. */
	double *dldp;
	
	/** The class constructor defines the constants and dynamically allocates
	 * memory. */
	lvs_dop54(double *p_,double t0=0.): dopri54(10,p_,t0), m(2), l(4), loss(0.),
		dldq(new double[m]), dqdp(new double[m*l]), dldp(new double[l]) {
		// Zero the gradient of the loss with respect to the parameters.
		for (int i=0;i<l;i++) dldp[i]=0.;
	}

	/** The class destructor frees the dynamically allocated memory. */
	~lvs_dop54() {
		delete [] dldp;
		delete [] dqdp;
		delete [] dldq;
	}
	
	/** Sets the initial conditions, including the sensitivity equations. */
	virtual void init() {lvs_init(q,q_hat);}
	
	/** Evaluates the right-hand side of the Lotka-Volterra problem. */
	virtual void f(double t_,double *in,double *out,double *p_) {
		lvs_f(t_,in,out,p_);
	}

	/** Updates the parameters using the sensitivity equations and gradient
	 * descent. */
	void assemble_grad(double *q_ref,int n) {
		// Calculate the error between the numerical solution and the true
		// solution and update the loss.
		for (int i=0;i<m;i++) {
			dldq[i]=q[i]-q_ref[i];
			loss+=(2./n)*dldq[i]*dldq[i];
		}

		// The gradient of the solution with respect to the parameters comes
		// from the solution to the sensitivity problem.
		for (int i=0;i<m*l;i++) dqdp[i]=q[i+2];
	
		// Assemble the gradient of the loss with respect to the parameters
		// using the chain rule. We will make use of BLAS's DGEMM routine for
		// the matrix multiplication. Keep in mind that BLAS stores matrices in
		// column-major order.
		const char trans_n='n',trans_t='t';
		const double alpha=1.;
		const int k=1;
		dgemm_(&trans_t,&trans_n,&k,&l,&m,&alpha,dldq,&m,dqdp,&m,&alpha,dldp,&k);
	}

	/** Solves the Lotka-Volterra problem including the sensitiviy equations
	 * with a fixed time step integration scheme and estimates the parameters of
	 * the system given a reference solution. */
	void solve(double t_end,int n_steps,std::string fname,double lr,
		bool verbose=false) {
		// Set the initial conditions.
		init();
		
		// Calculate the step size.
		double dt=(t_end-t)/n_steps;

		// Calculate k1 and update the function count.
		f(t,q,k1,p);
		fcount++;

		// Read data from the reference solution.
		std::ifstream f(fname);
		double t_ref,*q_ref=new double[m];
		f>>t_ref>>q_ref[0]>>q_ref[1];
		
		// Update the loss and the gradient.
		assemble_grad(q_ref,n_steps);
	
		for (int i=0;i<n_steps;i++) {
			// Perform an integration step and assemble the gradient of the loss
			// with respect to the parameters.
			step(dt);
			f>>t_ref>>q_ref[0]>>q_ref[1];
			assemble_grad(q_ref,n_steps);
		}
		
		// Print the loss to the terminal, if specified by the user.
		if (verbose) printf("%g\n",loss);

		// Update the parameters.
		for (int i=0;i<l;i++) p[i]=p[i]-lr*dldp[i];
	}
};

#endif
