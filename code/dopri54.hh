#ifndef DOPRI54_HH
#define DOPRI54_HH

#include <cstdio>
#include <cmath>

class dopri54 {
protected:
	/** The degrees of freedom. */
	int dof;
	/** The current time. */
	double t;
	/** A counter to keep track of the function evaluations. */
	int fcount;
	/** The absolute error with respect to some reference solution. */
	double abs_err;
	/** A temporary array to store intermediate results. */
	double *dq;
	/** The intermediate steps. */
	double *k1,*k2,*k3,*k4,*k5,*k6,*k7;
public:
	/** The solution vector. */
	double *q;
	/** An alternative solution used for error estimation. */
	double *q_hat;
	/** The parameters of the system. */
	double *p;
	dopri54(int dof_,double *p_,double t0=0.);
	~dopri54();
	void print(double t_,double *in);
	virtual void init()=0;
	virtual void f(double t_,double *in,double *out,double *p_)=0;
	void step(double dt);
	void update_params();
	void fixed_solve(double t_end,int n_steps,bool verbose=false);
	void error(double *q_,double *qh_);
	void work_precision();
};

#endif
