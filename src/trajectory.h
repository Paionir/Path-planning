#ifndef TRAJ_H
#define TRAJ_H

#include <cstdio>
#include <cassert>
#include <vector>
#include <math.h>
#include <algorithm>


using namespace std;

class Traj 
{
	// This class describes a unidimension trajectory in s or d
public:
	Traj() {};
	Traj(vector<double> coeffs, double T_f); // builds a traj using coefficients and final time T_f
	~Traj();
	void setCost(double C); // set cost for Traj
	vector<double> getCoeffs(); // returns coefficinets of trajectory
	double getCost(); // returns cost
	double getDis(double t); // returns position at time t
	double getVel(double t); // returns velocity at time t
	double getAcc(double t); // returns acceleration at time t
	double getJer(double t); // returns Jerk at time t
	double getJerAvg(double t1, double t2); // returns average Jerk between time t1 and t2

	void operator = ( const Traj &Tr ); // initialize a Traj as equal to another Traj

	double a0;
	double a1;
	double a2;
	double a3;
	double a4;
	double a5;
	double T;
	double max_d;
	double max_Dd;
	double max_v;
	double max_a;
	double max_J;
	double avg_J;
	double Cost;
	bool CostIsSet;
	bool isVelCost;

private:
	void findMax (); // find max of position, velocity, acceleration and jerk
	
};

class combiTraj
{
	// This class describes a bi dimensional trajectory created combining two unidimensional trajectories
	public:
		combiTraj() {};
		combiTraj ( Traj T1s, Traj T1d, double Tf ) ; // combines two Traj
		~combiTraj(){}; 

		bool collision(vector<double> params, double Tf); // checks for collisions, return true or false
		int dynamic(double vel_limit, double acc_limit, double Jer_limit); // checks for dynamic limits trespass (1 velocity, 2 acceleration, 3 Jerk)

		Traj Trs;
		Traj Trd;
		double T;
		
		double max_v;
		double max_a;
		double max_J;
		double Cost;


	private:
		void findMax();  // find max of position, velocity, acceleration and jerk
}; 

combiTraj::combiTraj(Traj T1s, Traj T1d, double Tf)
{
	Trs = T1s;
	Trd = T1d;
	T = Tf;
	Cost = T1s.Cost + T1d.Cost;   // set combined cost as a sum

	findMax();
}

int combiTraj::dynamic(double vel_limit, double acc_limit, double Jer_limit)
{
	double dt = 0.02;
	// check for velocity limit
	for (int i = 0; i < T/dt; i++)
	{
		if (sqrt(pow(Trs.getVel(i*dt),2.) + pow(Trd.getVel(i*dt),2.)) > 1.05*vel_limit) 
		{
		
			return 1;
		}
	}

	double avga = 0;
	int numa=0;
	dt = 0.02;
	// check for acceleration limit
	for (int i = 0; i < T/dt; i++)
	{	
		avga+=(sqrt(pow(Trs.getAcc(i*dt),2.) + pow(Trd.getAcc(i*dt),2.)));
		numa++;
		if (numa==10)
		{
			if (avga/numa>acc_limit)
			{
				return 2;
			}
			numa--;
			avga-= (sqrt(pow(Trs.getAcc((i-numa)*dt),2.) + pow(Trd.getAcc((i-numa)*dt),2.)));
		}
	}

	
	dt = 1;
	//check for Jerk limit
	for (int i = 0; i < T/dt -1 ; i++)
	{
		if ( sqrt( pow(Trs.getJerAvg(i*dt, (i+1)*dt),2.) + pow(Trd.getJerAvg(i*dt, (i+1)*dt),2.) ) > Jer_limit) 	
		{
			return 3;
		}	
	}
		
	
	return 0;
}

bool combiTraj::collision(vector<double> params, double Tf) 
{
	double max = 6945.554;

	//our position at time 0
	double s_0 = Trs.getDis(0.);
	double d_0 = Trd.getDis(0.);

	// our position at final time
	double s_h = Trs.getDis(Tf);
	double d_h = Trd.getDis(Tf);

	double sc_0 = params[0]; // near car position s
	double dc_0 = params[1]; // near car position d
	double vc_0 = params[2]; // near car velocity

	if (abs(s_0  - sc_0) > max/2.)  // check for lap end problems
	{
		if (s_0 > max/2.)
		{
			sc_0 += max;
		}
		else
		{
			sc_0 -= max;
		}
	}

	double dd_h = abs(dc_0 - d_h); // diff of intial and final d position

	double dt = 0.02; 

	// check for future collisions
	for (int i=0; i < Tf/dt ; i++)
	{
		//our velocity at time i*dt
		double vs_p = Trs.getVel(i*dt); 
		double vd_p = Trd.getVel(i*dt);
		
		//our posiiton at time i*dt
		double d_p = Trd.getDis(i*dt); 
		double s_p = Trs.getDis(i*dt);
		
		// our heading at time i*dt
		double head = atan2(vs_p, vd_p); 
		
		double v_p = sqrt( vs_p*vs_p + vd_p*vd_p  );
		
		double dd_max = abs(d_p - d_h);

		// position of near car at time i*dt
		double sc_p = sc_0 + vc_0*(i*dt);
		double dc_p = dc_0;

		// distance from near car at time dt*i
		double ds_p = s_p - sc_p;
		double dd_p = abs(d_p - dc_p);

		if (dd_p <= 2)	//car is in our lane
		{
			if (ds_p <= 0) // car is ahead of us
			{
				if (abs(ds_p) <= Tf/3.*vs_p && (dd_max < 2.)) // we want to stay in this lane and distance is too little
				{
					
					return true;
				}
				if (abs(ds_p) <= Tf/20.*vs_p && (dd_max > 2.)) // we want to change lane and distance is too short to change
				{
					
					return true;
				}
			}

		}
		else if (dd_p > 2.) // car is not in our lane
		{
			if(ds_p > 0) //car is behind
			{
				if (abs(ds_p) <= 0.5*vs_p && (dd_h < 2. || dd_max > 7.) ) // we are changing lane toward the cars
				{	
					return true;
				}
			}
			else if(ds_p <= 0) //car is ahead
			{
				if ( abs(ds_p) <= Tf/3.*vs_p && (dd_h < 2. || dd_max > 7.) ) // we are changing lane toward the cars
				{
					return true;	
				}
			}
		}
	}
	return false;
}

void combiTraj::findMax()
{
	double v_m = 0.;
	double a_m = 0.;
	double J_m = 0.; 

	for (int i=0; i<int(T/0.02); i++)
	{
		double t = i*0.02;
		double vel = sqrt( pow(Trs.getVel(t),2.) + pow(Trd.getVel(t),2.) );
		double acc = sqrt( pow(Trs.getAcc(t),2.) + pow(Trd.getAcc(t),2.) );
		double Jer = sqrt( pow(Trs.getJer(t),2.) + pow(Trd.getJer(t),2.) );
		
		if (v_m < vel) {v_m = vel;}
		if (a_m < acc) {a_m = acc;}
		if (J_m < Jer) {J_m = Jer;}
		
		
	}
	max_v = v_m;
	max_a = a_m;
	max_J = J_m;
	
}



Traj::Traj(std::vector<double> coeffs, double T_f)
{
	a0 = coeffs[0];
	a1 = coeffs[1];
	a2 = coeffs[2];
	a3 = coeffs[3];
	a4 = coeffs[4];
	if (coeffs.size()>5)
	{
		a5 = coeffs[5];
		isVelCost = false;	
	} else {
		
		a5=0.;
		isVelCost = true;
	}
	

	T=T_f;

	avg_J = sqrt((3.*a3*a3*T + 12.*a3*a4*T*T + (16.*a4*a4 + 20.*a3*a5)*T*T*T + 
			60.*a4*a5*T*T*T*T + 60.*a5*a5*T*T*T*T*T )/T);
	
	findMax();
	CostIsSet = false;
	Cost = -1.;
	
}

Traj::~Traj()
{

}

void Traj::findMax() 
{	
	double d_m = a0;
	double Dd_m = 0.;
	double v_m = 0.;
	double a_m = 0.;
	double J_m = 0.; 

	for (int i=0; i<int(T/0.02); i++)
	{
		double t = i*0.02;
		double dh = a0 + a1*t + a2*t*t + a3*t*t*t + a4*t*t*t*t + a5*t*t*t*t*t;
		double dis = abs(dh);
		double Dd = abs(dh - a0);
		double vel = abs(a1 + 2.*a2*t + 3.*a3*t*t + 4.*a4*t*t*t+ 5.*a5*t*t*t*t);
		double acc = abs(2.*a2 + 6.*a3*t + 12.*a4*t*t + 20.*a5*t*t*t);
		double Jer = abs(6.*a3 + 24.*a4*t + 60.*a5*t*t);
		if (d_m < dis) {d_m = dis;}
		if (v_m < vel) {v_m = vel;}
		if (a_m < acc) {a_m = acc;}
		if (J_m < Jer) {J_m = Jer;}
		if (abs(Dd_m) < Dd) {Dd_m = dh-a0;}
		
	}
	max_d = d_m;
	max_v = v_m;
	max_a = a_m;
	max_J = J_m;
	max_Dd = Dd_m;
}

void Traj::setCost(double C) 
{
	Cost = C;
	CostIsSet = true;	
}

double Traj::getCost()
{
	if (CostIsSet)
	{
		return Cost;
	}
	else 
	{
		return -1.;
	}
}

vector<double> Traj::getCoeffs()
{
	return{a0,a1,a2,a3,a4,a5};
}



double Traj::getDis(double t)
{
	if (t > T)
	{
		if (isVelCost)
		{
			double disf = a0 + a1*T + a2*T*T + a3*T*T*T + a4*T*T*T*T + a5*T*T*T*T*T;
			double velf = a1 + 2.*a2*T + 3.*a3*T*T + 4.*a4*T*T*T + 5.*a5*T*T*T*T;

			return (disf + velf*(t-T));	
		}
		else 
		{
			return (a0 + a1*T + a2*T*T + a3*T*T*T + a4*T*T*T*T + a5*T*T*T*T*T);	
		}	
	}
	else
	{
		return (a0 + a1*t + a2*t*t + a3*t*t*t + a4*t*t*t*t + a5*t*t*t*t*t);
	}
}

double Traj::getVel(double t)
{
	if (t > T)
	{
		if (isVelCost)
		{
			
			double velf = a1 + 2.*a2*T + 3.*a3*T*T + 4.*a4*T*T*T + 5.*a5*T*T*T*T;
			return (velf);	
		}
		else 
		{
			return (0.);	
		}	
	}
	else
	{
		return (a1 + 2.*a2*t + 3.*a3*t*t + 4.*a4*t*t*t+ 5.*a5*t*t*t*t);	
	}
	
}

double Traj::getAcc(double t)
{
	if (t > T)
	{
		return (0.);	
	}
	else
	{
		return (2.*a2 + 6.*a3*t + 12.*a4*t*t + 20.*a5*t*t*t);
	}
	
}

double Traj::getJer(double t)
{
	if (t > T)
	{
		return(0.);
	}
	else
	{
		return (6.*a3 + 24.*a4*t + 60.*a5*t*t);
	}
}

double Traj::getJerAvg(double t1, double t2)
{
	if (t2==t1) {return 0.;}

	return (sqrt(
			(3.*a3*a3*t2 + 12.*a3*a4*t2*t2 + (16.*a4*a4 + 20.*a3*a5)*t2*t2*t2 + 
			60.*a4*a5*t2*t2*t2*t2 + 60.*a5*a5*t2*t2*t2*t2*t2 )
			-
			(3.*a3*a3*t1 + 12.*a3*a4*t1*t1 + (16.*a4*a4 + 20.*a3*a5)*t1*t1*t1 + 
			60.*a4*a5*t1*t1*t1*t1 + 60.*a5*a5*t1*t1*t1*t1*t1 )
			/
			(t2-t1)
		));
}

void Traj::operator = (const Traj &Tr )
{
	a0 = Tr.a0;
	a1 = Tr.a1;
	a2 = Tr.a2;
	a3 = Tr.a3;
	a4 = Tr.a4;
	a5 = Tr.a5;
	T = Tr.T;
	max_d = Tr.max_d;
	max_Dd = Tr.max_Dd;
	max_v = Tr.max_v;
	max_a = Tr.max_a;
	max_J = Tr.max_J;
	avg_J = Tr.avg_J;
	CostIsSet = Tr.CostIsSet;
	Cost = Tr.Cost;
	isVelCost = Tr.isVelCost;

}

#endif