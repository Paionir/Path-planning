#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "Eigen-3.3/Eigen/Dense"
#include "spline.h"
#include "trajectory.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2.*pi() - angle, angle);

  if(angle > 2.*pi()/3.)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}



// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}
   

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// ADDITIONAL HELPER FUNCTIONS

double calcPoly (std::vector<double> coeffs, double t) 
{
  // This functions calculates the value of a polynomial at time t and returns its value

  //INPUT:
  // a vector of all coefficients for the polynomial sorted from lowest degree to highest
  // the time t at which evaluate the polynomial

  //OUTPUT:
  // the value of the polynomial at time t
  
  double pol = 0.;
  for (int i = 0; i < coeffs.size(); i++)
  {
    pol += coeffs[i] * pow(t, i);
  }
  return pol;
}


std::vector<double> parabolicInterpol(std::vector<double> X, std::vector<double> Y, int center, double ds, double d) 
{
  // This functions interpolates a 2nd grade polynomial between 3 waypoints X,Y and then uses ds and d to estimate the x,y position
  // of a point between the center waypoint and the next

  //INPUT:
  // vector of X coordinates of 3 waypoints
  // vector of Y coordinates of 3 waypoints
  // index of the central waypoint
  // arc lenght along the parabola ds
  // coordinate d 

  //OUTPUT:
  // a vector of (x,y) coordinate for point (ds,d)



  // transform to reference system of center point
  double x0 = X[center-1] - X[center]; 
  double x1 = X[center] - X[center];
  double x2 = X[center+1] - X[center];
  
  double y0 = Y[center-1] - Y[center];
  double y1 = Y[center] - Y[center];
  double y2 = Y[center+1] - Y[center];
    
  double den_X = (x0-x1)*(x0-x2)*(x1-x2);
  double den_y = (y0-y1)*(y0-y2)*(y1-y2);
  double disc_x = (x0-x1)*(x1-x2);
  double disc_y = (y0-y1)*(y1-y2);
  bool rotate = false;
  
  if (disc_x <= 0 )  
  {
    //rotate reference system, so that (x,y) -> (-y,x)
 
  double tx0 = -y0;
  double tx1 = -y1;
  double tx2 = -y2;
  
  y0 = x0;
  y1 = x1;
  y2 = x2;

  x0 = tx0;
  x1 = tx1;
  x2 = tx2;

  std::vector<double> TX;

  for (int i =0;i<X.size();i++)
  {
    TX.push_back(-Y[i]);
    Y[i]=X[i];
    X[i]=TX[i];
  }

  rotate = true;

  }

  // Calculate 3 parameters of the parabola passing by the 3 waypoints y=ax^2+bx+c
  double den = (x0-x1)*(x0-x2)*(x1-x2);
  double a = ( x2*(y1-y0) + x1*(y0-y2) + x0*(y2-y1) )/den;
  double b = ( x2*x2*(y0-y1) + x1*x1*(y2-y0) +x0*x0*(y1-y2) )/den;
  double c = ( x1*x2*(x1-x2)*y0 + x2*x0*(x2-x0)*y1 +x0*x1*(x0-x1)*y2 )/den;
  
  
  double sum = 0.;
  int idx = 0;

  double X1 = X[center]-X[center]; // transform to reference system of center point
  double X2 = X[center+abs(ds)/ds]-X[center]; // second integration limit is the previous or successive point of center, according to ds sign    
  
  double h = (X2-X1)/50.;

  // the arc lenght of a parabola is the definite integral of sqrt(1+f'(x)^2) in dx
  double u1 = 2.*a*X1 +b; // helper variable 
  double g1 = u1*sqrt(1+u1*u1) + log(abs(sqrt(1+u1*u1) + u1)); // primitive of sqrt(1+f'(x)^2) calculated in X1
  
  double xe2=X1;

  // EVALUATE xe2 at which the arc lenght equals |ds| with 1e-11 tolerance or with 10000000 max iterations, whatever happens first
  while (( abs(abs(ds) - sum) > 1e-11) && (idx < 10000000))
  {
    
    xe2 += h;
    double u2 = 2.*a*xe2 + b;
    double g2 = (u2*sqrt(1+u2*u2) + log(abs(sqrt(1+u2*u2) + u2))); // primitive of sqrt(1+f'(x)^2) calculated in xe2
    
    sum = abs((g2 - g1)/(4.*a)); // arc lenght from X1 to xe2
    if (sum > abs(ds) ) // if arc lenght is greater than |ds| go back one step and divide h by 2
    {
      xe2 -= h;  
      h = h/2.;  
    }
    idx++;
  }
  double xp = xe2;
  double yp = calcPoly({c,b,a},xp);
  double heading = atan2(2.*a*xp + b, 1.); //calculate heading of parabola at point (xp, yp=2axp+b)
 
  // transform back to global reference system
  xp += X[center];
  yp += Y[center];

  if (rotate)
  {
    //rotate back
    double txp= xp;
    xp = yp;
    yp = -txp;


    if (x1-x0 > 0.)
    {
      heading = heading + pi();
    }

    // add d offset using heading
    
    xp += d * cos(heading);
    yp += d * sin(heading);


  } else {
      

    if (x1-x0 < 0.)
    {
      heading = heading + pi();
     }
    heading = heading-pi()/2.;
    xp += d * cos(heading);
    yp += d * sin(heading);  
  }

  return{xp,yp};
}


vector<double> parabolicGetXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
  // This functions transform from s,d coordinates to global x,y coordinates using a waypoint maps of the highway
  // Instead of a linear interpolation, it uses two parabolic interpolation and then calculates a weighted mean. 
  // The first interpolation is made using the previous waypoint and the immidiately successive and previous waypoint,
  // the second interpolation is made using the previous waypoint and the immidiately 2 successive waypoints.
  // Then a weighted mean of the two points is calculated using, as weights, the inverse of the squared distance from the 
  // previous waypoint and the next waypoint

  //INPUT:
  // s coordinate
  // d coordinate
  // s values of waypoints
  // x values of waypoints
  // y values of waypoints

  //OUTPUT:
  // a vector of (x,y) coordinate for point (s,d)

  double max_s = 6945.554; //max s value for waypoints
  while (s > max_s)
  {
    s -= max_s;
  }

  int prev_wp = -1;

  // find the previous waypoint

  while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
  {
    prev_wp++;
  }

 
  int cubic_num; //helper index of the previous waypoint

  vector<double> X; // X coordinates of nearest waypoints used for interpolation
  vector<double> Y; // Y coordinates of nearest waypoints used for interpolation

  // fill X and Y with 1 previous waypoint and 2 successive waypoints,
  // if previous waypoint is 0 then start from the last waypoint in the map

  if (prev_wp >=1) 
  {
    cubic_num = 1;
    for (int i = -1; i < 3; i++)
    {
      X.push_back( maps_x[(prev_wp + i)%maps_x.size()] );
      Y.push_back( maps_y[(prev_wp + i)%maps_x.size()] );
    }
  } 
  else 
  {
    cubic_num = 1;
    for (int i = maps_x.size() -1 ; i < maps_x.size() + 3; i++)
    {
      X.push_back( maps_x[i%maps_x.size()] );
      Y.push_back( maps_y[i%maps_x.size()] );
    }
  }

  double ds_p = s - maps_s[prev_wp]; //distance in s from previous waypoint

  std::vector<double> XYp = parabolicInterpol(X,Y, cubic_num, ds_p, d); // calculate x,y using the previous waypoint as the central waypoint for interpolation
  
  double ds_s; // distance in s from the next waypoint
  if (prev_wp == maps_s.size() - 1 )
  {
    ds_s = s - max_s;  
  }
  else
  {
    ds_s = s - maps_s[(prev_wp+1)];  
  }
  

  std::vector<double> XYs = parabolicInterpol(X,Y, cubic_num+1, ds_s, d); // calculate x,y using the next waypoint as the central waypoint for interpolation

  // calculate the weighted mean of the two interpolations using the inverse sqaure of the distance from previous and next waypoint
  int n_exp=-2;
  double p1 = pow(ds_p,n_exp);
  double p2 = pow(ds_s,n_exp);
  double norm =p1+p2;
  double x = (XYp[0]*p1 + XYs[0]*p2)/(norm);
  double y = (XYp[1]*p1 + XYs[1]*p2)/(norm);  
 
  return {x,y};
   
}  


vector<double> keepVelPoly(vector<double> conds, double T)
{
  // This function calculates the jerk-minimizing trajectory that tries to keep a desired speed. 
  // The optimal trajectory is a quartic polynomial in this case.

  // INPUT:
  //  conds is a vector containing 5 boundary conditions (s,s' and s'' at time 0 and s', s'' at time T)
  //  T is the final time for the end of trajectory

  // OUTPUT:
  // a vector of 5 coefficients for the quartic polynomial

  double s_i = conds[0];
  double s1_i = conds[1];
  double s2_i = conds[2];
  double s1_f = conds[3];
  double s2_f = conds[4];

  double A1 = s1_f - (s1_i + s2_i*T);
  double A2 = s2_f-s2_i;

  vector<double> coeffs;
  
  coeffs.push_back(s_i);
  coeffs.push_back(s1_i);
  coeffs.push_back(s2_i*0.5);
  coeffs.push_back( (3.*A1 - A2*T)/(3.*T*T) );
  coeffs.push_back( (A2*T - 2.*A1)/(4.*T*T*T) );

  return coeffs;
}

vector<double> minJerkPoly(vector<double> conds, double T)
{
  // This function calculates the jerk-minimizing trajectory that tries to reach a final value d_f for d. 
  // The optimal trajectory is a quintic polynomial in this case.

  // INPUT:
  //  conds is a vector containing 6 boundary conditions (d,d' and d'' at time 0 and d, d', d'' at time T)
  //  T is the final time for the end of trajectory

  // OUTPUT:
  // a vector of 6 coefficients for the quintic polynomial

  double d_i = conds[0];
  double d1_i = conds[1];
  double d2_i = conds[2];
  double d_f = conds[3];
  double d1_f = conds[4];
  double d2_f = conds[5];

  double A1 = (d_f - d_i - d1_i*T - 0.5*d2_i*T*T);
  double A2 = (d1_f - d1_i - d2_i*T);
  double A3 = d2_f - d2_i;

  vector<double> coeffs;
  
  coeffs.push_back(d_i);
  coeffs.push_back(d1_i);
  coeffs.push_back(d2_i*0.5);
  coeffs.push_back( (20.*A1 - 8.*A2*T + A3*T*T)/(2.*T*T*T) );
  coeffs.push_back( (-15.*A1 + 7.*A2*T - A3*T*T)/(T*T*T*T) );
  coeffs.push_back( (12.*A1 - 6.*A2*T + A3*T*T)/(2.*T*T*T*T*T) );

  return coeffs;
}


// generate set of trajectories s(t) and d(t)
vector< vector<Traj> > genTrajSet (vector<double> conds_s, vector<double> conds_d,  double time_horizon, 
                      double s_goal, double l_desired, vector<double> limits, vector<double> & max_min)
{

  // This function generates a full set of unidimensional trajectories for s and d using jerk minimizing polynomials
  // and different values of final boundary conditions.
  // 
  // This function uses the class Traj that can be found in trajectory.h header file

  // INPUT:
  //  conds_s is the set of boundary conditions for s trajectories
  //  conds_d is the set of boundary conditions for d trajectories
  //  time_horizion is the final time for the end of trajectories
  //  s_goal is the desired speed for the car 
  //  l_desired is the desired lane for the car
  //  limits is a vector of dynamic limits (velocity, acceleration and jerk)
  //  max_min is a vector to store useful variables for cost normalization

  // OUTPUT:
  // a vector of 2 sets of trajectories, the first is the s trajectory set and the second is the d trajectory set

  vector<double> conds;
  vector<double> d_conds;

  double speed_limit = limits[0];
  double acc_limit = limits[1];
  double Jer_limit = limits[2];
  double speed_goal = s_goal;
  double speed_minimum = 8.;
  double lane_desired = l_desired;


  vector<Traj> longSet; // return set
  vector<Traj> lateSet; // return set

  // a set of useful variables for cosst normalization
  double max_ds = 0.;
  double min_ds = 10e10;
  double max_avgJ = 0.; 
  double min_avgJ = 10e10; 
  double max_T = 0.;
  double min_T = 10e10;

  double max_dd = 0.;
  double min_dd = 10e10;
  double d_max_avgJ = 0.; 
  double d_min_avgJ = 10e10; 
  double d_max_T = 0.;
  double d_min_T = 10e10;


  // create s trajectories using different time of manouver and final speed
  for (int i=0 ; i < 7 ; i++)
  {
    for (int j=0; j < 6 ; j++)
    {
      conds.push_back(conds_s[0]);
      conds.push_back(conds_s[1]);
      conds.push_back(conds_s[2]);
      conds.push_back(speed_limit - (speed_limit - speed_minimum)*i/6. );
      conds.push_back(conds_s[4]);
      
      double Tj = time_horizon;
      double dtj = (Tj - 1.5)/5.;
      
      if (Tj-dtj*j > 0)
      {
        Traj s_traj = Traj(keepVelPoly(conds, Tj - j*dtj), Tj -j*dtj);
        conds.clear();
        
          longSet.push_back(s_traj);
          // useful variables for cost normalization
          double ds = abs(s_traj.getVel(s_traj.T) - speed_limit);
          max_ds = max(max_ds,ds);
          min_ds = min(min_ds,ds);
          max_avgJ = max(max_avgJ, s_traj.avg_J);
          min_avgJ = min(min_avgJ, s_traj.avg_J);
          max_T = max(max_T, s_traj.T);
          min_T = min(min_T, s_traj.T);                      
      }
    }  
  }

  // create d trajectories using different time of manouver and final lane
  vector<double> lanes = {2., 6., 9.5};
  for (int i=0 ; i < 3 ; i++)
  {
    for (int j=0; j < 3 ; j++)
    {
         
      d_conds.push_back(conds_d[0]);
      d_conds.push_back(conds_d[1]);
      d_conds.push_back(conds_d[2]);
      d_conds.push_back(lanes[i]);
      d_conds.push_back(0.);
      d_conds.push_back(0.);

      double Tj = time_horizon-1;
      double dtj = (Tj - 2)/2.;
      
      if (Tj -dtj*j > 0)
      {
        Traj d_traj = Traj(minJerkPoly(d_conds, Tj - dtj*j), Tj - dtj*j);
        d_conds.clear();

          lateSet.push_back(d_traj);
          // useful variables for cost normalization
          double dd = abs(d_traj.getDis(d_traj.T) - lane_desired*4. - 2.);
          max_dd = max(max_dd,dd);
          min_dd = min(min_dd,dd);
          d_max_avgJ = max(d_max_avgJ, d_traj.avg_J);
          d_min_avgJ = min(d_min_avgJ, d_traj.avg_J);
          d_max_T = max(d_max_T, d_traj.T);
          d_min_T = min(d_min_T, d_traj.T);                      
      }
    }  
  }

  // store the max and min of d and s trajectories set for cost normalization
  max_min.push_back(max_ds);
  max_min.push_back(min_ds);
  max_min.push_back(max_avgJ);
  max_min.push_back(min_avgJ);
  max_min.push_back(max_T);
  max_min.push_back(min_T);

  max_min.push_back(max_dd);
  max_min.push_back(min_dd);
  max_min.push_back(d_max_avgJ);
  max_min.push_back(d_min_avgJ);
  max_min.push_back(d_max_T);
  max_min.push_back(d_min_T);
  
  return {longSet, lateSet};
}





vector<combiTraj> combineTrajectories(vector<Traj> longSet, vector<Traj> lateSet, double time_horizon, 
                      double s_goal, double l_desired, vector<double> limits, vector<double> & max_min, vector< vector<double> > & near_cars)
{

  // This function combines two set of trajectories into a set of combined bidimensional trajectories and assign to every combined trajectory a cost.
  // In addition this function checks for collisions and dynamic limits of every combined trajectory.

  // This function heavily uses the class combiTraj and its methods which can be found in trajectory.h header file

  // INPUT:
  //  longSet is the set of s trajectories
  //  lateSet is the set of d trajectories
  //  time_horizion is the final time for the end of combined trajectories
  //  s_goal is the desired speed for the car 
  //  l_desired is the desired lane for the car
  //  limits is a vector of dynamic limits (velocity, acceleration and jerk)
  //  max_min is a vector to store useful variables for cost normalization
  //  near_cars is the set of the nearest cars against which trajectories are checked for collisions



  double speed_limit = limits[0];
  double acc_limit = limits[1];
  double Jerk_limit = limits[2];
  double speed_goal = s_goal;
  double lane_desired = l_desired;

  // setting of max and min variables for cost normalization
  double max_ds = max_min[0];
  double min_ds = max_min[1];
  double max_avgJ = max_min[2]; 
  double min_avgJ = max_min[3]; 
  double max_T = max_min[4];
  double min_T = max_min[5];

  double max_dd = max_min[6];
  double min_dd = max_min[7];
  double d_max_avgJ = max_min[8]; 
  double d_min_avgJ = max_min[9]; 
  double d_max_T = max_min[10];
  double d_min_T = max_min[11];

  vector<combiTraj> combSet; // return set
  vector<int> dyn_rej = {0,0,0}; // helepr variable to count dynamic rejections
  int coll_rej =0; // helepr variable to count collision rejections

  // set cost for every accepted trajectory and combine them
  for (int k = 0; k<longSet.size(); k++)
  {
    for (int h = 0; h<lateSet.size(); h++)
    {
      if (h==0)
      {
        double ds = abs(longSet[k].getVel(longSet[k].T) - speed_limit); 
        ds = (ds - min_ds) / (max_ds - min_ds);   // normalized distance from desired speed
        double Tc = (longSet[k].T - min_T) / (max_T - min_T); // normalized time to complete manouver
        double Jc = (longSet[k].avg_J - min_avgJ) / (max_avgJ - min_avgJ); // normalized average Jerk
        longSet[k].setCost( (10.)*Jc + (10.)*Tc + (100.)*ds ); // normalized s cost
      }
      
      if (k==0)
      {
        double dd = abs(lateSet[h].getDis(lateSet[h].T) - lane_desired*4. - 2.);
        dd = (dd - min_dd) / (max_dd - min_dd); // normalized distance from desired lane
        double Tc = (lateSet[h].T - d_min_T) / (d_max_T - d_min_T); // normalized time to complete manouver
        double Jc = (lateSet[h].avg_J - d_min_avgJ) / (d_max_avgJ - d_min_avgJ); // normalized average Jerk
        lateSet[h].setCost((10.)*Jc + (10.)*Tc + (10.)*dd); // normalized d cost 
      }

      combiTraj comb_traj = combiTraj(longSet[k], lateSet[h], time_horizon); // combine 2 trajectories into 1 combiTraj
      int wrong_dyn = comb_traj.dynamic(speed_limit, acc_limit, Jerk_limit); // check for dynamic limit trespass

      if (wrong_dyn>0) {dyn_rej[wrong_dyn-1]++;}
      
      if ( wrong_dyn == 0 )  // if combiTraj survives dynamic check, check for collisions
      {
        bool coll = false;
        int idx = 0;

        while ( (idx < near_cars.size()) && (!coll) )
        {
          // for every near car get s,d and velocity
          double nc_s = near_cars[idx][5];
          double nc_d = near_cars[idx][6];
          double nc_vx = near_cars[idx][3];
          double nc_vy = near_cars[idx][4];
          double nc_v = sqrt(nc_vx*nc_vx + nc_vy*nc_vy);

          coll = (comb_traj.collision({nc_s,nc_d,nc_v}, time_horizon)); // use collision method to check for future collisions
          idx++;
        }
        if (!coll) 
        { 
          combSet.push_back(comb_traj); // if combiTraj survives the collisions check, add to the eligible trajectories set

        } else 
        {
         
          coll_rej++;

        }  
      }
      else
      {

      }
    }
  }

  return combSet;
}

bool funcia (combiTraj i, combiTraj j) { return (i.Cost < j.Cost); } // small helper function to sort a set of combiTraj by their costs
    
int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;
  
  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
  
  int size_prev_plan = 0; // number of points already reached since last planned path
  int size_prev_path = 0; // total size of last path passed to simulator
  int size_kept = 0; // points of previous path added to current path
  double speed_goal = 10e2; //desired speed 
  Traj longTrajectory; // trajectory along s
  Traj lateralTrajectory; //trajectory along d

    
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &size_prev_plan, &size_prev_path, &size_kept, &max_s, &speed_goal, 
                &longTrajectory, &lateralTrajectory](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            // Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            std::vector<std::vector<double>> near_cars;

            
            // cout <<"#near_cars: "<< near_cars.size()<<endl;
          	json msgJson;
            
          	vector<double> next_x_vals;
          	vector<double> next_y_vals;
           
            


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds 
              
              // const short int lane_desired = 2; //preferred lane to drive (range between 0 and 2 where 0 is the left lane)
              int lane_desired = floor(car_d/4);
              const double speed_limit = 22.352-2; // 22.352 ms/ is equal to 50Mph in a smarter units system (sorry USA) 2.2352
              const double acc_limit = 10.0; // max acceleration in m/2^2
              const double Jerk_limit = 10.0; // max Jerk in m/s^3
              const int size_horizon = 250; // size of path to pass to simulator for each new path 
              const int size_plan = 10; // size of path already driven after which a new path must be planned 
              const int size_keep = 0; // points of previous path to add to new path
              
             
              double time_horizon = (size_horizon - 1) * 0.02; // seconds for time horizon of path
              double time_plan =  (size_plan -1) * 0.02; // seconds between each path plannings
              int ncf = 0;
              int ncb = 0;
              speed_goal = min(speed_limit, speed_goal); // desired velocity for car

              
              // select the nearest cars from sensor fusion data
              for (int i=0; i< sensor_fusion.size(); i++)
              {
                double sc = sensor_fusion[i][5];
                if (abs(car_s  - sc) > max_s/2.) // check if the car is near the start of the track and adjust s values
                {
                  if (car_s > max_s/2.)
                  {
                    sc += max_s;
                  }
                  else
                  {
                    sc -= max_s;
                  }
                }
                double dc = sensor_fusion[i][6];
                double vx = sensor_fusion[i][3];
                double vy = sensor_fusion[i][4];
                double vc = sqrt(vx*vx + vy*vy);

                if (car_s > sc) // car is behind us
                {
                  // check if the car is at a distance such that it can reach our car in time_horizon*0.7 seconds at current speed
                  if (vc*time_horizon + sc >= car_s + car_speed*0.44704*time_horizon*0.7 - 4.) 
                  {
                    near_cars.push_back(sensor_fusion[i]);
                    ncb++;
                  }
                }
                else //car is ahead of us
                {
                  // check if the car is at a distance that our car can reach in time_horizon seconds at max speed, even if the car decelerates to 1/3 of its current velocity
                  if (vc*time_horizon/3. + sc <= car_s + car_speed*0.44704*time_horizon + 4.)
                  {
                    near_cars.push_back(sensor_fusion[i]);
                    ncf++;
                  } 
                }
              }

              
              double pos_x;
              double pos_y;
              double angle;
              int path_size = previous_path_x.size(); // how many points where passed in previous plan

              double next_s; // next s value
              double next_d; // next d value
              vector<double> sxy;
              
              std::vector<double> conds; // boundary conditions for s
              std::vector<double> d_conds; // boundary conditions for d


             
              if(path_size == 0) // intialize first path
              {
                // set intial s and d conditions
                conds.clear();
                conds.push_back(car_s);
                conds.push_back(0.);
                conds.push_back(0.);
                conds.push_back(speed_goal);
                conds.push_back(0.);
              
                d_conds.clear();
                d_conds.push_back(car_d);
                d_conds.push_back(0.);
                d_conds.push_back(0.);
                d_conds.push_back(0*4.+ 2.);
                d_conds.push_back(0.);
                d_conds.push_back(0.);


                vector<Traj> longSet;
                vector<Traj> lateSet;
                vector<double> max_min;
                
                //generate set of unidimensional trajectories                
                vector < vector<Traj> > sets = genTrajSet (conds, d_conds, time_horizon, speed_goal, lane_desired, {speed_limit,10.,10.}, max_min);
                
                // set cost for every trajectory and combine them. Check also for collisions and dynamic limits
                vector <combiTraj> combSet = combineTrajectories(sets[0], sets[1], time_horizon, speed_goal, lane_desired, {speed_limit,10.,10.}, max_min, near_cars); 

                // find minimal cost trajectory
                double min_Comb_Cost = 10e10;
                int min_Comb_idx = 0;
                for (int k = 0; k < combSet.size(); k++)
                {
             
                  if (min_Comb_Cost > combSet[k].Cost) 
                  {
                    min_Comb_Cost = combSet[k].Cost;
                    min_Comb_idx = k;
                  }
                }
                                
                longTrajectory = combSet[min_Comb_idx].Trs;  // set s trajectory 
                lateralTrajectory = combSet[min_Comb_idx].Trd; // set d trajectory
              
                size_prev_path=0;

                //generate next points using selected trajectory with a time pace of 0.02 seconds
                for (int i =0; i < size_horizon; i++)
                {
                  next_s = longTrajectory.getDis(i*0.02); // get s value at time i*0.02
                  next_d = lateralTrajectory.getDis(i*0.02); // get d value at time d*0.02
                  
                  // convert  to  global coordinates
                  sxy = parabolicGetXY(next_s, next_d,  map_waypoints_s, map_waypoints_x, map_waypoints_y); 
                  
                  // pass to simulator
                  next_x_vals.push_back(sxy[0]); 
                  next_y_vals.push_back(sxy[1]);
                
                  size_prev_path++;
                }
                size_kept = 0;
              } 
              else 
              {
                  // CURRENT PATH IS VALID UNTIL NEXT PLAN SO CHECK TIME ELAPSED FROM PREVIOUS PLANNING
            
                  size_prev_plan = size_prev_path - path_size;
                 
                 if (size_prev_plan >= size_plan) { 
                    
                    // PLAN AGAIN AND RESET time_prev_path
                    
	                  size_prev_path = 0;

                    // KEEP points of previous path
                     for(int i = 0; i < size_keep ; i++)
                    {
                      next_x_vals.push_back(previous_path_x[i]);
                      next_y_vals.push_back(previous_path_y[i]);
                    
                      size_prev_path++;
                      
                    }

                    // prepare new boundary conditions
                    conds.clear();
                    d_conds.clear();
                    // point from wich to start new plan
                    int point_i = size_prev_plan + size_keep - size_kept + 1;
                    size_kept = size_keep;
                    double t_i = (point_i-1) * 0.02;
              

                    double ss_i = longTrajectory.getDis(t_i); // s position at time t_i
                    while(ss_i>max_s)
                    {
                      ss_i -= max_s;
                    }
                    double vs_i = longTrajectory.getVel(t_i); // velocity along s at time t_ti
                    double as_i = longTrajectory.getAcc(t_i); // acceleration along s at time t_i

                    

                    double dd_i = lateralTrajectory.getDis(t_i); // d position at time t_i
                    double vd_i = lateralTrajectory.getVel(t_i); // velocity along d at time t_ti
                    double ad_i = lateralTrajectory.getAcc(t_i); // acceleration along d at time t_i
                    
                    // push conditions
                    conds.push_back(ss_i);
                    conds.push_back(vs_i);
                    conds.push_back(as_i);
                    conds.push_back(speed_goal);
                    conds.push_back(0.);
                    
                    d_conds.push_back(dd_i);
                    d_conds.push_back(vd_i);
                    d_conds.push_back(ad_i);
                    d_conds.push_back(lane_desired*4.+ 2.);
                    d_conds.push_back(0.);
                    d_conds.push_back(0.);

                    

                    double t_s =  (size_horizon - size_keep -1)*0.02;
                    double t_d =  (size_horizon - size_keep -1)*0.02;
                    
                    vector <combiTraj> combSet; // set of combined trajectories
                    double time_manouver = t_s;
                    int min_Comb_idx = 0;

                    // continue to generate trajectories if no suitable trajectory is found, until time_manouver is lower than 2 seconds
                    while(combSet.size() < 1 && time_manouver > 2. ) 
                    {
                      vector<Traj> longSet;
                      vector<Traj> lateSet;
                      vector<double> max_min;
                      
                      //generate set of unidimensional trajectories                
                      vector < vector<Traj> > sets = genTrajSet (conds, d_conds, time_manouver, speed_goal, lane_desired, {speed_limit,10.,10.}, max_min);
                      
                      // set cost for every trajectory and combine them. Check also for collisions and dynamic limits
                      combSet = combineTrajectories(sets[0], sets[1], time_manouver, speed_goal, lane_desired, {speed_limit,10.,10.}, max_min, near_cars); 
                      
                      
                      // find minimal cost trajectory
                      if (combSet.size() > 0)
                      {
                        std::sort (combSet.begin(), combSet.end(), funcia);
                      }
                      else 
                      {
                        time_manouver *= 0.9; // if no trajectory is found, repeate trajectories generation with a smaller time horizon
                      }
                    }
                    
                    longTrajectory = combSet[0].Trs;
                    lateralTrajectory = combSet[0].Trd;
                    
                       

                    for (int i = 0; i < (size_horizon - size_keep); i++)
                    {
                      //generate next points using selected trajectory with a time pace of 0.02 seconds
                      next_s = longTrajectory.getDis(i*0.02);
                      next_d = lateralTrajectory.getDis(i*0.02);
                      
                      // convert  to  global coordinates
                      sxy = parabolicGetXY(next_s, next_d,  map_waypoints_s, map_waypoints_x, map_waypoints_y);
            
                      // pass points to simulator
                      next_x_vals.push_back(sxy[0]);
                      next_y_vals.push_back(sxy[1]);
                      
                      size_prev_path++;
                    }
                    
                    size_prev_plan = 0.;

                  } else {
                    // NO PLANNING BECAUSE LAST PATH IS NOT EXPIRED 
                    for(int i = 0; i < path_size; i++)
                    {
                      next_x_vals.push_back(previous_path_x[i]);
                      next_y_vals.push_back(previous_path_y[i]);
                    }
                  }
              }
              
             
          	
          	//END TODO
            msgJson["next_x"] = next_x_vals;
            msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
