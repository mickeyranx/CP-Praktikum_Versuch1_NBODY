#pragma once
#include "vector.h"
#include "TimeStep.h"
#include "Body.h"
#include <vector>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>

class NbodyIntegrator{


private:
	
	//type of parameter
	TimeStep tsf;
	double maximum_timestep = 0;

public:
	NbodyIntegrator(TimeStep time_step_function);

	//overloaded constructor to set minimum time step
	NbodyIntegrator(TimeStep time_step_function, double max_step);
	//calculates the next picture

	
	//starts the Integration with 
	virtual void startIntegration(std::vector<Body> bodies, double timeStep, double max_integration_time, std::string outputPath);
	

	//virtual void simulate(int Iterations, int N);
	//List of bodies, N = number of bodies, j ist the current body
	static Customvectors::Vector calculateAcceleration(std::vector<Body> &bodies, Body &reference_body, int N, int j);
	//List of bodies, N = number of bodies, j ist the current body
	static Customvectors::Vector calculateJerk(std::vector<Body> &bodies, Body &reference_body, int N, int j);

	TimeStep getTimeStep() { return tsf; }


	//------------------------------------------------------
	//			time step
	//------------------------------------------------------
	double getMaxTimeStep() {
		return this->maximum_timestep;
	}
			
	static const int output_precision = 17;

	static double timeStepLinear(double t) { return t; };

	static double timeStepQuadratic(double t) { return pow(t, 2); };

	double timeStepCurvature(std::vector<Body>& bodies, std::vector<Customvectors::Vector> &accelerations, int N ,double eta);

	//------------------------------------------------------
	//			conserved quantities
	//------------------------------------------------------
	//calculates the Energy of the current state, N is the number of bodies present
	static double calculateEnergy(std::vector<Body>& bodies, int N);

	//calculates the total momentum (vector) of a system with N bodies
	Customvectors::Vector calculateMomentum(std::vector<Body>& bodies, int N);

	//calculates the angular momentum (vector) of a system with N bodies
	Customvectors::Vector calculateAngularMomentum(std::vector<Body>& bodies, int N);

	//2 Body:
	//calculates the angulare momentum of the 2 Body system
	Customvectors::Vector calculateSpecificAngularMomentum(std::vector<Body>& bodies);

	//calculates the Runge Lenz vector of the 2 Body system
	Customvectors::Vector calulateRungeLenz(std::vector<Body>& bodies, Customvectors::Vector j);

	//calculate the MajorSemiAxis of the orbit with the Runge-Lenz vector e and the angular momentum j
	double calculateMajorSemiAxis(Customvectors::Vector j, Customvectors::Vector e);

	


};

