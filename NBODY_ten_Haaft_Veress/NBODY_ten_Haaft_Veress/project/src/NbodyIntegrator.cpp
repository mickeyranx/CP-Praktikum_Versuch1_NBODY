#include "NbodyIntegrator.h"
#include <iostream>
#include <fstream>

using namespace Customvectors;



NbodyIntegrator::NbodyIntegrator(TimeStep timeStepFunction)
{
	tsf = timeStepFunction;

}

NbodyIntegrator::NbodyIntegrator(TimeStep time_step_function, double max_step)
{
	tsf = time_step_function;
	maximum_timestep = max_step;
}


//this function will be overriden by the child classes and will never be used
void NbodyIntegrator::startIntegration(std::vector<Body> bodies, double timeStep, double max_integration_time, std::string outputPath)
{
	
	std::cout << "test" << std::endl;

}

Customvectors::Vector NbodyIntegrator::calculateAcceleration(std::vector<Body> &bodies, Body &reference_body, int N, int j) {
	Customvectors::Vector a = Customvectors::Vector(0,0,0);

	//iterate over all bodies
	for (int i = 0; i < N; i++)
	{
		//no self-interaction 
		if (i == j) {
			continue;
		}
		else {
			Body current_body = bodies[i];
			Customvectors::Vector r_ij = current_body.getPosition() - reference_body.getPosition();
			//calculate acceleration on body j
			a = a + r_ij * (1 / pow(r_ij.getLength(), 3)) * current_body.getMass();

		}



	}
	return a;
}



Customvectors::Vector NbodyIntegrator::calculateJerk(std::vector<Body> &bodies, Body &reference_body, int N, int j) {
	Customvectors::Vector a_dot = Customvectors::Vector(0, 0, 0);
	//iterate over all bodies
	for (int i = 0; i < N; i++)
	{
		//no self-interaction 
		if (i == j) {
			continue;
		}
		else {
			Body current_body = bodies[i];
			Customvectors::Vector r_ij = current_body.getPosition() - reference_body.getPosition();
			Customvectors::Vector v_ij = current_body.getVelocity() - reference_body.getVelocity();
			double r_len = r_ij.getLength();

			//calculate acceleration on body j
			a_dot = a_dot + ((v_ij * pow((1 / r_len), 3)) - (r_ij * (1 / pow(r_len, 5)) * (r_ij * v_ij) * 3)) * current_body.getMass();
		}



	}
	return a_dot;
}


double NbodyIntegrator::calculateEnergy(std::vector<Body> &bodies, int N) {
	
	double E = 0;
	for (int i = 0; i < N; i++)
	{
		Body b1 = bodies[i];
		double m1 = b1.getMass();
		Vector v1 = b1.getVelocity();
		
		E += 0.5 * m1 * (v1 * v1); //kinetic part
		for (int j = i + 1; j < N; j++) 
		{
			Body b2 = bodies[j];
			E += -1 * m1 * b2.getMass() / (b1.getPosition() - b2.getPosition()).getLength(); //potential part, G is set to 1
			
		}



	}

	return E;

}

Customvectors::Vector NbodyIntegrator::calculateMomentum(std::vector<Body>& bodies, int N)
{
	Customvectors::Vector p = Customvectors::Vector(0, 0, 0);
	for (int i = 0; i < N; i++)
	{
		Body b = bodies[i];
		p = p + b.getVelocity() * b.getMass();

	}
	return p;
}

Customvectors::Vector NbodyIntegrator::calculateAngularMomentum(std::vector<Body>& bodies, int N)
{
	Customvectors::Vector l = Customvectors::Vector(0, 0, 0);
	for (int i = 0; i < N; i++)
	{
		Body b = bodies[i];
		// % is vector product
		l = l + b.getPosition() % b.getVelocity() * b.getMass();

	}
	return l;
}


Customvectors::Vector NbodyIntegrator::calculateSpecificAngularMomentum(std::vector<Body>& bodies) {
	Body b1 = bodies[0];
	Body b2 = bodies[1];
	double m1 = b1.getMass();
	double m2 = b2.getMass();
	Vector r1 = b1.getPosition();
	Vector r2 = b2.getPosition();
	Vector v1 = b1.getVelocity();
	Vector v2 = b2.getVelocity();
	//reduced mass
	//double mu = m1 * m2 / (m1 + m2);
	//relative position and velocity
	Vector r = r2 - r1;
	Vector v = v2 - v1;
	// % is vector product
	return (r % v);

}


Customvectors::Vector NbodyIntegrator::calulateRungeLenz(std::vector<Body>& bodies, Customvectors::Vector j) {
	Body b1 = bodies[0];
	Body b2 = bodies[1];
	double m1 = b1.getMass();
	double m2 = b2.getMass();
	Vector r1 = b1.getPosition();
	Vector r2 = b2.getPosition();
	Vector v1 = b1.getVelocity();
	Vector v2 = b2.getVelocity();
	Vector r = r2 - r1;
	//Runge Lenz vector
	Vector e = ((v2 - v1) % j) - r * (1/ r.getLength());
	return e;
}

double NbodyIntegrator::calculateMajorSemiAxis(Customvectors::Vector j, Customvectors::Vector e) {

	return (j * j) / (1 - e * e);

}

double NbodyIntegrator::timeStepCurvature(std::vector<Body> &bodies, std::vector<Customvectors::Vector> &accelerations ,int N,double eta) {
	
	std::vector<double> vals = {};
	for (int i = 0; i < N; i++)
	{
		Body b = bodies[i];
		//Customvectors::Vector acc = calculateAcceleration(bodies, b, N, i);
		Customvectors::Vector jerk = calculateJerk(bodies, b, N, i);
		double jerk_norm = jerk.getLength();
		if (jerk_norm == 0.0) continue;
		vals.push_back(accelerations[i].getLength() / jerk.getLength());
	}

	if (vals.empty()) return -1.0; //if the list of jerks is empty return -1.0 
	
	return eta * *std::min_element(vals.begin(), vals.end());


}

