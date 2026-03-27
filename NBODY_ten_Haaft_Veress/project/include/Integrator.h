#pragma once
using namespace std;
#include "body.h"
#include <vector>
#include "vector.h"
class Integrator {
public:
	//TODO:
	//takes present bodies as a List? and computes next positions and velocities
	virtual void integrate(std::vector<Body> bodies) const = 0;

	virtual Customvectors::Vector calculateAcceleration();

};





