#pragma once
#include "NbodyIntegrator.h"
#include <iostream>
#include <tuple>
#include <fstream>

class EulerChromer : public NbodyIntegrator {
public:
	
	EulerChromer(TimeStep tsf) : NbodyIntegrator(tsf) {};
	EulerChromer(TimeStep tsf, double max_step) : NbodyIntegrator(tsf, max_step) {};
	std::tuple<std::vector<Body>, double> integrate(std::vector<Body> &image, double eta , int N);
	void startIntegration(std::vector<Body> inital_image, double eta, double max_integration_time, std::string output_file);


};


