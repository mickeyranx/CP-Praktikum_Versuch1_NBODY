#pragma once
#include <iostream>
#include <fstream>
#include <tuple>

#include "NbodyIntegrator.h"


class Heun : public NbodyIntegrator {
public:

	Heun(TimeStep tsf) : NbodyIntegrator(tsf) {};
	Heun(TimeStep tsf, double max_step) : NbodyIntegrator(tsf, max_step) {};

	std::tuple<std::vector<Body>, double> integrate(std::vector<Body> &image, double eta, int N);

	void startIntegration(std::vector<Body> inital_image, double eta, double max_integration_time, std::string output_file) override;

};