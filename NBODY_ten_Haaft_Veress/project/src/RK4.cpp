#include "RK4.h"
#include <iostream>


std::tuple<std::vector<Body>, double> RK4::integrate(std::vector<Body> &image,  double eta, int N)
{
	std::vector<Customvectors::Vector> a_ns; //accelerations
	//------------------------------------------------
	//         calculate accelerations
	//------------------------------------------------
	for (int i = 0; i < N; i++)
	{
		a_ns.push_back(calculateAcceleration(image, image[i], N, i));
	}

	//------------------------------------------------
	//          determine next time step
	//------------------------------------------------
	double dt = 0;
	switch (getTimeStep())
	{
	case TimeStep::LINEAR:
		dt = eta;
		break;
	case TimeStep::QUADRATIC:
		dt = eta * eta;
		break;
	case TimeStep::DYNAMIC:
		double new_time_step = timeStepCurvature(image, a_ns,N, eta);
		double max_step = getMaxTimeStep();
		dt = (new_time_step > max_step || new_time_step < 0) ? max_step : new_time_step;
		break;

	}


	//------------------------------------------------
	//                 integration
	//------------------------------------------------
	std::vector<Body> new_image = image;
	std::vector<Body> new_temp_image_1 = {};
	std::vector<Body> new_temp_image_2 = {};

	for (int i = 0; i < N; i++)
	{
		Body current_body = image[i];
		
		Customvectors::Vector v_n = current_body.getVelocity();
		Customvectors::Vector r_n = current_body.getPosition();

		Customvectors::Vector v_1 = a_ns[i] * dt;
		Customvectors::Vector r_1 = v_n * dt;
		new_image[i].alterPosition(r_1*(1.0 / 6.0)); //add first order term
		new_image[i].alterVelocity(v_1*(1.0 / 6.0)); //add first order term
		current_body.setPosition(r_n + r_1 * 0.5);
		new_temp_image_1.push_back(current_body); //for calculation of v_2
		Customvectors::Vector r_2 = (v_n + v_1 * 0.5) * dt;
		new_image[i].alterPosition(r_2 * (1.0 / 3.0)); //add second order term
		current_body.setPosition(r_n + r_2 * 0.5);
		new_temp_image_2.push_back(current_body); //for calculation of v_3
		
		

	}
	std::vector<Body> new_temp_image_3 = {};
	for (int i = 0; i < N; i++)
	{
		Body current_body_1 = new_temp_image_1[i];
		Body old_body = image[i];
		Customvectors::Vector v_n = old_body.getVelocity();
		Customvectors::Vector r_n = old_body.getPosition();
		Customvectors::Vector v_2 = calculateAcceleration(new_temp_image_1, current_body_1, N, i) * dt;
		new_image[i].alterVelocity(v_2 * (1.0 / 3.0)); //add second order term
		Customvectors::Vector r_3 = (v_n + v_2 * 0.5) * dt;
		new_image[i].alterPosition(r_3 * (1.0 / 3.0)); //add third order term
		Body current_body_2 = new_temp_image_2[i];
		Customvectors::Vector v_3 = calculateAcceleration(new_temp_image_2, current_body_2, N, i) * dt;
		new_image[i].alterVelocity(v_3 * (1.0 / 3.0)); //add third order term
		Customvectors::Vector r_4 = (v_n + v_3) * dt;
		new_image[i].alterPosition(r_4 * (1.0 / 6.0)); //add fourth order term
		old_body.setPosition(r_n + r_3);
		new_temp_image_3.push_back(old_body);
		
	}

	for (int i = 0; i < N; i++)
	{
		Customvectors::Vector v_4 = calculateAcceleration(new_temp_image_3, new_temp_image_3[i], N, i) * dt;
		new_image[i].alterVelocity(v_4 * (1.0 / 6.0)); //add fourth order term
	}
	

	return std::make_tuple(new_image, dt);
	
}

void RK4::startIntegration(std::vector<Body> initial_image, double eta, double max_integration_time, std::string output_file)
{
	int N = initial_image.size();
	double time_passed = 0;
	double time_step;
	//-------------------------------------
	//        file setup
	//-------------------------------------
	
	//to test CPU time of one integration comment from here to ......
	
	std::ofstream File(output_file);
	File << "#this data set is generated with the RK4 integrator with eta = " + std::to_string(eta) + " and t_max = " + std::to_string(max_integration_time) + "\n";
	File << "t" << "\t";
	
	for (int i = 0; i < initial_image.size(); i++) //output file setup
	{
		File << "x_" + std::to_string(i) << "\t" << "y_" + std::to_string(i) << "\t" << "z_" + std::to_string(i) << "\t";
	}
	File << "E" << "\t" << "p_x" << "\t" << "p_y" << "\t" << "p_x" << "\t" << "p_z" << "\t" << "l_x" << "\t" << "l_y" << "\t" "l_z";
	if (N == 2) {
		File << "\t" << "|j|" << "\t" << "|e|" << "\t" << "a_maj";
	}
	File << "\n";

	File << std::fixed << std::setprecision(output_precision);
	//-------------------------------------
	//          inital image
	// ------------------------------------
	
	File << time_passed << "\t";
	for (Body b : initial_image) {
		File << b.getPosition().getX() << "\t" << b.getPosition().getY() << "\t" << b.getPosition().getZ() << "\t";

	}
	double E = calculateEnergy(initial_image, N);
	Customvectors::Vector p = calculateMomentum(initial_image, N);
	Customvectors::Vector l = calculateAngularMomentum(initial_image, N);
	File << E << "\t" << p.getX() << "\t" << p.getY() << "\t" << p.getZ() << "\t" << l.getX() << "\t" << l.getY() << "\t" << l.getZ() << "\t";
	if (N == 2) {
		Vector j = calculateSpecificAngularMomentum(initial_image);
		Vector e = calulateRungeLenz(initial_image, j);
		double a = calculateMajorSemiAxis(j, e);
		File << j.getLength() << "\t" << e.getLength() << "\t" << a;
	}
	
	File << "\n";
	
	//... here disable writing to file
	
	//-------------------------------------
	//            integration
	//-------------------------------------
	std::vector<Body> previous_image = initial_image;
	while (time_passed < max_integration_time) {
		std::vector<Body> new_image;

		std::tie(new_image, time_step) = integrate(previous_image, eta, N);
		
		
		time_passed += time_step;
		//comment from here to ...
		
		File << time_passed << "\t";
		
		for (Body b : new_image) {
			File << b.getPosition().getX() << "\t" << b.getPosition().getY() << "\t" << b.getPosition().getZ() << "\t";

		}
		
		
		//-------------------------------------
		//        conserved quantities
		//-------------------------------------
		double E = calculateEnergy(new_image, N);
		Customvectors::Vector p = calculateMomentum(new_image, N);
		Customvectors::Vector l = calculateAngularMomentum(new_image, N);
		File << E << "\t" << p.getX() << "\t" << p.getY() << "\t" << p.getZ() << "\t" << l.getX() << "\t" << l.getY() << "\t" << l.getZ() << "\t";

		if (N == 2) {
			Vector j = calculateSpecificAngularMomentum(new_image);
			Vector e = calulateRungeLenz(new_image, j);
			double a = calculateMajorSemiAxis(j, e);
			File << j.getLength() << "\t" << e.getLength() << "\t" << a;
		}
		File << "\n";
		
		//here ... to disable writing to file


		
		previous_image = new_image; //prepare for next iteration
	}
	
	//File.close();

}
