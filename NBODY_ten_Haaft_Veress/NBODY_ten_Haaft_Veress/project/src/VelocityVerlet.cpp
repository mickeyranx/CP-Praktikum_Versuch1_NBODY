#include "VelocityVerlet.h"




std::tuple<std::vector<Body>, double> VelocityVerlet::integrate(std::vector<Body> &image,  double eta, int N)
{
	
	
	
	std::vector<Customvectors::Vector> a_ns = {};
	//------------------------------------------------
	//          calculate accelerations
	//------------------------------------------------
	for (int i = 0; i < N; i++)
	{
		a_ns.push_back(calculateAcceleration(image ,image[i], N, i));
	}
	//------------------------------------------------
	//          determine time step
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
		double new_time_step = timeStepCurvature(image, a_ns ,N, eta);
		double max_step = getMaxTimeStep();
		dt = (new_time_step > max_step || new_time_step < 0) ? max_step : new_time_step;
		break;

	}
	std::vector<Body> new_temp_image = {};
	//-------------------------------------------------
	//                  integration
	//-------------------------------------------------

	//kick-drift for dynamic time step is not yet available

	for (int i = 0; i < N; i++)
	{
		Body current_body = image[i];
		Customvectors::Vector v_n = current_body.getVelocity();
		Customvectors::Vector r_n = current_body.getPosition();

		Customvectors::Vector r_nn = r_n + v_n * dt + a_ns[i] * 0.5 * dt * dt;
		current_body.setPosition(r_nn);
		new_temp_image.push_back(current_body);
		
	}
	std::vector<Body> new_image = {};
	for (int i = 0; i < N; i++)
	{
		Body current_body = new_temp_image[i];
		Customvectors::Vector a_nn = calculateAcceleration(new_temp_image, current_body, N, i);
		Customvectors::Vector v_n = current_body.getVelocity();
		Customvectors::Vector v_nn = v_n + (a_ns[i] + a_nn) * 0.5 * dt;
		current_body.setVelocity(v_nn);
		new_image.push_back(current_body);
	}


	return std::make_tuple(new_image, dt);
}

void VelocityVerlet::startIntegration(std::vector<Body> initial_image, double eta, double max_integration_time, std::string output_file)
{
	int N = initial_image.size();
	double time_passed = 0;
	double time_step;
	//-------------------------------------
	//        file setup
	//-------------------------------------
	std::ofstream File(output_file);
	File << "#this data set is generated with the velocity-verlet integrator with eta = " + std::to_string(eta) + " and t_max = " + std::to_string(max_integration_time) + "\n";
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
	//-------------------------------------
	//            integration
	//-------------------------------------
	std::vector<Body> previous_image = initial_image;
	while (time_passed < max_integration_time) {
		std::vector<Body> new_image;


		std::tie(new_image, time_step) = integrate(previous_image, eta, N);
		time_passed += time_step;
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


		
		previous_image = new_image; //prepare for next iteration

	}


	File.close();

}
