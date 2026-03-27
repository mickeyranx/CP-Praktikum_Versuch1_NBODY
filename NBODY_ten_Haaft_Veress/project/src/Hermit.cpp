#include "Hermit.h"





std::tuple<std::vector<Body>, double> Hermit::integrate(std::vector<Body> &image,  double eta, int N)
{
	std::vector<Customvectors::Vector> a_ns = {};
	//------------------------------------------------
	//          calculate accelerations
	//------------------------------------------------
	for (int i = 0; i < N; i++)
	{
		a_ns.push_back(calculateAcceleration(image, image[i], N, i));
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
	case TimeStep::CURVATUREHERMIT:
		double new_time_step = timeStepCurvature(image, a_ns, N, eta);
		double max_step = getMaxTimeStep();
		dt = (new_time_step > max_step || new_time_step < 0) ? max_step : new_time_step;
		break;
	}

	//------------------------------------------------
	//          integration
	//------------------------------------------------
	std::vector<Body> new_temp_image = {};
	
	std::vector<Customvectors::Vector> a_n_dots = {};
	std::vector<Customvectors::Vector> a_3s = {};
	std::vector<Customvectors::Vector> a_2s = {};
	for (int i = 0; i < N; i++)
	{
		Body current_body = image[i];
		Customvectors::Vector a_n = a_ns[i];
		Customvectors::Vector a_n_dot = calculateJerk(image, current_body, N, i);
		Customvectors::Vector v_n = current_body.getVelocity();
		Customvectors::Vector r_n = current_body.getPosition();
		Customvectors::Vector r_nn_p = r_n + v_n * dt + a_n * 0.5 * dt * dt + a_n_dot * (1.0/6.0) * pow(dt, 3);
		Customvectors::Vector v_nn_p = v_n + a_n * dt + a_n_dot * 0.5 * dt * dt;
		current_body.setPosition(r_nn_p);
		current_body.setVelocity(v_nn_p);
		a_n_dots.push_back(a_n_dot);
		new_temp_image.push_back(current_body);

	}
	std::vector<Body> new_image = {};
	for (int i = 0; i < N; i++)
	{
		Body current_body = new_temp_image[i];
		Customvectors::Vector a_nn_p = calculateAcceleration(new_temp_image, current_body, N, i);
		Customvectors::Vector a_nn_p_dot = calculateJerk(new_temp_image, current_body, N, i);
		Customvectors::Vector a_n = a_ns[i];
		Customvectors::Vector a_n_dot = a_n_dots[i];

		Customvectors::Vector r_nn_p = current_body.getPosition();
		Customvectors::Vector v_nn_p = current_body.getVelocity();
		//Hermit-Interpolation
		Customvectors::Vector a_2 = (a_n - a_nn_p) * (-3.0 / (dt * dt)) - (a_n_dot * 2.0 + a_nn_p_dot) * (1/dt);
		Customvectors::Vector a_3 = (a_n - a_nn_p) * (2.0 / (dt * dt * dt)) + (a_n_dot + a_nn_p_dot) * (1 / (dt * dt));
		a_2s.push_back(a_2); //with this we can consider higher order to calculate the next timestep dynamically
		a_3s.push_back(a_3); 
		//correction

		/*-------------------------------------------------------------------------------
		if the time step is calculated dynamically with higher orders of the acceleration
		a_2 and a_3 will be calculated with a standart dynamical time step 
		-------------------------------------------------------------------------------*/

		if (!(getTimeStep() == TimeStep::CURVATUREHERMIT)) {
			Customvectors::Vector v_nn_c = v_nn_p + a_2 * (1.0 / 3.0) * pow(dt, 3) + a_3 * (1.0 / 4.0) * pow(dt, 4);
			Customvectors::Vector r_nn_c = r_nn_p + a_2 * (1.0 / 12.0) * pow(dt, 4) + a_3 * (1.0 / 20.0) * pow(dt, 5);
			current_body.setVelocity(v_nn_c);
			current_body.setPosition(r_nn_c);
			new_image.push_back(current_body);


		}
	}
	if (getTimeStep() == TimeStep::CURVATUREHERMIT) {
		double new_time_step = calculateTimeStepHermit(a_ns, a_n_dots, a_2s, a_3s, N, eta);
		double max_step = getMaxTimeStep();
		dt = (new_time_step > max_step || new_time_step < 0) ? max_step : new_time_step;
		for (int i = 0; i < N; i++)
		{
			Body current_body = new_temp_image[i];
			Customvectors::Vector r_nn_p = current_body.getPosition();
			Customvectors::Vector v_nn_p = current_body.getVelocity();
			Customvectors::Vector a_2 = a_2s[i];
			Customvectors::Vector a_3 = a_3s[i];
			Customvectors::Vector v_nn_c = v_nn_p + a_2 * (1.0 / 3.0) * pow(dt, 3) + a_3 * (1.0 / 4.0) * pow(dt, 4);
			Customvectors::Vector r_nn_c = r_nn_p + a_2 * (1.0 / 12.0) * pow(dt, 4) + a_3 * (1.0 / 20.0) * pow(dt, 5);
			current_body.setVelocity(v_nn_c); //this doesnt break symmetry of the forces, because "current_body" isnt saved for the next iteration of the for loop
			current_body.setPosition(r_nn_c);
			new_image.push_back(current_body);

		}



	}


	

	return std::make_tuple(new_image, dt); //this line of code should never run
}

void Hermit::startIntegration(std::vector<Body> initial_image, double eta, double max_integration_time, std::string output_file)
{
	int N = initial_image.size();
	double time_passed = 0;
	double time_step;
	//-------------------------------------
	//        file setup
	//-------------------------------------
	std::ofstream File(output_file);
	File << "#this data set is generated with the hermit integrator with eta = " + std::to_string(eta) + " and t_max = " + std::to_string(max_integration_time) + "\n";
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

double Hermit::calculateTimeStepHermit(
	std::vector<Customvectors::Vector> &as, 
	std::vector<Customvectors::Vector> &a_dots, 
	std::vector<Customvectors::Vector> &a_2s, 
	std::vector<Customvectors::Vector> &a_3s, int N, double dt)
{
	std::vector<double> acc(N, 0);
	for (int i = 0; i < N; i++)
	{
		double a_dot_i = a_dots[i].getLength();
		double a_2 = a_2s[i].getLength();
		acc[i] = sqrt((as[i].getLength() * a_2 + a_dot_i * a_dot_i) / (a_dot_i * a_3s[i].getLength() + a_2 * a_2));
	}

	return dt * *std::min_element(acc.begin(), acc.end());
}
