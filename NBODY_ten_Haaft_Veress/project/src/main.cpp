
#include <iostream>
#include "vector.h"
#include <fstream>
#include "Body.h"
#include <vector>
#include <string>
#include "NbodyIntegrator.h"
#include <sstream>
#include <iomanip>
#include <format>
// integrators -------------------------
#include "Euler.h"
#include "Hermit.h"
#include "IteratedHermit.h"
#include "EulerChromer.h"
#include "Heun.h"
#include "RK4.h"
#include "VelocityVerlet.h"
//--------------------------------------
#include "TimeStep.h"
#include <tuple>


using namespace std;
using namespace Customvectors;





//this function reads a file containing the initial values of the bodies in the x y z vx vy vz m format
//             and parses the COM-transformed data into a list of bodies
static std::vector<Body> setup(std::string file_name) {
	//--------------------------------------------------
	//          parse starting data from file 
	//               into a list of bodies
	//--------------------------------------------------

	cout << "starting the setup with the input data" << endl;

	//open file with given file name
	ifstream file(file_name);
	
	//initialize list  
	std::vector<Body> bodies;
	
	double x;
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double m;
	//read doubles from file until there are no more values left

	for (int i = 0; (file >> x >> y >> z >> vx >> vy >> vz >> m); i++)
	{
		Vector pos = Vector(x, y, z);
		Vector vel = Vector(vx, vy, vz);
		Body body = Body(m, pos, vel);
		//add created body to the list
		bodies.push_back(body);

	}
	

	
	//--------------------------------------------------
	//          calculate COM reference
	//--------------------------------------------------
	Vector r_s = Vector(0, 0, 0);
	Vector v_s = Vector(0, 0, 0);
	double M = 0;
	for (int i = 0; i < (int) bodies.size(); i++)
	{
		Body b = bodies[i];
		double m = b.getMass();
		r_s = r_s + b.getPosition() * m;
		v_s = v_s + b.getVelocity() * m;
		M += m; //should add up to 1 in the end
	}
	cout << "checking the total mass: M = " << M << endl;

	vector<Body> bodies_COM = {};
	for (int i = 0; i < (int) bodies.size(); i++)
	{
		Body b = bodies[i];
		Vector vel_COM =  (b.getVelocity() - v_s) * (1 / M);
		Vector pos_COM = (b.getPosition() - r_s) * (1 / M);
		bodies_COM.push_back(Body(b.getMass(), pos_COM, vel_COM));
	}
	
	cout << "setup done" << endl;
	return bodies_COM;

}

//main function 
int main() {
	clock_t start = clock();


	//------------------------------------------------
	//       select integrator and time step
	//				  calculation
	//------------------------------------------------
	/*
	   available Time steps:
			-LINEAR
			-QUADRATIC
			-CURVATURE
			-HERMIT (exclusive for hermit and iterated hermit)

		available Integrators:
			-Euler
			-EulerChromer
			-VelocityVerlet
			-Hermit
			-IteratedHermit
			-Heun
			-RK4
		examples:
		     RK4 integrator = RK4(TimeStep::LINEAR);
			Euler integrator = Euler(TimeStep::DYNAMIC, 0.04); (overloaded constructor with TimeStep and max time_step)

	*/

	
	

	//create integrator objects with linear Timesteps
	Euler euler_integrator = Euler(TimeStep::LINEAR);
	EulerChromer ec_integrator = EulerChromer(TimeStep::LINEAR);
	Heun heun_integrator = Heun(TimeStep::LINEAR);
	Hermit hermit_integrator = Hermit(TimeStep::LINEAR);
	IteratedHermit it_hermit_integrator = IteratedHermit(TimeStep::LINEAR);
	VelocityVerlet vv_integrator = VelocityVerlet(TimeStep::LINEAR);
	RK4 rkfour_integrator = RK4(TimeStep::LINEAR);


	double eta;
	double t_max;

	//------------------------------------------------
	//          exercise 2 : 2-Body Problem
	//------------------------------------------------

	// !!!output precision can be changed in "NbodyIntegrator.h", default is 17!!!
	 
	//----------------------------------------------------------------------------------------
	//    exercise 2, part 1: try out different timesteps and investigate the different 
	//                                 integrators
	//----------------------------------------------------------------------------------------

	

	//list of constant timesteps
	std::vector<double> etas = {0.5, 0.1, 0.05, 0.01, 0.005, 0.001};

	//Setup and loading files:
	//inputfile-name for body data 
	string input_name = "2body.txt"; //filename of constellation input

	//filename of integration-info input
	string input_name2 = "in2.txt";

	string output_path = "output_files/2.1/"; //output directory path
	std::vector<Body> starting_image = setup("input_files/" + input_name);
	t_max = 1000.0; //since a_maj is 1 and T = 2pi sqrt(a_maj^3) = 2pi this will result in about 160 periosds
	

	
	//generate outputfiles for the different timesteps in the "output_files/2.1/" subdirectory
	for (double eta : etas){
		stringstream stream;
		stream << fixed << setprecision(3) << eta;
		string eta_as_string = stream.str();
		euler_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_euler_eta=" + eta_as_string + ".txt");
		//ec_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_eulerchromer_eta=" + eta_as_string + ".txt");
		//heun_integrator.startIntegration(starting_image, eta, t_max, output_path  + "2body_heun_eta=" + eta_as_string + ".txt");
		//hermit_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_hermit_eta=" + eta_as_string + ".txt");
		//it_hermit_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_iteratedhermit_eta=" + eta_as_string + ".txt");
		//vv_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_velocityverlet_eta=" + eta_as_string + ".txt");
		//rkfour_integrator.startIntegration(starting_image, eta, t_max, output_path + "2body_rk4_eta=" + eta_as_string + ".txt");

	}

	

	//-----------------------------------------------------------------------------
	//    exercise 2, part 2: alter starting velocity for different orbits
	//-----------------------------------------------------------------------------
	/*
	cout << "starting exercise 2 part 2" << endl;
	
	//eta = 0.005;
	eta = 0.001;
	t_max = 1000.0;

	stringstream stream;
	stream << fixed << setprecision(3) << eta;
	string eta_as_string = stream.str();
	string output_path = "output_files/2.2/";
	string input_name;
	std::vector<Body> starting_image;
	
	int variations_count = 7; //number of input files with different velocities
	
	
	for (int i = 1; i < variations_count + 1; i++)
	{
		cout << "doing simulations with variation " << i << endl;
		input_name = "2body_var" + std::to_string(i) + ".txt";
		starting_image = setup("input_files/" + input_name);
		euler_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_euler_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		ec_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_ec_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		hermit_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_hermit_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		heun_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_heun_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		it_hermit_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_ithermit_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		vv_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_vv_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");
		rkfour_integrator.startIntegration(starting_image, eta, t_max,
			output_path + "2body_rk4_var" + std::to_string(i) + "_eta = " + eta_as_string + ".txt");

	}
	

	cout << "exercise 2 part 2 done" << endl;
	*/
	//-----------------------------------------------------------------------------
	//    exercise 3:  100- and 1000-Body simulation
	//-----------------------------------------------------------------------------

	std::vector<Body> starting_image_2N = setup("input_files/2body.txt");
	std::vector<Body> starting_image_3N = setup("input_files/3body.txt");
	std::vector<Body> starting_image_100N = setup("input_files/100body.txt");
	std::vector<Body> starting_image_1000N = setup("input_files/1kbody.txt");

	//-----------------------------------------------------------------------------
	//             3.1: checking the CPU-Time for one integration step
	//-----------------------------------------------------------------------------
	

	//check that writing to file is disabled in rk4.cpp via commenting!
	/*
	eta = 0.01;
	t_max = 1.0; //100 integrations -> divide result by 100
	clock_t start_1 = clock();
	rkfour_integrator.startIntegration(starting_image_2N, eta, t_max, "output_files/testing_time_elapsed.txt");
	clock_t end_1 = clock();
	
	printf("\n time for one integration of 2-Body system: %.16f sec", (double(end_1 - start_1) / CLOCKS_PER_SEC) / 100.0);
	clock_t start_2 = clock();
	rkfour_integrator.startIntegration(starting_image_3N, eta, t_max, "output_files/testing_time_elapsed.txt");
	clock_t end_2 = clock();
	printf("\n time for one integration of 3-Body system: %.16f sec", (double(end_2 - start_2) / CLOCKS_PER_SEC) / 100.0);
	clock_t start_3 = clock();
	rkfour_integrator.startIntegration(starting_image_100N, eta, t_max, "output_files/testing_time_elapsed.txt");
	clock_t end_3 = clock();
	printf("\n time for one integration of 100-Body system: %.8f sec", (double(end_3 - start_3) / CLOCKS_PER_SEC) / 100.0);
	t_max = 0.2;
	clock_t start_4 = clock();
	rkfour_integrator.startIntegration(starting_image_1000N, eta, t_max, "output_files/testing_time_elapsed.txt");
	clock_t end_4 = clock();
	printf("\n time for one integration of 1000-Body system: %.8f sec", (double(end_4 - start_4) / CLOCKS_PER_SEC) / 20.0);
	*/


	//--------------------------------------------------------------------------------------------------------------
	//           3.2: running 3, 100 and 1000-Body simulations with rk4 with different timesteps
	//--------------------------------------------------------------------------------------------------------------
	
	
	/*
	
	//std::vector<double> etas = {0.1, 0.05, 0.025, 0.02, 0.0125, 0.01};
	std::vector<double> etas = {0.02};
	//std::vector<double> etas = {0.0075, 0.005, 0.006};
	t_max = 0.1;
	
	//std::vector<double> etas = {0.0001};
	 // tmax is 2.0 in all input files
	stringstream stream1;
	stream1 << fixed << setprecision(3) << t_max;
	string output_path_front = "output_files/3.0/rk4/";
	string output_path_back = "_tmax=" + stream1.str() + ".txt";

	for (double eta : etas)
	{
		
		stringstream stream2;
		stream2 << fixed << setprecision(5) << eta;
		string eta_string = stream2.str();
		cout << "\n now integrating with eta: " + eta_string;
		
		string output_path_100N = output_path_front + "100Body_rk4_eta=" + eta_string + output_path_back;
		rkfour_integrator.startIntegration(starting_image_100N, eta, t_max, output_path_100N);
		cout << "\n 100-body done";
		
		clock_t t0 = clock();
		string output_path_1000N = output_path_front + "1000Body_rk4_eta=" + eta_string + output_path_back;
		rkfour_integrator.startIntegration(starting_image_1000N, eta, t_max, output_path_1000N);
		cout << "\n 1000-body done";
		clock_t t1 = clock();
		cout << "time: " << double(t1 - t0) / CLOCKS_PER_SEC << endl;

	}
	
	
	*/
	//----------------------------------------------------------------------------------------------------------------
	clock_t end = clock();
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	printf("\n execution time of the entire program: %.3f sec", elapsed);
}



