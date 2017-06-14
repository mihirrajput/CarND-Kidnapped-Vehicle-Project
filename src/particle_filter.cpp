/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <map>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	// Initialize the number of particles
	num_particles = 100;
	// Make a Generator
	default_random_engine gen;
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std[1]);
	// This line creates a normal (Gaussian) distribution for theta
	normal_distribution<double> dist_theta(theta, std[2]);
	// Sample particles from the distribution
	for (int i = 0; i < num_particles; ++i) {
		// Sample from the normal distrubtions 
		// where "gen" is the random engine initialized earlier.
		Particle particle;
		particle.id = i;
		particle.x = dist_x(gen);
		particle.y = dist_y(gen);
		particle.theta = dist_theta(gen);
		particle.weight = 1.0;
		particles.push_back(particle);
		weights.push_back(1.0);

		// debug
		cout << particles[i].x << "\t" << particles[i].y << "\t" << particles[i].theta << "\n";
	}
	// To make sure init executes only once
	is_initialized = true;
	// debug
	cout << "init complete" << "\n";
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	default_random_engine gen;
	normal_distribution<double> dist_x(0, std_pos[0]);
	normal_distribution<double> dist_y(0, std_pos[1]);
	normal_distribution<double> dist_theta(0, std_pos[2]);

	double yawrate_times_dt = yaw_rate*delta_t;
	cout << "72"<<"\n";

	for (auto& particle:particles) {
		// Sample from the normal distrubtions 
		// where "gen" is the random engine initialized earlier.
		double theta = particle.theta;
		cout << "78" << "\n";

		// expected pose value from deterministic motion model
		if (yaw_rate>0.001)
		{
			double ratio_vel_yawrate = velocity / yaw_rate;
			particle.x = particle.x + ratio_vel_yawrate*(sin(theta + yawrate_times_dt) - sin(theta)) + dist_x(gen);
			particle.y = particle.y + ratio_vel_yawrate*(-cos(theta + yawrate_times_dt) + cos(theta)) + dist_y(gen);
			cout << "86" << "\n";
		}
		else
		{
			particle.x = particle.x + velocity*cos(theta)*delta_t + dist_x(gen);
			particle.y = particle.y + velocity*sin(theta)*delta_t + dist_y(gen);
			cout << "92" << "\n";
		}
		particle.theta = theta + yawrate_times_dt + dist_theta(gen);
		cout << "95" << "\n";

		// debug
		cout << particle.x << "\t" << particle.y << "\t" << particle.theta << "\n";
	}
	// debug
	cout << "Predict complete" << "\n";
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	for (int i = 0; i < observations.size(); ++i)
	{
		double min = 999.9;
		int min_id = 0;
		for (int j = 0; j < predicted.size(); ++j)
		{
			double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
			if ( distance < min) {
				min = distance;
				min_id = predicted[j].id;
			}
		}
		observations[i].id = min_id;
	}
}

double ParticleFilter::prob(double x, double y, float xm, float ym, double std_landmark[]) {
	double xpart = -0.5*(x - xm)*(x - xm) / (std_landmark[0] * std_landmark[0]);
	double ypart = -0.5*(y - ym)*(y - ym) / (std_landmark[1] * std_landmark[1]);
	return exp(xpart + ypart) / (2 * M_PI*std_landmark[0] * std_landmark[1]);
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	Map temp_landmarks;
	weights.clear(); // clear the previous weights vector

	for (auto& particle : particles)
	{
		double ctheta = cos(particle.theta);
		double stheta = sin(particle.theta);
		double xtrans = particle.x;
		double ytrans = particle.y;

		for (int z = 0; z < map_landmarks.landmark_list.size(); ++z)
		{
			if (dist(xtrans, ytrans, map_landmarks.landmark_list[z].x_f, map_landmarks.landmark_list[z].y_f) < sensor_range) {
				temp_landmarks.landmark_list.push_back(map_landmarks.landmark_list[z]);
			}
		}
		cout << "pf cpp line 159" << "\n";

		LandmarkObs LMObs;
		int j = 0;
		for (auto& LMObs : observations)
		{
			double temp_x = ctheta*LMObs.x - stheta*LMObs.y + xtrans;
			double temp_y = stheta*LMObs.x + ctheta*LMObs.y + ytrans;
			particle.sense_x[j] = temp_x; //obs rotated then translated 
			particle.sense_y[j] = temp_y; //obs rotated then translated

			double min = 999.9;
			int min_id = 0; 
			cout << "pf cpp line 171" << "\n";
			for (int k = 0; k < temp_landmarks.landmark_list.size(); ++k)
			{
				double distance = dist(temp_x, temp_y, temp_landmarks.landmark_list[k].x_f, temp_landmarks.landmark_list[k].y_f);
				if (distance < min)
				{
					min = distance;
					min_id = temp_landmarks.landmark_list[k].id_i;
				}
			}
			cout << "pf cpp line 179" << "\n";
			particle.associations[j] = min_id;
			particle.weight *= prob(temp_x, temp_y, map_landmarks.landmark_list[min_id].x_f, map_landmarks.landmark_list[min_id].y_f, std_landmark);
			j += 1;
		}
		cout << "pf cpp line 183" << "\n";
		weights.push_back(particle.weight);
		temp_landmarks.landmark_list.clear();
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	vector<Particle> new_particles;
	std::random_device rd;
	std::mt19937 gen(rd());
	std::discrete_distribution<> d(weights.begin(), weights.end());
	std::map<int, int> m;
	for (int n = 0; n<num_particles; ++n) {
		new_particles.push_back(particles[d(gen)]);
	}
	particles = new_particles; // resampled set of particles
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
