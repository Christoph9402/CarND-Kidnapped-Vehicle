/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"
using namespace std;
using std::string;
using std::vector;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //Set to 100
  num_particles = 100;
  
  //Generator to create random numbers
  std::random_device rand; 
  std::mt19937 generator(rand());
  
  //use x, y and theta as mean values
  normal_distribution<double> x_noise(x, std[0]);
  normal_distribution<double> y_noise(y, std[1]);
  normal_distribution<double> t_noise(theta, std[2]);
  
  //Initialize all particles posisiton x and y, as well as theta
  //Loop to assign values to each particle
for(int i=0;i<num_particles;i++){
  	//Initialize Values for Structure of each particle
  	//Set weight for each particle to one
	particles[i].weight=1.0;
  	particles[i].id=i;
  	//Initialize position x, y and theta for each particle using a normal distribution
    particles[i].x=x_noise(generator);
    particles[i].y=y_noise(generator);
  	particles[i].theta=t_noise(generator);
  	particles.push_back(particles[i]);
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
//Generator to create random numbers
  std::random_device rand; 
  std::mt19937 generator(rand());
  
  //use 0 as mean values
  normal_distribution<double> x_noise(0, std_pos[0]);
  normal_distribution<double> y_noise(0, std_pos[1]);
  normal_distribution<double> t_noise(0, std_pos[2]);
  
  //Formulas for predicition (See Lesson 5.8) --> differentiate, if the yaw rate is zero or not
  
  for(int i=0;i<num_particles;i++){
  // If yaw rate is not 0:
  if(yaw_rate!=0.0){
  particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta+(yaw_rate * delta_t)) - sin(particles[i].theta));
  particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta)-cos(particles[i].theta+yaw_rate * delta_t));
    particles[i].theta += yaw_rate * delta_t;
  }
  //If yaw rate is 0, use formulas from video in lesson 3.4
  else {
    particles[i].x += velocity * delta_t * cos(particles[i].theta);
    particles[i].y += velocity * delta_t * sin(particles[i].theta);
    //Final yaw is equal to initial yaw --> add nothing
    particles[i].theta += 0;  
  }
  //add random gaussian noise
   particles[i].x += x_noise(generator);
   particles[i].y += y_noise(generator);
   particles[i].theta += t_noise(generator);

 }
  
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
	for (unsigned int i=0;i<observations.size();i++){
    //initialize the minimum distance
    	double min_dist = numeric_limits<double>::max();
    // initialize the elementnumber for the lowest distance as -1
    	int no_min_dist = -1;
    //while looping through the obserations vector, loop through the predicted vector
    	for (unsigned int j=0;j<predicted.size();j++){
    //calculate the distance
    		double distance = dist(observations[i].x,observations[i].y,predicted[j].x,predicted[j].y);
    //double distance = sqrt(pow(observations[i].x - predicted[j].x,2) + pow(observations[i].y - prediction[j].y, 2));
    //check, if the calulated current distance is smaller than the min_dist
    		if(distance<min_dist){
      //replace min_dist with the current distance
      			min_dist = distance;
      //replace no_min_dist with the id of the predicted element
      			no_min_dist = predicted[j].id;
    		}
    //give the element in the observations vector the same id as the element in the predicted vector, with which it shares the shortest distance
   		observations[i].id=no_min_dist;
    	}
  	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
    
for(int i=0;i<num_particles;i++){
  //create vector vor transformed points
  vector<LandmarkObs> Obs_transformed;
  vector<LandmarkObs> landmark_Range;
  particles[i].weight=1.0;
  
  //Transform Coordinates (formulas lesson5.17)
  double x_transformed, y_transformed;
  for(unsigned int j=0;j<observations.size();j++){
    
  	x_transformed=particles[i].x+cos(particles[i].theta)*observations[j].x-sin(particles[i].theta)*observations[j].y;
  y_transformed=particles[i].y+sin(particles[i].theta)*observations[j].x+cos(particles[i].theta)*observations[j].y;
    Obs_transformed.push_back(LandmarkObs{observations[j].id,x_transformed,y_transformed});                           
  }
  
  for (unsigned int j=0;j<map_landmarks.landmark_list.size();j++){
  	if(dist(map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f,particles[i].x,particles[i].y)<=sensor_range){
	landmark_Range.push_back(LandmarkObs{map_landmarks.landmark_list[j].id_i,map_landmarks.landmark_list[j].x_f,map_landmarks.landmark_list[j].y_f});
	} // end if
} //end j
  
dataAssociation(landmark_Range,Obs_transformed);
  
  //cout << landmark_Range.size() << endl;
  //cout << Obs_transformed.size() << endl;
for(unsigned int j=0;j<Obs_transformed.size();j++){
 x_transformed=Obs_transformed[j].x;
 y_transformed=Obs_transformed[j].y;
   
 for(unsigned int k=0;k<landmark_Range.size();k++){
   double x_landmark,y_landmark;
  if(landmark_Range[k].id==Obs_transformed[j].id){
    x_landmark=landmark_Range[k].x;
    y_landmark=landmark_Range[k].y;
  }//end if
  
	 particles[i].weight *= (1/2*M_PI*std_landmark[0]*std_landmark[1]) * exp(-(pow(x_transformed-x_landmark,2)/2*std_landmark[0]*std_landmark[0]+pow(y_transformed-y_landmark,2)/2*std_landmark[1]*std_landmark[1]));
   } //end k
   }//end j
   
  Obs_transformed.clear();
  landmark_Range.clear();
 } //end i

} // end function

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
   //Max weight
  double mw=numeric_limits<double>::min();
  vector<double> weights;
  
  for(unsigned int i=0;i<num_particles;i++){
    weights.push_back(particles[i].weight);
  	if(particles[i].weight > mw){
    	mw=particles[i].weight;
    } //end if
  } //end for
  //Distributions
  std::random_device rand; 
  std::mt19937 generator(rand());
  std::discrete_distribution<int> int_distrib(weights.begin(),weights.end());
  std::uniform_real_distribution<double> double_distrib(0.0,mw);
    
  //using code explained by sebastian in lesson 4.20
  vector<Particle> new_p(num_particles);
  
  int index;
  double b;
  
  for (int i=0;i<num_particles;i++){
    index = int_distrib(generator);
    b=0.0;
    
    b+=2.0*double_distrib(generator);
    while(b>weights[index]){
      b-=weights[index];
      index=(1+index) % num_particles;
    }
    new_p[i]=particles[index];
   //new_p.push_back(particles[index]);
  }
  particles=new_p;
 
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}