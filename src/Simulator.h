#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include <string>
#include "MprgpFemSolver.h"

class Simulator{
	
public:
  Simulator() {fem_solver = boost::shared_ptr<MprgpFemSolver>(new MprgpFemSolver());}
  void init(const string &json_file);
  void run();

  double timeStep()const{
	return time_step;
  }
  size_t totalFrames()const{
	return total_frames;
  }
  const string &saveResultsTo()const{
	return save_results_to;
  }
  
private:
  string save_results_to;
  double time_step;
  size_t total_frames;
  boost::shared_ptr<MprgpFemSolver> fem_solver;
};
  
#endif /*_SIMULATOR_H_*/
