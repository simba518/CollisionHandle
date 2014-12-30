#ifndef _MPRGPFEMSOLVER_H_
#define _MPRGPFEMSOLVER_H_

#include <boost/filesystem.hpp>
#include <FEMSystem.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

class MprgpFemSolver:public FEMSolver{
	
public:
  MprgpFemSolver();
  void useSimpleSimulation(const bool use){
	this->use_simple_sim = use;
  }
  void setVel(const Vector3d &vel, const int body_id);
  void advance(const double dt);
  void setLinearSolverParameters(const double mprgp_tol, const int mprgp_it);
  void setFriction(const double mu_s, const double mu_k){
	collider->setFriction(mu_s, mu_k);
  }
  const vector<size_t> &getVarOffset()const{return off_var;}
  const VVVec4d &getLinearCon()const{return collider->getLinearCon();}
  int currentFrame()const{
	return current_frame;
  }

  void setTargetFold(const string &fold_for_saving_results){
	save_results_to = fold_for_saving_results;
	boost::filesystem::create_directory(save_results_to);
	boost::filesystem::create_directory(save_results_to+"/QP/");
	boost::filesystem::create_directory(save_results_to+"/collisions/");
  }
  const string &saveResultsTo()const{
	return save_results_to;
  }
  void print()const;

protected:
  void handleCollDetection();
  void buildVarOffset();
  void initVelPos(const double dt);
  void updateMesh(const double dt);
  void forward(const double dt);
  void forwardSimple(const double dt);
  void buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt);
  double updateVelPos(const VectorXd &new_pos, const double dt);

private:
  boost::shared_ptr<LinearConCollider> collider;
  size_t num_var;
  vector<size_t> off_var;
  VectorXd pos0, pos1, vel0, vel1;
  TRIPS HTrips, UTrips;
  double mprgp_tol;
  int mprgp_max_it;
  int newton_inner_max_it;
  double newton_inner_tol;
  int current_frame;
  string save_results_to;
  bool use_simple_sim;
};
  
#endif /*_MPRGPFEMSOLVER_H_*/
