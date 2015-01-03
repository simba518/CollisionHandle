#ifndef _MPRGPFEMSOLVER_H_
#define _MPRGPFEMSOLVER_H_

#include <boost/filesystem.hpp>
#include <FEMSystem.h>
#include <MPRGPSolver.h>
#include "LinearConCollider.h"
using namespace MATH;
USE_PRJ_NAMESPACE

class FemSolverExt:public FEMSolver{
	
public:
  FemSolverExt(sizeType dim=3,sizeType cOption=2): FEMSolver(dim, cOption){
	current_frame = 0;
	use_simple_sim = false;
	setTargetFold( "./tempt");
  }
  void useSimpleSimulation(const bool use){
	this->use_simple_sim = use;
  }
  void setVel(const Vector3d &vel, const int body_id){
	// assert_in(body_id, 0, _mesh->nrB()-1);
	// assert( _mesh->getB(body_id)._system );

	// FEMSystem& sys = *(_mesh->getB(body_id)._system);
	// VectorXd velB;
	// sys.getVel(velB);
	// for(int j = 0; j < velB.size(); j += 3){
	// 	velB.segment<3>(j) = vel;
	// }
	// sys.setVel(velB);
  }
  void setTargetFold(const string &fold_for_saving_results){
	save_results_to = fold_for_saving_results;
	boost::filesystem::create_directory(save_results_to);
	boost::filesystem::create_directory(save_results_to+"/QP/");
	boost::filesystem::create_directory(save_results_to+"/collisions/");
  }

  virtual void advance(const double dt){
	FEMSolver::advance(dt);
	current_frame ++;
  }
  virtual void setLinearSolverParameters(const double mprgp_tol, const int mprgp_it){}
  virtual void setFriction(const double mu_s, const double mu_k){}

  const string &saveResultsTo()const{
	return save_results_to;
  }
  int currentFrame()const{
	return current_frame;
  }
  virtual void print()const{
	int num_var = 0;
	for(int i = 0; i < _mesh->nrB(); i++){
	  assert(_mesh->getB(i)._system);
	  num_var += _mesh->getB(i)._system->size();
	}
	INFO_LOG("number of nodes: "<<num_var/3);
	INFO_LOG("FEM Solver: FemSolverExt");
	// INFO_LOG("newton max outter it: "<<_tree.get<int>("maxIter"));
	// INFO_LOG("newton outter tol: "<<_tree.get<double>("eps"));
  }

protected:
  int current_frame;
  string save_results_to;
  bool use_simple_sim;
};

class MprgpFemSolver:public FemSolverExt{
	
public:
  MprgpFemSolver();
  void advance(const double dt);
  void setLinearSolverParameters(const double mprgp_tol, const int mprgp_it){
	assert_gt(mprgp_it, 0);
	assert_gt(mprgp_tol, 0.0);
	this->mprgp_max_it = mprgp_it;
	this->mprgp_tol = mprgp_tol;
  }
  void setFriction(const double mu_s, const double mu_k){
	collider->setFriction(mu_s, mu_k);
  }

  const VVVec4d &getLinearCon()const{return collider->getLinearCon();}
  void print()const;

protected:
  void buildVarOffset();
  void initPos(const double dt);
  void handleCollDetection();
  void initVel(const double dt);
  double updatePos();
  void updateMesh(const double dt);
  void forward(const double dt);
  void solve(const SparseMatrix<double> &LHS_mat, VectorXd &RHS, 
			 PlaneProjector<double> &projector, PlaneProjector<double> &projector_no_con);
  void buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt);

private:
  boost::shared_ptr<LinearConCollider> collider;
  VectorXd x0, x1, X0, X1, PHI, PSI, feasible_pos, new_pos;
  double mprgp_tol;
  int mprgp_max_it;
};

class FemSolverExtDebug:public FemSolverExt{
	
public:
  FemSolverExtDebug():FemSolverExt(){}
  void advance(const double dt);
};
  
#endif /*_MPRGPFEMSOLVER_H_*/
