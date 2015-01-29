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
  FemSolverExt(sizeType cOption=2): FEMSolver(3, cOption){

	_tree.put<bool>("denseSystem",false);
	current_frame = 0;
	setTargetFold( "./tempt");
	solver_name = "FemSolverExt";
	debug_coll = false;
	use_iterative_solver = true;
  }
  void setVel(const Vector3d &vel, const int body_id);
  void setTargetFold(const string &fold_for_saving_results){
	save_results_to = fold_for_saving_results;
	boost::filesystem::create_directory(save_results_to);
	boost::filesystem::create_directory(save_results_to+"/QP/");
	boost::filesystem::create_directory(save_results_to+"/collisions/");
  }

  virtual void advance(const double dt);
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
	INFO_LOG("FEM Solver: "<< solver_name);
	INFO_LOG("number of nodes: "<<num_var/3);
	// INFO_LOG("number of tets: "<<); /// @todo
	INFO_LOG("newton it: "<<_tree.get<int>("maxIter"));
	INFO_LOG("newton tol: "<<_tree.get<double>("eps"));
	INFO_LOG("debug collision: "<< (debug_coll ? "true":"false"));
	INFO_LOG("use iterative solver: "<< (use_iterative_solver ? "true":"false"));
  }

protected:
  void solve(const FEMSystemMatrix &LHS, const Vec &RHS, Vec &DELTA);

protected:
  string solver_name;
  int current_frame;
  string save_results_to;
  bool debug_coll;
  bool use_iterative_solver;
  SparseMatrix<double> LHS_mat;
};

class MprgpFemSolver:public FemSolverExt{
	
public:
  MprgpFemSolver(int cOption=2);
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
  virtual void forward(const double dt);

  void buildVarOffset();
  void initPos(const double dt);
  void handleCollDetection();
  void initVel(const double dt);
  double updatePos();
  void updateMesh(const double dt);
  void solve(const SparseMatrix<double> &LHS_mat, VectorXd &RHS, 
			 PlaneProjector<double> &projector, PlaneProjector<double> &projector_no_con);
  void buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt);

protected:
  boost::shared_ptr<LinearConCollider> collider;
  VectorXd x0, x1, X0, X1, PHI, PSI, feasible_pos, new_pos;
  double mprgp_tol;
  int mprgp_max_it;
};

class MoseckFemSolver:public MprgpFemSolver{

public:
  MoseckFemSolver(int cOption=2):MprgpFemSolver(cOption){
	solver_name = "MoseckFemSolver";
  }

protected:
  void forward(const double dt);

};

class ICAFemSolver:public MprgpFemSolver{

public:
  ICAFemSolver(int cOption=2):MprgpFemSolver(cOption){
	solver_name = "ICAFemSolver";
  }

protected:
  void forward(const double dt);

};

#endif /*_MPRGPFEMSOLVER_H_*/
