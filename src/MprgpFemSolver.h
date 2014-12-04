#ifndef _MPRGPFEMSOLVER_H_
#define _MPRGPFEMSOLVER_H_

#include <FEMSystem.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

class MprgpFemSolver:public FEMSolver{
	
public:
  MprgpFemSolver();
  void advance(const double dt);
  void setLinearSolverParameters(double mprgp_tol, int mprgp_it);
  const vector<size_t> &getVarOffset()const{return off_var;}
  const VVVec4d &getLinearCon()const{return collider->getLinearCon();}

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
};
  
#endif /*_MPRGPFEMSOLVER_H_*/
