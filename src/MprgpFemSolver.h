#ifndef _MPRGPFEMSOLVER_H_
#define _MPRGPFEMSOLVER_H_

#include <FEMSystem.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

class MprgpFemSolver:public FEMSolver{
	
public:
  MprgpFemSolver();
  void advance(const double dt);

protected:
  void handleCollDetection();
  void buildVarOffset();
  void initVelPos(const double dt);
  void updateMesh(const double dt);
  void forward(const double dt);
  void buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt);
  double updateVelPos(const VectorXd &new_pos, const double dt);

private:
  size_t num_var;
  vector<size_t> off_var;
  VectorXd pos0, pos1, vel0, vel1;
  VVVec4d linear_con;
  vector<SelfConCache> self_con;
  TRIPS HTrips, UTrips;
};
  
#endif /*_MPRGPFEMSOLVER_H_*/
