#ifndef _MPRGPFEMSOLVER_H_
#define _MPRGPFEMSOLVER_H_

#include <FEMSystem.h>

class MprgpFemSolver:public FEMSolver{
	
public:
  MprgpFemSolver(){
	
	const int dim = 3;
	boost::shared_ptr<FEMCollision> coll(new SBVHFEMCollision);
	_mesh.reset(new FEMMesh(dim,coll));
	resetImplicitEuler();
	setCollK(1E5f);
	_mesh->setCellSz(1.0f);
	setSelfColl(true);
  }
  void advance(const double dt);

protected:
  void handleCollDetection();
  void buildVarOffset();
  void initVelPos();
  void updateMesh(const double dt);
  void forward(const double dt);
  void buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS);
  double updateVelPos(const VectorXd &new_pos);

private:
  size_t num_var;
  vector<size_t> off_var;
  VectorXd pos0, pos1, vel0, vel1;
  VVVec4d linear_con;
  vector<SelfConCache> self_con;
  TRIPS HTrips, UTrips;
};
  
#endif /*_MPRGPFEMSOLVER_H_*/
