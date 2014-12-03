#include <FEMCollision.h>
#include <MPRGPSolver.h>
#include "MprgpFemSolver.h"
using namespace MATH;
USE_PRJ_NAMESPACE

MprgpFemSolver::MprgpFemSolver():FEMSolver(3){
	
  boost::shared_ptr<FEMCollision> coll(new SBVHFEMCollision);
  _mesh.reset(new FEMMesh(3,coll));
  resetImplicitEuler();
  setCollK(1E5f);
  _mesh->setCellSz(1.0f);
  setSelfColl(false);
}

void MprgpFemSolver::advance(const double dt){
	
  buildVarOffset();
  initVelPos(dt);
  handleCollDetection();
  forward(dt);
  updateMesh(dt);
}

void MprgpFemSolver::buildVarOffset(){

  assert(_mesh);
  _mesh->buildOffset();
  num_var = 0;
  off_var.resize(_mesh->nrB());
  for(size_t i=0;i<_mesh->nrB();i++){
	assert(_mesh->getB(i)._system);
	FEMSystem& sys=*(_mesh->getB(i)._system);
	off_var[i] = num_var;
	num_var += sys.size();
  }
}

void MprgpFemSolver::initVelPos(const double dt){

  pos0.resize(num_var);
  vel0.resize(num_var);
  for ( size_t i=0; i<_mesh->nrB(); i++ ){

	FEMSystem& sys = *(_mesh->getB(i)._system);
	Vec posB,velB,accelB,XB;
	sys.getPos(posB);
	sys.getVel(velB);
	sys.getAccel(accelB);
	_mesh->getB(i).getPos(XB);

	posB=posB+dt*velB+(1.0f-2.0f*_beta2)*dt*dt*0.5f*accelB;
	velB=velB+(1.0f-_gamma)*dt*accelB;
	sys.setPos(posB);
	sys.setVel(velB);

	pos0.block(off_var[i],0,sys.size(),1)=posB;
	vel0.block(off_var[i],0,sys.size(),1)=velB;
  }
}

void MprgpFemSolver::handleCollDetection(){
  
  for ( size_t i=0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->beforeCollision();

  linear_con.resize(_mesh->nrV());
  
  static int frame = 0;
  ostringstream oss3;
  oss3 << "./tempt/coll_"<< frame++ << ".vtk";
  DebugFEMCollider coll_debug(oss3.str(),3);

  LinearConCollider coll( linear_con, self_con, geom_con, pos0);

  if ( _geom ){
	_mesh->getColl().collideGeom( *_geom,coll,true );
	_mesh->getColl().collideGeom( *_geom,coll_debug,true );
  }

  if ( _selfColl ) {

	_mesh->getColl().collideMesh(coll, true);
	_mesh->getColl().collideMesh(coll_debug, true);

	SelfCollHandler self_coll(self_con);
	self_coll.addSelfConAsLinearCon(linear_con);
  }

  for ( size_t i=0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->afterCollision();
}

void MprgpFemSolver::forward(const double dt){

  VectorXd RHS;
  Eigen::SparseMatrix<double> LHS;

  pos1 = pos0;
  vel1 = vel0;

  VectorXd new_pos;
  PlaneProjector<double> projector(linear_con, pos0);
  for ( size_t i = 0; i < _maxIter ; i++ ) {

	buildLinearSystem(LHS, RHS, dt);
	const FixedSparseMatrix<double> A(LHS);
	new_pos = pos0;
	MPRGPPlane<double>::solve( A, RHS, projector, new_pos );
	if ( updateVelPos (new_pos, dt) < _eps )
	  break;
  }

  if (_selfColl){
	
	SelfCollHandler self_coll(self_con);
	const VectorXd diff = LHS*new_pos-RHS;
	self_coll.addJordanForce(diff, projector.getFaceIndex(), linear_con, RHS);

	_sol.compute(LHS);
	ERROR_LOG_COND("Factorization Fail!", _sol.info() == Eigen::Success);
	new_pos = _sol.solve(RHS);
  }

  updateVelPos (new_pos, dt);
}

void MprgpFemSolver::buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt){
  
  RHS.resize(num_var);
  LHS.resize(num_var,num_var);

  //build LHS and RHS
  HTrips.clear();
  UTrips.clear();

  for ( int i = 0; i < _mesh->nrB(); i++ ) {

	//tell system to build: M+(beta*dt*dt)*K+(gamma*dt)*C with off
	//tell system to build: M(dSn+((1-gamma)*dt)*ddSn-dS_{n+1})+(gamma*dt)*f_I-C*dS_{n+1}
	FEMSystem& sys=*(_mesh->getB(i)._system);

	Vec RHSB=Vec::Zero(sys.size());
	Vec MRHSB=(vel0-vel1).block(off_var[i],0,sys.size(),1);
	Vec KRHSB=Vec::Zero(sys.size());
	Vec CRHSB=-vel1.block(off_var[i],0,sys.size(),1)*_gamma*dt;

	sys.buildSystem(1.0,_beta*dt*dt,_gamma*dt,HTrips,UTrips,
					MRHSB,KRHSB,CRHSB,_gamma*dt,RHSB,off_var[i]);

	RHS.block(off_var[i],0,sys.size(),1)=RHSB;
  }

  LHS.setFromTriplets(HTrips.begin(),HTrips.end());

  // using pos as variables instead of delta_velocity
  LHS *= (_gamma/(_beta*dt));
  RHS += (LHS*pos1);

}

double MprgpFemSolver::updateVelPos(const VectorXd &new_pos, const double dt){
  
  const VectorXd DELTA = (new_pos-pos1)*(_gamma/(_beta*dt));
  vel1 += DELTA;
  pos1 = new_pos;
  for(int i=0;i<_mesh->nrB();i++){

	FEMSystem& sys=*(_mesh->getB(i)._system);
	sys.setPos(pos1.block(off_var[i],0,sys.size(),1));
	sys.setVel(vel1.block(off_var[i],0,sys.size(),1));
  }
  return DELTA.cwiseAbs().maxCoeff();
}

void MprgpFemSolver::updateMesh(const double dt){
  
  vel1=(vel1-vel0)/(_gamma*dt);
  for(int i=0;i<_mesh->nrB();i++){
	FEMSystem& sys= *( _mesh->getB(i)._system );
	sys.setAccel( vel1.block(off_var[i],0,sys.size(),1) );
  }
  _mesh->updateMesh();
}
