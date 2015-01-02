#include <FEMCollision.h>
#include <MPRGPSolver.h>
#include "MprgpFemSolver.h"
using namespace MATH;
USE_PRJ_NAMESPACE

MprgpFemSolver::MprgpFemSolver():FEMSolver(3){

  // boost::shared_ptr<FEMCollision> coll(new BVHFEMCollision);
  boost::shared_ptr<FEMCollision> coll(new SBVHFEMCollision);
  collider = boost::shared_ptr<LinearConCollider>(new LinearConCollider(feasible_pos));

  _mesh.reset(new FEMMesh(3,coll));
  resetImplicitEuler();
  setCollK(1E5f);
  _mesh->setCellSz(1.0f);
  setSelfColl(false);
  mprgp_max_it = 1000;
  mprgp_tol = 1e-4;
  newton_inner_max_it = 1;
  newton_inner_tol = 1e-4;

  current_frame = 0;
  setTargetFold( "./tempt");
  use_simple_sim = false;
}

void MprgpFemSolver::advance(const double dt){

  FUNC_TIMER();
  buildVarOffset();
  initVelPos(dt);
  handleCollDetection();
  if(use_simple_sim){
	forwardSimple(dt);
  }else{
	forward(dt);
  }
  updateMesh(dt);
  current_frame ++;
}

void MprgpFemSolver::buildVarOffset(){

  assert(_mesh);
  _mesh->buildOffset();
  num_var = 0;
  off_var.resize(_mesh->nrB());
  for(int i=0;i<_mesh->nrB();i++){
	assert(_mesh->getB(i)._system);
	FEMSystem& sys=*(_mesh->getB(i)._system);
	off_var[i] = num_var;
	num_var += sys.size();
  }
}

void MprgpFemSolver::initVelPos(const double dt){

  pos0.resize(num_var);
  vel0.resize(num_var);
  for ( int i=0; i<_mesh->nrB(); i++ ){

	const FEMSystem& sys = *(_mesh->getB(i)._system);
	Vec posB,velB,accelB,XB;
	sys.getPos(posB);
	sys.getVel(velB);
	sys.getAccel(accelB);
	_mesh->getB(i).getPos(XB);

	posB=posB+dt*velB+(1.0f-2.0f*_beta2)*dt*dt*0.5f*accelB;
	velB=velB+(1.0f-_gamma)*dt*accelB;

	// sys.setPos(posB);
	// sys.setVel(velB);

	pos0.block(off_var[i],0,sys.size(),1)=posB;
	vel0.block(off_var[i],0,sys.size(),1)=velB;
  }
}

void MprgpFemSolver::handleCollDetection(){

  feasible_pos = pos0;
  for ( int i = 0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->beforeCollision();
  collider->reset();

  ostringstream oss3;
  oss3 << saveResultsTo()+"/collisions/coll_"<< currentFrame() << ".vtk";
  DebugFEMCollider coll_debug( oss3.str(), 3 );

  if ( _geom ){
	DEBUG_FUN( _mesh->getColl().collideGeom( *_geom,coll_debug,true ) );
	_mesh->getColl().collideGeom( *_geom,*collider,true );
  }

  if ( _selfColl ) {
	DEBUG_FUN( _mesh->getColl().collideMesh(coll_debug, true) );
	_mesh->getColl().collideMesh(*collider, true);
  }

  for ( int i = 0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->afterCollision();

  const bool find_feasible = MATH::findFeasible(getLinearCon(), feasible_pos, true);
  assert(find_feasible);

  for ( int i=0; i<_mesh->nrB(); i++ ){
  	FEMSystem& sys = *(_mesh->getB(i)._system);
  	sys.setPos( pos0.block(off_var[i],0,sys.size(),1) );
  	sys.setVel( vel0.block(off_var[i],0,sys.size(),1) );
  }
}

void MprgpFemSolver::forward(const double dt){

  VectorXd RHS;
  Eigen::SparseMatrix<double> LHS;

  pos1 = pos0;
  vel1 = vel0;

  VVVec4d empty_con(num_var/3);
  PlaneProjector<double> projector(getLinearCon(), feasible_pos);
  PlaneProjector<double> projector_no_con(empty_con, feasible_pos);

  VectorXd new_pos;  
  for(int i = 0; i < _maxIter ; i++){

	// build linear system
	buildLinearSystem(LHS, RHS, dt);
	const FixedSparseMatrix<double> A(LHS);

	// use constraints, no frictional and collision forces.
	{
	  new_pos = feasible_pos;
	  typedef DiagonalPlanePreconSolver<double,FixedSparseMatrix<double>, true > NoPreconditioner;
	  const int rlst_code = MPRGPPlane<double>::solve<FixedSparseMatrix<double>, NoPreconditioner>
		( A, RHS, projector, new_pos, mprgp_tol, mprgp_max_it);
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<<rlst_code<<endl, rlst_code == 0);
	  DEBUG_FUN( MPRGPPlane<double>::checkResult(LHS, RHS, projector, new_pos, mprgp_tol) );
	}

	// compute frictional and collision forces
	{
	  const VectorXd diff = LHS*new_pos-RHS;
	  collider->computeAllLambdas( diff, projector.getFaceIndex() );
	  collider->addJordanForce(RHS);
	  collider->addFrictionalForce(vel1, RHS); ///@todo no friction.
	}

	// use frictional and collision forces, and no constraints.
	{
	  const int rlst_code = MPRGPPlane<double>::solve( A, RHS, projector_no_con, new_pos, mprgp_tol, mprgp_max_it);
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<<rlst_code<<endl, rlst_code == 0);
	}

	// check tolerance
	if ( updateVelPos (new_pos, dt) < _eps ){
	  break;
	}

  }
}

void MprgpFemSolver::forwardSimple(const double dt){

  // a simple forward method for test.
  VectorXd RHS;
  Eigen::SparseMatrix<double> LHS;

  pos1 = pos0;
  vel1 = vel0;

  VectorXd new_pos;
  PlaneProjector<double> projector(getLinearCon(), feasible_pos);
  for ( int i = 0; i < _maxIter ; i++ ) {

	buildLinearSystem(LHS, RHS, dt);
	const FixedSparseMatrix<double> A(LHS);
	new_pos = feasible_pos;

	DEBUG_FUN({ostringstream ossm_qp;
		ossm_qp << saveResultsTo() << "/QP/frame_" << currentFrame() << "_it_"<< i << ".b";
		writeQP<double>(LHS, RHS, projector.getPlanes(), new_pos, ossm_qp.str());} );

	const int code=MPRGPPlane<double>::solve(A,RHS,projector,new_pos,mprgp_tol,mprgp_max_it);

	ERROR_LOG_COND("MPRGP is not convergent, result code is "<<code<<endl,code==0);
	DEBUG_FUN( MPRGPPlane<double>::checkResult(LHS, RHS, projector, new_pos, mprgp_tol));

	if ( updateVelPos (new_pos, dt) < _eps )
	  break;
  }
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
	Vec CRHSB=-vel1.block(off_var[i],0,sys.size(),1)*(_gamma*dt);

	sys.buildSystem(1.0,_beta*dt*dt,_gamma*dt,HTrips,UTrips,
					MRHSB,KRHSB,CRHSB,_gamma*dt,RHSB,off_var[i]);

	RHS.block(off_var[i],0,sys.size(),1)=RHSB;
  }

  LHS.setFromTriplets(HTrips.begin(),HTrips.end());

  // using pos as variables instead of delta_velocity
  RHS *= ((_beta*dt)/_gamma);
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

void MprgpFemSolver::setLinearSolverParameters(const double mprgp_tol,const int mprgp_it){
  
  assert_gt(mprgp_it, 0);
  assert_gt(mprgp_tol, 0.0);
  this->mprgp_max_it = mprgp_it;
  this->mprgp_tol = mprgp_tol;
}

void MprgpFemSolver::print()const{

  int num_var = 0;
  for(int i=0;i<_mesh->nrB();i++){
	assert(_mesh->getB(i)._system);
	num_var += _mesh->getB(i)._system->size();
  }

  INFO_LOG("number of nodes: "<<num_var/3);
  // INFO_LOG("number of tets: "<<); /// @todo
  INFO_LOG("newton max outter it: "<<_maxIter);
  INFO_LOG("newton outter tol: "<<_eps);
  INFO_LOG("newton max inner it: "<<newton_inner_max_it);
  INFO_LOG("newton inner tol: "<<newton_inner_tol);
  INFO_LOG("mprgp max it: "<<mprgp_max_it);
  INFO_LOG("mprgp tol: "<<mprgp_tol);
  INFO_LOG("use simple simulation: "<< (use_simple_sim ? "true":"false"));
  collider->print();
}

void MprgpFemSolver::setVel(const Vector3d &vel, const int body_id){

  assert_in(body_id, 0, _mesh->nrB()-1);
  assert( _mesh->getB(body_id)._system );

  FEMSystem& sys = *(_mesh->getB(body_id)._system);
  VectorXd velB;
  sys.getVel(velB);
  for(int j = 0; j < velB.size(); j += 3){
	velB.segment<3>(j) = vel;
  }
  sys.setVel(velB);
}
