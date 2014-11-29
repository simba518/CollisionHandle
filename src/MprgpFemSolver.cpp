#include "MprgpFemSolver.h"

void MprgpFemSolver::advance(const double dt){
	
  buildVarOffset();
  initVelPos();
  handleCollDetection();
  forward(dt);
  updateMesh(dt);
}

void MprgpFemSolver::handleCollDetection(){
  
  for ( size_t i=0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->beforeCollision();

  LinearConCollider coll( linear_con, self_con, _mesh->nrV() );

  if ( _geom )
	_mesh->getColl().collideGeom( *_geom,coll,true );

  if ( _selfColl ) {
	_mesh->getColl().collideMesh(coll,true);
	SelfCollHandler self_coll(self_con);
	self_coll.addSelfConAsLinearCon(linear_con);
  }

  for ( size_t i=0; i < _mesh->nrB(); i++ )
	_mesh->getB(i)._system->afterCollision();
}

void MprgpFemSolver::buildVarOffset(){
  
  num_var = 0;
  off_var.resize(_mesh->nrB());
  for(size_t i=0;i<_mesh->nrB();i++){
	FEMSystem& sys=*(_mesh->getB(i)._system);
	off_var[i] = num_var;
	num_var += sys.size();
  }
}

void MprgpFemSolver::initVelPos(){
  
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

void MprgpFemSolver::updateMesh(const double dt){
  
  vel1=(vel1-vel0)/(_gamma*dt);
  for(size_t i=0;i<_mesh->nrB();i++){
	FEMSystem& sys= *( _mesh->getB(i)._system );
	sys.setAccel( vel1.block(off_var[i],0,sys.size(),1) );
  }
  _mesh->updateMesh();
}

void MprgpFemSolver::forward(const double dt){

  VectorXd RHS;
  Eigen::SparseMatrix<double> LHS;

  pos1 = pos0;
  vel1 = vel0;
  for ( size_t i = 0; i < _maxIter ; i++ ) {

	buildLinearSystem(LHS, RHS);
	VectorXd new_pos;
	MPRGPPlane::solve( LHS, RHS, linear_con, new_pos );
	if ( updateVelPos (new_pos) < _eps )
	  break;
  }
}

void MprgpFemSolver::buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS){
  
  RHS.resize(num_var);
  LHS.resize(num_var,num_var);

  //build LHS and RHS
  HTrips.clear();
  UTrips.clear();

  for ( size_t i = 0; i < _mesh->nrB(); i++ ) {

	//tell system to build: M+(beta*dt*dt)*K+(gamma*dt)*C with off
	//tell system to build: M(dSn+((1-gamma)*dt)*ddSn-dS_{n+1})+(gamma*dt)*f_I-C*dS_{n+1}
	FEMSystem& sys=*(_mesh->getB(i)._system);

	const Vec RHSB=Vec::Zero(sys.size());
	const Vec MRHSB=(vel0-vel1).block(off_var[i],0,sys.size(),1);
	const Vec KRHSB=Vec::Zero(sys.size());
	const Vec CRHSB=-vel1.block(off_var[i],0,sys.size(),1)*_gamma*dt;

	sys.buildSystem(1.0f,_beta*dt*dt,_gamma*dt,HTrips,UTrips,
					MRHSB,KRHSB,CRHSB,_gamma*dt,RHSB,off_var[i]);

	RHS.block(off_var[i],0,sys.size(),1)=RHSB;
  }

  LHS.setFromTriplets(HTrips.begin(),HTrips.end());

  // using pos as variables instead of delta_velocity
  LHS *= (_gamma/(_beta*dt));
  RHS += (LHS*pos1);

}

double MprgpFemSolver::updateVelPos(const VectorXd &new_pos){
  
  const VectorXd DELTA = (new_pos-pos1)*(_gamma/(_beta*dt));
  vel1 += DELTA;
  pos1 = new_pos;
  for(size_t i=0;i<_mesh->nrB();i++){

	FEMSystem& sys=*(_mesh->getB(i)._system);
	sys.setPos(pos1.block(off_var[i],0,sys.size(),1));
	sys.setVel(vel1.block(off_var[i],0,sys.size(),1));
  }
  return DELTA.cwiseAbs().maxCoeff();
}
