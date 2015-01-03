#include <FEMCollision.h>
#include <MPRGPSolver.h>
#include "MprgpFemSolver.h"
using namespace MATH;
USE_PRJ_NAMESPACE

#define BLK(IN,I) (IN).block(_offVar[i],0,_offVar[i+1]-_offVar[i],1)
#define BLKL(IN,I) (IN).block(_offVarL[i],0,_offVarL[i+1]-_offVarL[i],1)
#define BLKF(IN,I) (IN).block(_offVarF[i],0,_offVarF[i+1]-_offVarF[i],1)

MprgpFemSolver::MprgpFemSolver():FemSolverExt(3,2){

  collider = boost::shared_ptr<LinearConCollider>(new LinearConCollider(feasible_pos));
  mprgp_max_it = 1000;
  mprgp_tol = 1e-4;
}

void MprgpFemSolver::advance(const double dt){

  FUNC_TIMER();
  buildVarOffset();
  initPos(dt);
  handleCollDetection();
  initVel(dt);
  forward(dt);
  updateMesh(dt);
  current_frame ++;
}

void MprgpFemSolver::buildVarOffset(){

  _mesh->buildOffset();
  assembleOffVar();
}

void MprgpFemSolver::initPos(const double dt){

  x0.resize(nrVar());
  x1.resize(nrVar());
  X0.resize(nrVarL());
  X1.resize(nrVarL());

  for(int i=0; i<_mesh->nrB(); i++) {

	const FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec xB,XB;
	sys.getPos(xB);
	sys.getPosL(XB);
	BLK(x0,i)=xB;
	BLKL(X0,i)=XB;
  }
}

void MprgpFemSolver::handleCollDetection(){

  feasible_pos = x0;
  collider->reset();

  ostringstream oss3;
  oss3 << saveResultsTo()+"/collisions/coll_"<< currentFrame() << ".vtk";
  DebugFEMCollider coll_debug( oss3.str(), 3 );

  for(int i=0; i<_mesh->nrB(); i++){
	_mesh->getB(i)._system->beforeCollision();
  }
  if(_geom){
	DEBUG_FUN( _mesh->getColl().collideGeom( *_geom,coll_debug,true ) );
	_mesh->getColl().collideGeom(*_geom,*collider,true);
  }
  if(_tree.get<bool>("selfColl")){
	DEBUG_FUN( _mesh->getColl().collideMesh(coll_debug, true) );
	_mesh->getColl().collideMesh(*collider,true);
  }

  const bool find_feasible = MATH::findFeasible(getLinearCon(), feasible_pos, true);
  assert(find_feasible);
}

void MprgpFemSolver::initVel(const double dt){

  PHI.resize(nrVarL());
  PSI.resize(nrVarL());

  const double gamma=_tree.get<double>("gamma");
  const double beta2=_tree.get<double>("beta2");
  for(int i=0; i<_mesh->nrB(); i++) {

	const FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec xB,XB,VB,AB;
	sys.getPosL(XB);
	sys.getVelL(VB);
	sys.getAccelL(AB);

	BLKL(PHI,i)=VB+(1.0f-gamma)*dt*AB;
	BLKL(PSI,i)=XB+dt*VB+(1.0f-2.0f*beta2)*dt*dt*0.5f*AB;
  }
}

void MprgpFemSolver::forward(const double dt){

  const double eps=_tree.get<double>("eps");
  const int maxIter=_tree.get<int>("maxIter");

  Vec RHS(nrVar());
  _LHS.reset(nrVar(),nrVar(),false);
  _U.reset(nrVarF(),nrVar(),false);

  VVVec4d empty_con(nrVar()/3);
  PlaneProjector<double> projector(getLinearCon(), feasible_pos);
  PlaneProjector<double> projector_no_con(empty_con, feasible_pos);
  SparseMatrix<double> LHS_mat;

  for(int i = 0; i < maxIter; i++) {

	buildLinearSystem(LHS_mat, RHS, dt);
	solve(LHS_mat, RHS, projector, projector_no_con);
	if (updatePos() < eps) 
	  break;
  }
}

void MprgpFemSolver::buildLinearSystem(SparseMatrix<double>&LHS_mat,VectorXd&RHS,const double dt){

  _LHS.clear();
  _U.clear();

  for ( int i = 0; i < _mesh->nrB(); i++) {
	const FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec xB,XB;
	sys.getPos(xB);
	sys.getPosL(XB);
	BLK(x1,i)=xB;
	BLKL(X1,i)=XB;
  }

  const double beta=_tree.get<double>("beta");
  for ( int i = 0; i < _mesh->nrB(); i++) {

	const double bdt = beta*dt;
	const FEMSystem& sys=*(_mesh->getB(i)._system);
	const Vec MRHSLB=BLKL(PSI-X1,i);
	const Vec CRHSB=BLK(x0-x1,i)*bdt;
	Vec RHSB=Vec::Zero(sys.size());
	sys.buildSystem(1.0f,bdt*dt,bdt,_LHS,_U,MRHSLB,CRHSB,bdt*dt,RHSB,_offVar[i]);
	BLK(RHS,i)=RHSB;
  }

  const boost::shared_ptr<const TRIPS> trips = _LHS.getSparse();
  LHS_mat.resize(_LHS.rows(), _LHS.cols());
  LHS_mat.setFromTriplets(trips->begin(),trips->end());

  assert_eq(LHS_mat.rows(), LHS_mat.cols());
  assert_eq(LHS_mat.cols(), x1.size());
  RHS += LHS_mat*x1;
}

void MprgpFemSolver::solve(const SparseMatrix<double> &LHS_mat, VectorXd &RHS, 
						   PlaneProjector<double> &projector, 
						   PlaneProjector<double> &projector_no_con){

  const FixedSparseMatrix<double> A(LHS_mat);
  
  // use constraints, no frictional and collision forces.
  {
	new_pos = feasible_pos;
	typedef DiagonalPlanePreconSolver<double,FixedSparseMatrix<double>,true> NoPreconditioner;
	const int rlst_code=MPRGPPlane<double>::solve<FixedSparseMatrix<double>, NoPreconditioner>
	  ( A, RHS, projector, new_pos, mprgp_tol, mprgp_max_it);
	ERROR_LOG_COND("MPRGP is not convergent, result code is "<<rlst_code<<endl,rlst_code==0);
	DEBUG_FUN( MPRGPPlane<double>::checkResult(LHS_mat, RHS, projector, new_pos, mprgp_tol));
  }

  // compute frictional and collision forces
  {
	const VectorXd diff = LHS_mat*new_pos-RHS;
	collider->computeAllLambdas( diff, projector.getFaceIndex() );
	collider->addJordanForce(RHS);
	// collider->addFrictionalForce(vel1, RHS); ///@todo no friction.
  }

  // use frictional and collision forces, and no constraints.
  {
	const int rlst_code = MPRGPPlane<double>::solve( A, RHS, projector_no_con, new_pos, mprgp_tol, mprgp_max_it);
	ERROR_LOG_COND("MPRGP is not convergent, result code is "<<rlst_code<<endl,rlst_code==0);
  }
}

double MprgpFemSolver::updatePos(){

  const VectorXd DELTA = new_pos-x1;
  x1 = new_pos;
  for(int i=0; i<_mesh->nrB(); i++) {
	FEMSystem& sys=*(_mesh->getB(i)._system);
	sys.setPos(BLK(x1,i));
  }

  double tol_err = 0;
  if(DELTA.size() > 0){
	tol_err = DELTA.cwiseAbs().maxCoeff();
  }

  return tol_err;
}

void MprgpFemSolver::updateMesh(const double dt){

  const double beta=_tree.get<double>("beta");
  const double gamma=_tree.get<double>("gamma");
  for(int i=0; i<_mesh->nrB(); i++) {
	FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec XB,VB,AB;
	sys.getPosL(XB);
	AB=(XB-BLKL(PSI,i))/(beta*dt*dt);
	VB=BLKL(PHI,i)+AB*gamma*dt;

	sys.setVelL(VB);
	sys.setAccelL(AB);
  }
  _mesh->updateMesh();
}

void MprgpFemSolver::print()const{

  // INFO_LOG("number of nodes: "<<nrVar()/3);
  // INFO_LOG("number of tets: "<<); /// @todo
  INFO_LOG("newton it: "<<_tree.get<int>("maxIter"));
  INFO_LOG("newton tol: "<<_tree.get<double>("eps"));
  INFO_LOG("mprgp it: "<<mprgp_max_it);
  INFO_LOG("mprgp tol: "<<mprgp_tol);
  INFO_LOG("use simple simulation: "<< (use_simple_sim ? "true":"false"));
  collider->print();
}

void FemSolverExtDebug::advance(const double dt){

  //initialize
  _mesh->buildOffset();
  const double collK=_tree.get<double>("collK");
  const double beta2=_tree.get<double>("beta2");
  const double beta=_tree.get<double>("beta");
  const double gamma=_tree.get<double>("gamma");
  const double eps=_tree.get<double>("eps");
  const int maxIter=_tree.get<int>("maxIter");

  //handle collision detection
  ostringstream oss3;
  oss3 << saveResultsTo()+"/collisions/coll_"<< currentFrame() << ".vtk";
  DebugFEMCollider coll_debug( oss3.str(), 3 );

  const bool selfColl=_tree.get<bool>("selfColl");
  for(int i=0; i<_mesh->nrB(); i++)
	_mesh->getB(i)._system->beforeCollision();
  DefaultFEMCollider coll(*_mesh,collK,1.0f,1.0f);
  if(_geom){
	DEBUG_FUN( _mesh->getColl().collideGeom( *_geom,coll_debug,true ) );
	_mesh->getColl().collideGeom(*_geom,coll,true);
  }
  if(selfColl){
	DEBUG_FUN( _mesh->getColl().collideMesh(coll_debug, true) );
	_mesh->getColl().collideMesh(coll,true);
  }

  Vec fE=coll.getFE();
  Eigen::SparseMatrix<double,0,long int> HE;
  coll.getHE(HE);
  for(int i=0; i<_mesh->nrB(); i++)
	_mesh->getB(i)._system->afterCollision();

  //find variable offset
  assembleOffVar();

  //initialize velocity and position
  Vec x0(nrVar()),x1(nrVar());
  Vec X0(nrVarL()),X1(nrVarL()),PHI(nrVarL()),PSI(nrVarL());
  for(int i=0; i<_mesh->nrB(); i++) {
	const FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec xB,XB,VB,AB;
	sys.getPos(xB);
	sys.getPosL(XB);
	sys.getVelL(VB);
	sys.getAccelL(AB);

	BLK(x0,i)=xB;
	BLKL(X0,i)=XB;
	BLKL(PHI,i)=VB+(1.0f-gamma)*dt*AB;
	BLKL(PSI,i)=XB+dt*VB+(1.0f-2.0f*beta2)*dt*dt*0.5f*AB;
  }

  if(_tree.get<bool>("debug",false))
	INFO("Newton Iteration");
	  //main loop: we use Implicit Newmark Scheme
  Vec RHS(nrVar()),DELTA;
  _LHS.reset(nrVar(),nrVar(),_tree.get<bool>("denseSystem"));
  _U.reset(nrVarF(),nrVar(),_tree.get<bool>("denseSystem"));
  for(int i=0; i<maxIter; i++) {
	_LHS.clear();
	_U.clear();
	for(int i=0; i<_mesh->nrB(); i++) {
	  const FEMSystem& sys=*(_mesh->getB(i)._system);
	  Vec xB,XB;
	  sys.getPos(xB);
	  sys.getPosL(XB);
	  BLK(x1,i)=xB;
	  BLKL(X1,i)=XB;
	}

	//build LHS and RHS
	for(int i=0; i<_mesh->nrB(); i++) {
	  //tell system to build: M+(beta*dt*dt)*K+(gamma*dt)*C with off
	  //tell system to build: U with offU
	  //tell system to build: M(dSn+((1-gamma)*dt)*ddSn-dS_{n+1})+(gamma*dt)*f_I-C*dS_{n+1} with offU[1]
	  const FEMSystem& sys=*(_mesh->getB(i)._system);
	  Vec RHSB=Vec::Zero(sys.size());
	  Vec MRHSLB=BLKL(PSI-X1,i);
	  Vec CRHSB=BLK(x0-x1,i)*beta*dt;
	  sys.buildSystem(1.0f,beta*dt*dt,beta*dt,_LHS,_U,
					  MRHSLB,CRHSB,beta*dt*dt,RHSB,_offVar[i]);
	  BLK(RHS,i)=RHSB;
	}

	//add external term
	if(selfColl){
	  if(_U.isDense()) {
		Matd UDense;
		_U.getDense(UDense);
		RHS+=UDense.transpose()*(fE-HE*LtoF(X1-X0))*(beta*dt*dt);
		_LHS.addDense(0,0,UDense.transpose()*(HE*UDense).eval()*(beta*dt*dt));
	  } else {
		Eigen::SparseMatrix<double,0,long int> USparse;
		_U.getSparse(USparse);
		RHS+=USparse.transpose()*(fE-HE*LtoF(X1-X0))*(beta*dt*dt);
		_LHS.addSparse(0,0,USparse.transpose()*(HE*USparse).eval()*(beta*dt*dt));
	  }
	}

	//solve
	DELTA=_LHS.solve(RHS);
	x1+=DELTA;
	for(int i=0; i<_mesh->nrB(); i++) {
	  FEMSystem& sys=*(_mesh->getB(i)._system);
	  sys.setPos(BLK(x1,i));
	}

	//exit test
	if(DELTA.size() == 0)
	  break;
	if(_tree.get<bool>("debug",false))
	  INFOV("\tResidue: %f",DELTA.cwiseAbs().maxCoeff());
	if(DELTA.cwiseAbs().maxCoeff() < eps)
	  break;
  }

  //update velocity and acceleration
  for(int i=0; i<_mesh->nrB(); i++) {
	FEMSystem& sys=*(_mesh->getB(i)._system);
	Vec XB,VB,AB;
	sys.getPosL(XB);
	AB=(XB-BLKL(PSI,i))/(beta*dt*dt);
	VB=BLKL(PHI,i)+AB*gamma*dt;

	sys.setVelL(VB);
	sys.setAccelL(AB);
  }
  _mesh->updateMesh();

  current_frame ++;
}
#undef BLK
#undef BLKL
#undef BLKF
