#include <FEMCollision.h>
#include <MPRGPSolver.h>
#include "MprgpFemSolver.h"
using namespace MATH;
USE_PRJ_NAMESPACE

MprgpFemSolver::MprgpFemSolver():FemSolverExt(3){

  boost::shared_ptr<FEMCollision> coll(new BVHFEMCollision);
  // boost::shared_ptr<FEMCollision> coll(new SBVHFEMCollision);
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

void MprgpFemSolver::buildVarOffset(){}

void MprgpFemSolver::initVelPos(const double dt){}

void MprgpFemSolver::handleCollDetection(){}

void MprgpFemSolver::forward(const double dt){}

void MprgpFemSolver::forwardSimple(const double dt){}

void MprgpFemSolver::buildLinearSystem(Eigen::SparseMatrix<double> &LHS, VectorXd &RHS, const double dt){}

double MprgpFemSolver::updateVelPos(const VectorXd &new_pos, const double dt){}

void MprgpFemSolver::updateMesh(const double dt){}

void MprgpFemSolver::print()const{

  int num_var = 0;
  for(int i=0;i<_mesh->nrB();i++){
	assert(_mesh->getB(i)._system);
	num_var += _mesh->getB(i)._system->size();
  }

  INFO_LOG("number of nodes: "<<num_var/3);
  // INFO_LOG("number of tets: "<<); /// @todo
  INFO_LOG("newton max outter it: "<<_tree.get<int>("maxIter"));
  INFO_LOG("newton outter tol: "<<_tree.get<double>("eps"));
  INFO_LOG("newton max inner it: "<<newton_inner_max_it);
  INFO_LOG("newton inner tol: "<<newton_inner_tol);
  INFO_LOG("mprgp max it: "<<mprgp_max_it);
  INFO_LOG("mprgp tol: "<<mprgp_tol);
  INFO_LOG("use simple simulation: "<< (use_simple_sim ? "true":"false"));
  collider->print();
}

#define BLK(IN,I) (IN).block(_offVar[i],0,_offVar[i+1]-_offVar[i],1)
#define BLKL(IN,I) (IN).block(_offVarL[i],0,_offVarL[i+1]-_offVarL[i],1)
#define BLKF(IN,I) (IN).block(_offVarF[i],0,_offVarF[i+1]-_offVarF[i],1)
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
	INFO("Newton Iteration")
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
	if(selfColl)
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
	  INFOV("\tResidue: %f",DELTA.cwiseAbs().maxCoeff())
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
