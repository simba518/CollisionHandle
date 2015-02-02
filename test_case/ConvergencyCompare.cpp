#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MPRGPSolver.h>
#include <ICASolver.h>
#include <FlexiblePCG.h>
#include <EigenSolver.h>
using namespace Eigen;
using namespace MATH;

BOOST_AUTO_TEST_SUITE(ConvergencyCompare)

BOOST_AUTO_TEST_CASE(test_iterative_solvers){
  
  const string project_dir = "/home/simba/Workspace/CollisionHandle/";
  // const string qp = project_dir+"/data/cube/tempt_selfcon_mprgp60/QP/qp9.b";
  // const string qp = project_dir+"/data/cube/tempt_selfcon_mprgp300/QP/qp9.b";
  // const string qp = project_dir+"/data/cube/tempt_selfcon_mprgp2000/QP/qp9.b";
  const string qp = project_dir+"/data/cube/tempt_selfcon_mprgp14000/QP/qp9.b";
  // const string qp = project_dir+"/data/cube/tempt_selfcon_mprgp14000/QP/qp7.b";
  // const string qp = project_dir+"/data/dino/tempt_selfcon_mprgp/QP/qp29.b";
  // const string qp = project_dir+"/data/dragon/tempt_plane_mprgp/QP/qp26.b";
  // const string qp = project_dir+"/data/dragon/tempt_plane_mprgp/QP/qp60.b";
  // const string qp = project_dir+"/data/dragon/tempt_plane_mprgp/QP/qp110.b";
  
  SparseMatrix<double> A;
  VectorXd B, init_x;
  VVVec4d planes;
  UTILITY::Timer timer;

  if ( loadQP(A, B, planes, init_x, qp) ){

	// analysis matrix A
	{
	  VectorXd eig_vec;
	  double eig_val = 0.0f;
	  EIGEN3EXT::largestEigen(A,eig_vec,eig_val);
	  INFO_LOG("largest eigen value: " << eig_val);
	  EIGEN3EXT::smallestEigen(A,eig_vec,eig_val,(int)(A.rows()/10),1e-3*A.norm());
	  INFO_LOG("smallest eigen value: " << eig_val);
	}

	const double tol = 2e-4*B.norm();
	const int max_it = 10000;

	// mprgp with con
	if(true){
	  VectorXd x = init_x;
	  PlaneProjector<double> projector(planes, x);
	  const FixedSparseMatrix<double> FA(A);
	  const int code = MPRGPPlane<double>::solve(FA, B, projector, x, tol, max_it, "mprgp with con");
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
	  MPRGPPlane<double>::checkResult(A, B, projector, x, tol);
	}

	// ica with con
	if(true){
	  ConjugateGradient<SparseMatrix<double> > cg_solver;
	  timer.start();
	  const VectorXd uncon_x = cg_solver.compute(A).solve(B);
	  timer.stop("eigen-cg solving time: ");

	  SparseMatrix<double> J;
	  VectorXd c, p, x(uncon_x.size());
	  convert(planes, J, c);
	  p = c - J*uncon_x;

	  timer.start();
	  ICASolver ica_solver(max_it, tol);
	  ica_solver.setName("ica with con");
	  ica_solver.reset(A);
	  // x.setZero();
	  x = init_x - uncon_x;
	  const bool succ = ica_solver.solve(J, p, x);
	  timer.stop("ica with con solving time: ");
	  ica_solver.printSolveInfo(A, J, p, x);
	  ERROR_LOG_COND("ICA is not convergent, (iterations, residual) = " << 
					 ica_solver.getIterations() << ", " << ica_solver.getResidual(), succ);
	}

	// pcg wihout con: solve Ax = B
	if(true){
	  typedef DiagonalPlanePreconSolver<double,FixedSparseMatrix<double>,false>Preconditioner;
	  typedef FlexiblePCG<double, SparseMatrix<double>, VectorXd, Preconditioner, true> PCG;

	  const FixedSparseMatrix<double> SA(A);
	  VVVec4d planes_for_each_node(init_x.size()/3);
	  VectorXd x = init_x;
	  PlaneProjector<double> projector(planes_for_each_node, x);
	  Preconditioner precond(SA, projector.getFace(), projector.getPlanes(), projector.getFaceIndex());
	  PCG solver(A, precond, tol, max_it);
	  solver.setName("pcg without con");

	  timer.start();
	  const int code = solver.solve(B, x);
	  timer.stop("my cg solving time: ");

	  assert_eq(code ,0);
	}
	
	// mprgp wihout con
	if(false){
	  VectorXd x = init_x;
	  VVVec4d planes_for_each_node(init_x.size()/3);
	  PlaneProjector<double> projector(planes_for_each_node, x);
	  const FixedSparseMatrix<double> FA(A);
	  const int code = MPRGPPlane<double>::solve(FA, B, projector, x, tol, max_it, "mprgp without con");
	  ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
	  MPRGPPlane<double>::checkResult(A, B, projector, x, tol);
	}

	// gs without con
	if(true){
	  GaussSeidel<false> GS(max_it, tol);
	  timer.start();
	  GS.setName("gs without con");
	  GS.reset(A);
	  VectorXd x = init_x;
	  const bool succ = GS.solve(B, x);
	  timer.stop("GS solving time: ");
	  GS.printSolveInfo(A, B, x);
	  ERROR_LOG_COND("GS is not convergent\n", succ);
	}

	// block-gs without con
	if(false){
	  BlockGaussSeidel BlockGS(max_it, tol);
	  BlockGS.setName("block without con");
	  BlockGS.reset(A);
	  VectorXd x = init_x;
	  const bool succ = BlockGS.solve(B, x);
	  // BlockGS.printSolveInfo(A, B, x);
	  ERROR_LOG_COND("Block GS is not convergent\n", succ);
	}

  }
}

BOOST_AUTO_TEST_SUITE_END()
