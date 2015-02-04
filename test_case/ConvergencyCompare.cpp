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
  // const string qp = project_dir+"/data/cube/tempt_2cube_mprgp14000/QP/qp6.b";
  const string qp = project_dir+"/data/dino/tempt_selfcon_mprgp/QP/qp45.b";
  
  SparseMatrix<double> A, J;
  VectorXd B, c, init_x;
  UTILITY::Timer timer;

  const bool load_succ = loadQP(A, B, J, c, init_x, qp);
  assert(load_succ);

  // print information
  if(true){
	cout << "A.rows() = " << A.rows() << endl;
	cout << "J.rows() = " << J.rows() << endl;
	cout << "A.nonZeros() = " << A.nonZeros() << endl;
	cout << "J.nonZeros() = " << J.nonZeros() << endl;
  }

  // analysis matrix A
  if(false){
	VectorXd eig_vec;
	double eig_val = 0.0f;
	EIGEN3EXT::largestEigen(A,eig_vec,eig_val);
	INFO_LOG("largest eigen value: " << eig_val);
	EIGEN3EXT::smallestEigen(A,eig_vec,eig_val,(int)(A.rows()/10),1e-3*A.norm());
	INFO_LOG("smallest eigen value: " << eig_val);
  }

  const double tol = 1e-5*B.norm();
  const int max_it = 10000;

  // mprgp with decoupled general con
  if(true){

	const SparseMatrix<double> JJt_mat = J*J.transpose();
	assert_eq_ext(JJt_mat.nonZeros(), J.rows(), "Matrix J is not decoupled.\n" << J);
	VectorXd JJt;
	MATH::getDiagonal(JJt_mat, JJt);
	DecoupledConProjector<double> projector(J, JJt, c);
	  
	VectorXd x;
	projector.project(init_x, x);

	typedef FixedSparseMatrix<double> MAT;
	MAT FA(A);
	const int code = MPRGPDecoupledCon<double>::solve<MAT,true>(FA, B, J, c, x, tol, max_it, "decoupled mprgp");
	ERROR_LOG_COND("MPRGP is not convergent, result code is "<< code << endl, code==0);
  }

  // ica with con
  if(true){

	ConjugateGradient<SparseMatrix<double> > cg_solver;
	timer.start();
	const VectorXd uncon_x = cg_solver.compute(A).solve(B);
	timer.stop("eigen-cg solving time: ");

	VectorXd p, x(uncon_x.size());
	p = c - J*uncon_x;

	timer.start();
	ICASolver ica_solver(max_it, tol);
	ica_solver.setName("ica with con");
	ica_solver.reset(A);
	x.setZero();
	// x = init_x - uncon_x;
	const bool succ = ica_solver.solve(J, p, x);
	timer.stop("ica with con solving time: ");
	ica_solver.printSolveInfo(A, J, p, x);
	ERROR_LOG_COND("ICA is not convergent, (iterations, residual) = " << 
				   ica_solver.getIterations() << ", " << ica_solver.getResidual(), succ);
  }

  // mprgp wihout con
  if(false){
	  
  }

  // gs without con
  if(false){
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

BOOST_AUTO_TEST_SUITE_END()
