#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <Simulator.h>
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(SimulatorTest)

BOOST_AUTO_TEST_CASE(test_one_tet){
  
  Simulator simulator;
  simulator.init("./test_case/test_data/one_tet/collision_ball.ini");
  simulator.print();
  simulator.run();

  boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();

  double ek, ep;
  fem_solver->getSystemEnergy(ek, ep);

  ASSERT_EQ_TOL(ek, 18.06418368811757, 1e-10);
  ASSERT_EQ_TOL(ep, -28.73280337280945, 1e-10);

  const VVVec4d &linear_con = fem_solver->getLinearCon();

  ASSERT_EQ(linear_con.size(), 4);
  ASSERT_EQ(linear_con[0].size(), 0);
  ASSERT_EQ(linear_con[1].size(), 1);
  ASSERT_EQ(linear_con[2].size(), 1);
  ASSERT_EQ(linear_con[3].size(), 0);

  Vector4d p1,p2;
  p1 << -0.9161066774901623,  -0.289000350170286,  0.2778980983371203,   2.270172353019375;
  p2 << -0.9527175769240369, -0.2890043417724823, 0.09383873963692302,   2.551793027809889;
  ASSERT_EQ_SMALL_VEC_TOL(linear_con[1][0], p1, 4, 1e-10);
  ASSERT_EQ_SMALL_VEC_TOL(linear_con[2][0], p2, 4, 1e-10);

  // VectorXd pos, vel, accel;
  // fem_solver->getPos(pos);
  // fem_solver->getPos(vel);
  // fem_solver->getPos(accel);
  
}

BOOST_AUTO_TEST_CASE(test_coarse_beam_in_ball){
 
  Simulator simulator;
  simulator.init("./test_case/test_data/beam-coarse/collision_in_ball.ini");
  simulator.print();
  simulator.run();

  boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();

  double ek, ep;
  fem_solver->getSystemEnergy(ek, ep);

  ASSERT_EQ_TOL(ek, 3.192594279389013, 1e-10);
  ASSERT_EQ_TOL(ep, 879.7698309439567, 1e-10);

}

BOOST_AUTO_TEST_CASE(test_coarse_beam_out_ball){
 
  Simulator simulator;
  simulator.init("./test_case/test_data/beam-coarse/collision_out_ball.ini");
  simulator.print();
  simulator.run();

  boost::shared_ptr<const MprgpFemSolver> fem_solver = simulator.getFemSolver();
  double ek, ep;
  fem_solver->getSystemEnergy(ek, ep);

  ASSERT_EQ_TOL(ek, 44.56998214487999, 1e-10);
  ASSERT_EQ_TOL(ep, 69.75710437612339, 1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
