#include <float.h>
#include <limits>
#include <MPRGPSolver.h>
using namespace MATH;

float ScalarUtil<float>::scalar_max=FLT_MAX;
float ScalarUtil<float>::scalar_eps=1E-5f;

double ScalarUtil<double>::scalar_max=DBL_MAX;
double ScalarUtil<double>::scalar_eps=1E-9;

#include "Simulator.h"

int main(int argc, char *argv[]){

  assert_ge(argc,2);
  Simulator simulator;
  simulator.init(argv[1]);
  simulator.run();
  return 0;
}
