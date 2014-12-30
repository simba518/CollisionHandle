#ifndef _LINEARCONCOLLIDER_H_
#define _LINEARCONCOLLIDER_H_

#include <vector>
#include <FEMCollider.h>
#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <MPRGPSolver.h>
using namespace Eigen;
using namespace std;
USE_PRJ_NAMESPACE

typedef vector<Vector4d,Eigen::aligned_allocator<Vector4d> > VVec4d;
typedef vector<VVec4d > VVVec4d;
typedef Matrix<double,15,1> Vector15d;
typedef Matrix<double,16,1> Vector16d;

class GeomConCache{
  
public:
  GeomConCache(){
	vert_id = plane_id = -1;
	normal.setZero();
  }
  bool handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n,
			  VVVec4d &linear_con, VectorXd &pos0);
  void addJordanForce(const vector<vector<double> > &lambda, VectorXd &force)const;
  void addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &lambda, VectorXd &force, double mu_s, double mu_k)const;
  
  static bool addConPlane(VVec4d &con_planes, const Vector4d &p);
  
private:
  int vert_id;
  int plane_id;
  Vector3d normal;
};

// self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl),
// where a,b,c is the weight coordinates.
class SelfConCache{

public:
  bool handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,
			  boost::shared_ptr<FEMVertex> v,const Vec4& bary, VVVec4d &linear_con);
  void convertToLinearCon(VVVec4d &linear_con);
  void addJordanForce(const vector<vector<double> > &lambda, VectorXd &force)const;
  void addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &lambda, VectorXd &force, double mu_s, double mu_k)const;

protected:
  double computeLambda(const double lambda[5])const;
  void findFeasiblePoints(Vector15d &x)const;
  double fun(const Vector15d &x)const;
  void grad(const Vector15d &x, Vector15d &g)const;

private:  
  Vector15d initial_x;
  int v[5];
  double d[4];
  int plane_id[5];
  Vector3d n[5];
};

/**
 * Handle both self and geometric collisions, and return them as linear constraints.
 * 1. self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl).
 * 2. geometric collision constraints: c = n^t*xi+p.
 */
class LinearConCollider:public FEMCollider{
  
public:
  LinearConCollider(VectorXd &pos0):pos0(pos0){
	setFriction(0.5f, 0.4f);
	reset();
  }

  void reset(){

	const size_t num_verts = pos0.size()/3;
	linear_con.clear();
	linear_con.resize(num_verts);
	self_con.clear();
	geom_con.clear();
	all_lambdas.clear();
  }

  void setFriction(double mu_s, double mu_k){
	assert_ge(mu_k, 0.0f);
	assert_ge(mu_s, mu_k);
	friction_s = mu_s;
	friction_k = mu_k;
  }

  void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV) {}

  void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);

  void handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary);

  void computeAllLambdas(const VectorXd &K_x_f, const vector<vector<int> > &face){
	MATH::MPRGPPlane<double>::computeLagMultipliers(K_x_f, linear_con, face, all_lambdas);
  }

  void addJordanForce(VectorXd &force)const;

  void addFrictionalForce(const VectorXd &vel, VectorXd &force)const;

  const VVVec4d &getLinearCon()const{return linear_con;}

  void print()const;

private:
  VVVec4d linear_con;
  vector<SelfConCache> self_con;
  vector<GeomConCache> geom_con;
  vector<vector<double> > all_lambdas;
  VectorXd &pos0;
  double friction_s;
  double friction_k;
};

#endif /* _LINEARCONCOLLIDER_H_ */
