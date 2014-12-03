#ifndef _LINEARCONCOLLIDER_H_
#define _LINEARCONCOLLIDER_H_

#include <vector>
#include <FEMCollider.h>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace std;
USE_PRJ_NAMESPACE

typedef vector<Vector4d,Eigen::aligned_allocator<Vector4d> > VVec4d;
typedef vector<VVec4d > VVVec4d;

// self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl), 
// where a,b,c is the weight coordinates.
class SelfConCache{

public:
  SelfConCache(){
	i = j = k = l = -1;
	a = b = c = -1.0f;
	n.setZero();
	x0.setZero();
  }
  void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV);
  void convertToLinearCon(VVVec4d &linear_con);
  void addJordanForce(const vector<vector<double> > &lambda, VectorXd &force)const;
  double computeLambda(const double lambda[4])const;

private:  
  int i,j,k,l;
  int pi, pj, pk, pl;
  double a,b,c;
  Vector3d n, x0;
};

class GeomConCache{
  
public:
  GeomConCache(){
	vert_id = plane_id = -1;
	normal.setZero();
  }
  void addJordanForce(const vector<vector<double> > &lambda, VectorXd &force)const;
  void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);
  
private:
  int vert_id;
  int plane_id;
  Vector3d normal;
};

/**
 * Handle both self and geometric collisions, and return them as linear constraints.
 * 1. self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl).
 * 2. geometric collision constraints: c = n^t*xi+p.
 */
class LinearConCollider:public FEMCollider{
  
public:
  LinearConCollider(VVVec4d &linear_con, vector<SelfConCache> &self_con, 
					vector<GeomConCache> &geom_con, VectorXd &pos0):
	linear_con(linear_con),self_con(self_con),geom_con(geom_con), pos0(pos0){

	const size_t num_verts = pos0.size()/3;
	self_con.clear();
	linear_con.clear();
	linear_con.resize(num_verts);
  }
  
  void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV);

  void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n);

protected: 
  void addConPlane(VVec4d &con_planes, const Vector4d &p)const;
  
private:
  VVVec4d &linear_con;
  vector<SelfConCache> &self_con;
  vector<GeomConCache> &geom_con;
  VectorXd &pos0;
};

class SelfCollHandler{
  
public:
  SelfCollHandler(const vector<SelfConCache> &self_con):self_con(self_con){}
  void addSelfConAsLinearCon(VVVec4d &linear_con);
  void addJordanForce(const VectorXd &x, const vector<vector<int> > &face, const VVVec4d &linear_con, VectorXd &force);

protected:
  void computeLambdas(const Vector3d &K_x_f_i, const vector<int> &face_i, const VVec4d &planes, vector<double> &lambdas)const;

private:
  const vector<SelfConCache> &self_con;
};

#endif /* _LINEARCONCOLLIDER_H_ */
