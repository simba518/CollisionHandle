#ifndef _LINEARCONCOLLIDER_H_
#define _LINEARCONCOLLIDER_H_

#include <FEMCollider.h>

typedef vector<Vector4d,Eigen::aligned_allocator<Vector4d> > VVec4d;
typedef vector<VVec4d > VVVec4d;

// self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl), 
// where a,b,c is the weight coordinates.
struct SelfConCache{

  int i,j,k,l;
  Vector3d n;
  double a,b,c;
};

/**
 * Handle both self and geometric collisions, and return them as linear constraints.
 * 1. self collision constraints: c = n^t*(xi-a*xj-b*xk-c*xl).
 * 2. geometric collision constraints: c = n^t*xi+p.
 */
class LinearConCollider:public FEMCollider{
  
public:
  LinearConCollider(VVVec4d &geom_con, vector<SelfConCache> &self_con,const size_t num_verts):
	geom_con(geom_con),self_con(self_con){
	self_con.clear();
	geom_con.clear();
	geom_con.resize(num_verts);
  }
  
  void handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){}

  void handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){}
  
private:
  VVVec4d &geom_con;
  vector<SelfConCache> &self_con;
};

class SelfCollHandler{
  
public:
  SelfCollHandler(const vector<SelfConCache> &self_con):self_con(self_con){}
  void addSelfConAsLinearCon(VVVec4d &linear_con){}
  void addJordanForce(const VectorXd &x, const vector<vector<int> > &face, VectorXd &force){}

private:
  const vector<SelfConCache> &self_con;
  
};

#endif /* _LINEARCONCOLLIDER_H_ */
