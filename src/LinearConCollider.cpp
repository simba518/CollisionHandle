#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){

  SelfConCache one_con;
  one_con.handle(b, v, coef, nrV);
  self_con.push_back(one_con);
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){

  if (n.norm() < ScalarUtil<double>::scalar_eps){
	return;
  }

  Vector3d N, V;
  N[0] = n[0];
  N[1] = n[1];
  N[2] = n[2];

  V[0] = v->_pos[0];
  V[1] = v->_pos[1];
  V[2] = v->_pos[2];
  
  Vector4d plane;
  plane.segment<3>(0) = N;
  plane.segment<3>(0).normalize();

  const Vector3d x = V+N*(1.0+ScalarUtil<double>::scalar_eps);
  plane[3] = -(plane.segment<3>(0).dot(x));

  const int vert_id = b->_offset + v->_index;
  assert_in(vert_id, 0, linear_con.size());
  assert_eq_ext(plane, plane, "n:" << n.transpose());
  addConPlane(linear_con[vert_id], plane);
  
  pos0[vert_id*3+0] = x[0];
  pos0[vert_id*3+1] = x[1];
  pos0[vert_id*3+2] = x[2];
}

void LinearConCollider::addConPlane(VVec4d &planes, const Vector4d &p)const{

  assert_eq(p,p);
  size_t k = 0;
  for (; k < planes.size(); ++k){
	if ((p-planes[k]).norm() <= 1e-4)
	  break;
  }
  if (k == planes.size())
	planes.push_back(p);
}

void GeomConCache::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){
  
  
}

void GeomConCache::addJordanForce(const vector<vector<double> > &lambda, VectorXd &force)const{
  
  /// @todo
}

void SelfConCache::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){
  
  /// @todo
  
}

double SelfConCache::computeLambda(const double lambda[4])const{
  
  /// @todo 
  return 0;  
}

void SelfConCache::convertToLinearCon(VVVec4d &linear_con){
  
  /// @todo
}

void SelfConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{

  double lambda[4];
  lambda[0] = all_lambdas[i][pi];
  lambda[1] = all_lambdas[i][pj];
  lambda[2] = all_lambdas[i][pk];
  lambda[3] = all_lambdas[i][pl];
  const Vector3d N = computeLambda(lambda)*n;
  force.segment<3>(i*3) += N;
  force.segment<3>(j*3) -= N*a;
  force.segment<3>(k*3) -= N*b;
  force.segment<3>(l*3) -= N*c;
}

void SelfCollHandler::addSelfConAsLinearCon(VVVec4d &linear_con){
  
  for(size_t i = 0; i < self_con.size(); i++)
	self_con[i].convertToLinearCon(linear_con);
}

void SelfCollHandler::addJordanForce(const VectorXd &K_x_f, const vector<vector<int> > &face, 
									 const VVVec4d &linear_con, VectorXd &force){

  const int num_verts = face.size();
  assert_eq(K_x_f.size(), num_verts*3);
  assert_eq(force.size(), num_verts*3);

  vector<vector<double> > all_lambdas(num_verts);
  for(int i = 0; i < num_verts; i++){
	const Vector3d Kxf_i = K_x_f.segment<3>(i*3);
	computeLambdas(Kxf_i, face[i], linear_con[i], all_lambdas[i]);
  }

  for(size_t i = 0; i < self_con.size(); i++)
	self_con[i].addJordanForce(all_lambdas, force);
}

void SelfCollHandler::computeLambdas(const Vector3d &Kxf_i, const vector<int> &face_i, 
									 const VVec4d &planes, vector<double> &lambdas)const{
  
  /// @todo
  
}
