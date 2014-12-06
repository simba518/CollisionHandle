#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

bool GeomConCache::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,
						  const Vec3& n, VVVec4d &linear_con, VectorXd &pos0){

  if(n.norm() < ScalarUtil<double>::scalar_eps)
	return false;

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
  const bool added = addConPlane(linear_con[vert_id], plane);

  if(added){

	pos0[vert_id*3+0] = x[0];
	pos0[vert_id*3+1] = x[1];
	pos0[vert_id*3+2] = x[2];

	this->vert_id = vert_id;
	this->plane_id = linear_con[vert_id].size()-1;
	normal = plane.segment<3>(0);
  }

  return added;
}

void GeomConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{

  assert_eq(force.size(), all_lambdas.size()*3);
  assert_in(vert_id, 0, all_lambdas.size());
  assert_in(plane_id, 0, (int)all_lambdas[vert_id].size()-1);

  const double lambda = all_lambdas[vert_id][plane_id];
  force.segment<3>(vert_id*3) += normal*lambda;
}

void GeomConCache::addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &all_lambdas, 
									  VectorXd &force, double mu_s, double mu_k)const{

  assert_eq(vel.size(), force.size());
  assert_eq(force.size(), all_lambdas.size()*3);
  assert_in(vert_id, 0, all_lambdas.size());
  assert_in(plane_id, 0, (int)all_lambdas[vert_id].size()-1);

  const double lambda = all_lambdas[vert_id][plane_id];
  if(lambda <= 0){
	ERROR_LOG_COND(setprecision(10)<<"lambda should be >= 0, lambda = "<<lambda, (lambda>=0) );
	return;
  }

  const Vector3d ft = force.segment<3>(vert_id*3) - normal.dot(force.segment<3>(vert_id*3))*normal;
  const Vector3d vt = vel.segment<3>(vert_id*3) - normal.dot(vel.segment<3>(vert_id*3))*normal;
  assert_ge(mu_k, 0.0f);
  assert_ge(mu_s, 0.0f);

  if(vt.norm() > ScalarUtil<double>::scalar_eps ){
	const Vector3d fr = ( -min( lambda*mu_k, ft.norm() ) )*vt.normalized();
	force.segment<3>(vert_id*3) += fr;
  }else if(ft.norm() > ScalarUtil<double>::scalar_eps){
	const Vector3d fr = ( -min( lambda*mu_s, ft.norm() ) )*ft.normalized();
	force.segment<3>(vert_id*3) += fr;
  }

  /// @todo update velocity

}

bool GeomConCache::addConPlane(VVec4d &planes, const Vector4d &p)const{

  assert_eq(p,p);
  size_t k = 0;
  for (; k < planes.size(); ++k){
	if ((p-planes[k]).norm() <= 1e-4)
	  break;
  }
  if (k == planes.size()){
	planes.push_back(p);
	return true;
  }else{
	return false;
  }
}

bool SelfConCache::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){
  
  /// @todo
  return false;
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

void SelfConCache::addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &all_lambdas, 
									  VectorXd &force, double mu_s, double mu_k)const{
  
  /// @todo
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){

  GeomConCache one_con;
  if( one_con.handle(b, v, n, linear_con, pos0) )
	geom_con.push_back(one_con);
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){

  SelfConCache one_con;
  if ( one_con.handle(b, v, coef, nrV) ){
	self_con.push_back(one_con);
	one_con.convertToLinearCon(linear_con);
  }
}

void LinearConCollider::addJordanForce(VectorXd &force)const{

  for(size_t i = 0; i < geom_con.size(); i++)
	geom_con[i].addJordanForce(all_lambdas, force);

  for(size_t i = 0; i < self_con.size(); i++)
	self_con[i].addJordanForce(all_lambdas, force);
}

void LinearConCollider::addFrictionalForce(const VectorXd &vel, VectorXd &force)const{
  
  for(size_t i = 0; i < geom_con.size(); i++)
	geom_con[i].addFrictionalForce(vel, all_lambdas, force, friction_s, friction_k);

  for(size_t i = 0; i < self_con.size(); i++)
	self_con[i].addFrictionalForce(vel, all_lambdas, force, friction_s, friction_k);
}

void LinearConCollider::print()const{
  
  INFO_LOG("static friction: "<< friction_s);
  INFO_LOG("kinetic friction: "<< friction_k);
}
