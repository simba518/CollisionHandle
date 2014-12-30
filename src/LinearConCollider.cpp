#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
#include "TetrahedronVertexCon.h"
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

  const int vert_id = b->_offset/3 + v->_index;
  assert_in(vert_id, 0, (int)linear_con.size());
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

  assert_eq(force.size(), (int)all_lambdas.size()*3);
  assert_in(vert_id, 0, all_lambdas.size());
  assert_in(plane_id, 0, (int)all_lambdas[vert_id].size()-1);

  const double lambda = all_lambdas[vert_id][plane_id];
  force.segment<3>(vert_id*3) += normal*lambda;
}

void GeomConCache::addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &all_lambdas, 
									  VectorXd &force, double mu_s, double mu_k)const{

  assert_eq(vel.size(), force.size());
  assert_eq(force.size(), (int)all_lambdas.size()*3);
  assert_in(vert_id, 0, (int)all_lambdas.size());
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

bool GeomConCache::addConPlane(VVec4d &planes, const Vector4d &p){

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

bool SelfConCache::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,
						  boost::shared_ptr<FEMVertex> v,const Vec4& bary, VVVec4d &linear_con){
  
  // compute v[i], d[i]
  this->v[0] = bv->_offset + v->_index;
  initial_x.segment<3>(0) = v->_pos;
  for (int i = 0; i < 4; ++i){
	this->v[i+1] = bv->_offset + c->_v[i]->_index;
	d[i] = bv->_offset + c->_v[i]->_matDist;
	initial_x.segment<3>(i*3+3) = c->_v[i]->_pos;
  }
  
  // convert to linear constraints
  const size_t old_size = linear_con.size();
  convertToLinearCon(linear_con);
  return old_size < linear_con.size();
}

double SelfConCache::fun(const Vector15d &x)const{
  return TetrahedronVertexCon::functionValue(d, &x[0]);
}

void SelfConCache::grad(const Vector15d &x, Vector15d &g)const{
  TetrahedronVertexCon::gradient(d, &x[0], &g[0]);
}

void SelfConCache::findFeasiblePoints(Vector15d &x)const{
  
  Matrix<double, 16, 16> A;
  for (int i = 0; i < 15; ++i){
    A(i,i) = 1;
  }
  A(15,15) = 0;

  Vector16d x_old, x_new, b;
  x_old.head<15>() = initial_x;
  b.head<15>() = initial_x;

  Vector15d J;
  const double tol = 1e-3;
  const int max_it = 100;
  for (int i = 0; i < max_it; ++i){
	grad(x_old.head<15>(), J);
	A.block<15,1>(1,0) = J;
	A.block<1,15>(0,1) = J.transpose();
	b[15] = J.dot(x_old.head<15>()) - fun(x_old.head<15>());
	x_new = A.lu().solve(b);
	const double err = (x_new-x_old).norm();
	if(err <= tol)
	  break;
	x_old = x_new;
  }
  x = x_new.head<15>();
}

void SelfConCache::convertToLinearCon(VVVec4d &linear_con){

  Vector15d x, g;
  findFeasiblePoints(x);
  const double f = fun(x);
  grad(x,g);
  Vector4d plane;
  for (int i = 0; i < 5; ++i){
	n[i] = g.segment<3>(i*3);
	const double norm = n[i].norm();
	assert_gt(norm, 0);
	const double pi = (f-n[i].dot(x.segment<3>(i*3)))/norm;
	n[i] = n[i]/norm;
	plane.segment<3>(0) = n[i];
	plane[3] = pi;
	if( GeomConCache::addConPlane(linear_con[v[i]], plane) ){
	  plane_id[i] = linear_con[v[i]].size()-1;
	}
  }
}

double SelfConCache::computeLambda(const double lambda[5])const{

  return (lambda[0]+lambda[1]+lambda[2]+lambda[3]+lambda[4])/5.0f;
}

void SelfConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{

  double lambda[5];
  for (int i = 0; i < 5; ++i){
	lambda[i] = all_lambdas[v[0]][plane_id[i]];
  }
  const double la = computeLambda(lambda);
  assert_le(la, 0.0f);
  for (int i = 0; i < 5; ++i){
	force.segment<3>(v[i]*3) += la*n[i];
  }
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

void LinearConCollider::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,
							   boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary){

  SelfConCache one_con;
  if ( one_con.handle(bc, c, bv, v, bary, linear_con) ){
	self_con.push_back(one_con);
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
