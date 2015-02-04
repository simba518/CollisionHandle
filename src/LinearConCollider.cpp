#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
#include "TetrahedronVertexCon.h"
#include <ActiveSetQP3D.h>
USE_PRJ_NAMESPACE

bool GeomConCache::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,
						  const Vec3& n, VVVec4d &linear_con, VectorXd &feasible_pos){

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
  const int added = addConPlane(linear_con[vert_id], plane);

  if(added < 0){
	// assert_ext(MATH::isFeasible(linear_con[vert_id], x), "x: "<<x.transpose());
	feasible_pos.segment<3>(vert_id*3) = x;
	this->vert_id = vert_id;
	this->plane_id = (int)linear_con[vert_id].size()-1;
	normal = plane.segment<3>(0);
  }

  return added;
}

void GeomConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{

  assert_eq(force.size(), (int)all_lambdas.size()*3);
  assert_in(vert_id, 0, (int)all_lambdas.size());
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

int GeomConCache::addConPlane(VVec4d &planes, const Vector4d &p){

  assert_eq(p,p);
  int k = 0;
  for (; k < (int)planes.size(); ++k){
	if ((p-planes[k]).norm() <= 1e-4)
	  break;
  }
  if (k == (int)planes.size()){
	planes.push_back(p);
	k = -1;
  }
  return k;
}

void GeomConCache::getConstraints(TRIPS &trips, vector<double> &rhs, const VVVec4d &linear_con)const{
  
  assert_in(vert_id, 0, (int)linear_con.size()-1);
  assert_in(plane_id, 0, (int)linear_con[vert_id].size()-1);
  const double p = linear_con[vert_id][plane_id][3];
  const Vector3d& n = normal;

  const int row = rhs.size();
  const int col0 = vert_id*3;
  trips.push_back( Triplet<double,int>(row, col0+0, n[0]) );
  trips.push_back( Triplet<double,int>(row, col0+1, n[1]) );
  trips.push_back( Triplet<double,int>(row, col0+2, n[2]) );

  rhs.push_back(-p);
}

bool VolumeSelfConCache::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,boost::shared_ptr<FEMBody> bv,
						  boost::shared_ptr<FEMVertex> v,const Vec4& bary, VVVec4d &linear_con){

  // compute v[i], d[i]
  this->v[0] = bv->_offset/3 + v->_index;
  initial_x.segment<3>(0) = v->_pos;
  for (int i = 0; i < 4; ++i){
	this->v[i+1] = bv->_offset/3 + c->_v[i]->_index;
	d[i] = bv->_offset/3 + c->_v[i]->_matDist;
	initial_x.segment<3>(i*3+3) = c->_v[i]->_pos;
  }
  
  // convert to linear constraints
  return convertToLinearCon(linear_con);
}

double VolumeSelfConCache::fun(const Vector15d &x)const{
  return TetrahedronVertexCon::functionValue(d, &x[0]);
}

void VolumeSelfConCache::grad(const Vector15d &x, Vector15d &g)const{
  TetrahedronVertexCon::gradient(d, &x[0], &g[0]);
}

void VolumeSelfConCache::findFeasiblePoints(Vector15d &x)const{
  
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
  assert_eq(x_new, x_new);
  x = x_new.head<15>();
}

bool VolumeSelfConCache::convertToLinearCon(VVVec4d &linear_con){

  bool added = false;
  Vector15d x, g;
  findFeasiblePoints(x);
  const double f = fun(x);
  assert_eq(f, f);
  grad(x,g);
  assert_eq(g, g);
  Vector4d plane;
  for (int i = 0; i < 5; ++i){
	n[i] = g.segment<3>(i*3);
	const double norm = n[i].norm();
	assert_gt(norm, 0);
	const double pi = (f-n[i].dot(x.segment<3>(i*3)))/norm;
	n[i] = n[i]/norm;
	plane.segment<3>(0) = n[i];
	plane[3] = pi;
	if( GeomConCache::addConPlane(linear_con[v[i]], plane) < 0){
	  added = true;
	  plane_id[i] = (int)linear_con[v[i]].size()-1;
	}
  }
  return added;
}

double VolumeSelfConCache::computeLambda(const double lambda[5])const{

  return (lambda[0]+lambda[1]+lambda[2]+lambda[3]+lambda[4])/5.0f;
}

void VolumeSelfConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{

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

void VolumeSelfConCache::addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &all_lambdas, 
									  VectorXd &force, double mu_s, double mu_k)const{
  
  /// @todo
}

void VolumeSelfConCache::getConstraints(TRIPS &trips, vector<double> &rhs, const VVVec4d &linear_con)const{
  
  /// @todo
}

bool SurfaceSelfConCache::handle(boost::shared_ptr<FEMBody> body[5], boost::shared_ptr<FEMVertex> v[5],
								 const Vec3d coef[5], sizeType nrV, VVVec4d &linear_con, VectorXd &feasible_pos){

  i = body[0]->_offset/3 + v[0]->_index;
  j = body[1]->_offset/3 + v[1]->_index;
  k = body[1]->_offset/3 + v[2]->_index;
  l = body[1]->_offset/3 + v[3]->_index;

  const Vector3d &vj = v[1]->_pos;
  const Vector3d &vk = v[2]->_pos;
  const Vector3d &vl = v[3]->_pos;

  // n = coef[0]; /// @bug
  n = ((vk-vj).cross(vl-vj)).normalized();
  assert_eq(n,n);

  a = -coef[1].dot(n);
  b = -coef[2].dot(n);
  c = -coef[3].dot(n);

  x0 = vj*a + vk*b + vl*c;
  x0 -= (n.dot(x0-vj))*n; /// @bug

  convertToLinearCon(linear_con);
  feasible_pos.segment<3>(i*3) = x0;
  assert_eq(linear_con[i].size(), 1);

  // assert_ext(MATH::isFeasible(linear_con[i],x0), "con size: "<<linear_con[i].size());
  // assert_ext(MATH::isFeasible(linear_con[j],feasible_pos.segment<3>(j*3)), "con size: "<<linear_con[j].size()<< ", j="<<j);
  // assert_ext(MATH::isFeasible(linear_con[k],feasible_pos.segment<3>(k*3)), "con size: "<<linear_con[k].size()<< ", k="<<k);
  // assert_ext(MATH::isFeasible(linear_con[l],feasible_pos.segment<3>(l*3)), "con size: "<<linear_con[l].size()<< ", l="<<l);

  return true;
}

void SurfaceSelfConCache::addJordanForce(const vector<vector<double> > &all_lambdas, VectorXd &force)const{
  
  double lambda[4];
  lambda[0] = all_lambdas[i][pi];
  lambda[1] = all_lambdas[j][pj];
  lambda[2] = all_lambdas[k][pk];
  lambda[3] = all_lambdas[l][pl];

  const Vector3d N = computeLambda(lambda)*n;
  force.segment<3>(i*3) += N;
  force.segment<3>(j*3) -= N*a;
  force.segment<3>(k*3) -= N*b;
  force.segment<3>(l*3) -= N*c;
}

void SurfaceSelfConCache::addFrictionalForce(const VectorXd &vel, const vector<vector<double> > &lambda, 
											 VectorXd &force, double mu_s, double mu_k)const{
  /// @todo
}

bool SurfaceSelfConCache::convertToLinearCon(VVVec4d &linear_con){

  // add n*(xi-x0)
  Vector4d plane;
  plane.head(3) = n;
  plane[3] = -n.dot(x0);

  pi = GeomConCache::addConPlane(linear_con[i], plane);
  if(pi < 0) pi = (int)linear_con[i].size()-1;

  // add n*(xj,k,l - x0)
  plane.head(3) = -n;
  plane[3] = -plane[3];

  pj = GeomConCache::addConPlane(linear_con[j], plane);
  if(pj < 0) pj = (int)linear_con[j].size()-1;

  pk = GeomConCache::addConPlane(linear_con[k], plane);
  if(pk < 0) pk = (int)linear_con[k].size()-1;

  pl = GeomConCache::addConPlane(linear_con[l], plane);
  if(pl < 0) pl = (int)linear_con[l].size()-1;

  return true;
}

double SurfaceSelfConCache::computeLambda(const double lambda[4])const{

  const double result = (lambda[0]+lambda[1]+lambda[2]+lambda[3])/(1.0f+a+b+c);
  assert_eq(result, result);
  return result;
}

void SurfaceSelfConCache::getConstraints(TRIPS &trips, vector<double> &rhs, const VVVec4d &linear_con)const{
  
  const int row = rhs.size();
  for (int dim = 0; dim < 3; ++dim){
	trips.push_back( Triplet<double,int>(row, i*3+dim, n[dim]) );
	trips.push_back( Triplet<double,int>(row, j*3+dim, -a*n[dim]) );
	trips.push_back( Triplet<double,int>(row, k*3+dim, -b*n[dim]) );
	trips.push_back( Triplet<double,int>(row, l*3+dim, -c*n[dim]) );
  }
  rhs.push_back(0);
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){

  const int vert_id = b->_offset/3 + v->_index;
  if(!collided(vert_id)){
	GeomConCache one_con;
	if( one_con.handle(b, v, n, linear_con, feasible_pos) ){
	  geom_con.push_back(one_con);
	  coll_as_vert[vert_id] = true;
	}
  }
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> bc,boost::shared_ptr<FEMCell> c,
							   boost::shared_ptr<FEMBody> bv,boost::shared_ptr<FEMVertex> v,const Vec4& bary){

  bool coll = collided(bv->_offset/3 + v->_index);
  for (int i = 0; (!coll) && i < 4; ++i){
	coll = collided(bv->_offset/3 + c->_v[i]->_index);
  }
  if (!coll){
	VolumeSelfConCache one_con;
	if ( one_con.handle(bc, c, bv, v, bary, linear_con) ){
	  vol_self_con.push_back(one_con);
	}
	coll_as_vert[ collided(bv->_offset/3 + v->_index) ] = true;
	for (int i = 0; i < 4; ++i){
	  coll_as_vol[bv->_offset/3 + c->_v[i]->_index] = true;
	}
  }
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],
							   const Vec3d coef[5],sizeType nrV) {

  const int i = b[0]->_offset/3 + v[0]->_index;
  const int j = b[1]->_offset/3 + v[1]->_index;
  const int k = b[1]->_offset/3 + v[2]->_index;
  const int l = b[1]->_offset/3 + v[3]->_index;

  if( (!collided(i))&& (!collided(j))&& (!collided(k))&& (!collided(l)) ){
	SurfaceSelfConCache one_con;
	if ( one_con.handle(b, v, coef, nrV, linear_con, feasible_pos) ){
	  surface_self_con.push_back(one_con);
	  coll_as_vert[i] = true;
	  coll_as_face[j] = true;
	  coll_as_face[k] = true;
	  coll_as_face[l] = true;
	}
  }
}

void LinearConCollider::addJordanForce(VectorXd &force)const{

  for(size_t i = 0; i < geom_con.size(); i++)
	geom_con[i].addJordanForce(all_lambdas, force);

  for(size_t i = 0; i < vol_self_con.size(); i++)
	vol_self_con[i].addJordanForce(all_lambdas, force);

  for(size_t i = 0; i < surface_self_con.size(); i++)
	surface_self_con[i].addJordanForce(all_lambdas, force);
}

void LinearConCollider::addFrictionalForce(const VectorXd &vel, VectorXd &force)const{
  
  for(size_t i = 0; i < geom_con.size(); i++)
	geom_con[i].addFrictionalForce(vel, all_lambdas, force, friction_s, friction_k);

  for(size_t i = 0; i < vol_self_con.size(); i++)
	vol_self_con[i].addFrictionalForce(vel, all_lambdas, force, friction_s, friction_k);

  for(size_t i = 0; i < surface_self_con.size(); i++)
	surface_self_con[i].addFrictionalForce(vel, all_lambdas, force, friction_s, friction_k);
}

void LinearConCollider::getConstraints(SparseMatrix<double> &A, VectorXd &c, const bool decoupled)const{

  if (!decoupled){

	TRIPS trips;
	vector<double> rhs;
	trips.reserve(geom_con.size()*3+vol_self_con.size()*15+surface_self_con.size()*12);
	rhs.reserve(geom_con.size()+vol_self_con.size()+surface_self_con.size());

	for(size_t i = 0; i < geom_con.size(); i++)
	  geom_con[i].getConstraints(trips, rhs, linear_con);

	for(size_t i = 0; i < vol_self_con.size(); i++)
	  vol_self_con[i].getConstraints(trips, rhs, linear_con);

	for(size_t i = 0; i < surface_self_con.size(); i++)
	  surface_self_con[i].getConstraints(trips, rhs, linear_con);

	const int num_var = getLinearCon().size()*3;
	A.setZero();
	A.resize(rhs.size(), num_var);
	A.reserve(trips.size());
	A.setFromTriplets( trips.begin(), trips.end() );
	A.makeCompressed();
  
	c.resize(rhs.size());
	for (int i = 0; i < c.size(); ++i){
	  c[i] = rhs[i];
	}

  }else{
	
	MATH::convert(getLinearCon(), A, c);
  }
}

void LinearConCollider::print()const{
  
  INFO_LOG("static friction: "<< friction_s);
  INFO_LOG("kinetic friction: "<< friction_k);
}
