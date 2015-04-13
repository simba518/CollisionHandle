#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

bool GeomConCache::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n, VVVec4d &linear_con){

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
	this->vert_id = vert_id;
	this->plane_id = (int)linear_con[vert_id].size()-1;
	normal = plane.segment<3>(0);
	return true;
  }

  return false;
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

bool SurfaceSelfConCache::handle(boost::shared_ptr<FEMBody> body[5], boost::shared_ptr<FEMVertex> v[5],
								 const Vec3d coef[5], sizeType nrV, VVVec4d &linear_con){

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

  return true;
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

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b, const double ground_y, const double delta_y){
  
  Vec3d n;
  n << 0,1,0;
  Vector4d plane;
  plane.segment<3>(0) = n;
  plane[3] = -ground_y;

  for (sizeType i = 0; i < b->nrV(); ++i){
	boost::shared_ptr<FEMVertex> v = b->getVPtr(i);
	if (v->_pos[1] - delta_y <= ground_y){
	  const int vert_id = b->_offset/3 + v->_index;
	  const int added = GeomConCache::addConPlane(linear_con[vert_id], plane);
	  if (added >= 0){
		continue;
	  }
	  const int plane_id = linear_con[vert_id].size()-1;
	  if (!decouple_constraints){
		geom_con.push_back( GeomConCache(vert_id, plane_id, n) );
	  }else if ( !collided(vert_id) ){
		geom_con.push_back( GeomConCache(vert_id, plane_id, n) );
		coll_as_vert[vert_id] = true;
	  }
	}
  }
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b,boost::shared_ptr<FEMVertex> v,const Vec3& n){

  if (!decouple_constraints){
	GeomConCache one_con;
	if( one_con.handle(b, v, n, linear_con) )
	  geom_con.push_back(one_con);
  }else{
	const int vert_id = b->_offset/3 + v->_index;
	if( !collided(vert_id) ){
	  GeomConCache one_con;
	  if( one_con.handle(b, v, n, linear_con) ){
		geom_con.push_back(one_con);
		coll_as_vert[vert_id] = true;
	  }
	}
  }
}

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5], const Vec3d coef[5],sizeType nrV) {

  if (!decouple_constraints){
	SurfaceSelfConCache one_con;
	if ( one_con.handle(b, v, coef, nrV, linear_con) )
	  surface_self_con.push_back(one_con);
  }else {
	const int i = b[0]->_offset/3 + v[0]->_index;
	const int j = b[1]->_offset/3 + v[1]->_index;
	const int k = b[1]->_offset/3 + v[2]->_index;
	const int l = b[1]->_offset/3 + v[3]->_index;
	const bool coll = (collided(i))|| (collided(j))|| (collided(k))|| (collided(l));
	if(!coll){
	  SurfaceSelfConCache one_con;
	  if ( one_con.handle(b, v, coef, nrV, linear_con) ){
		surface_self_con.push_back(one_con);
		coll_as_vert[i] = true;
		coll_as_face[j] = true;
		coll_as_face[k] = true;
		coll_as_face[l] = true;
	  }
	}
  }
}

void LinearConCollider::getConstraints(SparseMatrix<double> &A, VectorXd &c)const{

	TRIPS trips;
	vector<double> rhs;
	trips.reserve(geom_con.size()*3+surface_self_con.size()*12);
	rhs.reserve(geom_con.size()+surface_self_con.size());

	for(size_t i = 0; i < geom_con.size(); i++)
	  geom_con[i].getConstraints(trips, rhs, linear_con);

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
}

void LinearConCollider::print()const{
  
  INFO_LOG("decouple constraints: "<< (decouple_constraints ? "true" : "false"));
}

void ContinueCollider::handle(boost::shared_ptr<ClothMesh::ClothVertex> V1,
							  boost::shared_ptr<ClothMesh::ClothTriangle> T2,
							  const Vec3d normal,const Vec4d& omg,scalarD t){

  // ClothCollision::CollisionHandler::handle(V1,T2,normal,omg,t);
  Vec3d n = T2->getNormal().normalized();
  if (n.dot(normal) >= 0){
	n = normal;
  }
  if(ClothMesh::CLOTH_MESH == V1->_type){
	const int body_1 = V1->body_index;
	const int vert_1 = V1->node_id_on_body;
	if(ClothMesh::CLOTH_MESH == T2->_type){
	  SurfaceSelfConCache one_con;
	  one_con.i = V1->node_id_on_body;
	  one_con.j = T2->getV0()->node_id_on_body;
	  one_con.k = T2->getV1()->node_id_on_body;
	  one_con.l = T2->getV2()->node_id_on_body;
	  const bool coll = (collided(one_con.i))|| (collided(one_con.j))
		|| (collided(one_con.k))|| (collided(one_con.l));
	  if( (!decouple_constraints) || (!coll) ){
		one_con.a = -omg[1];
		one_con.b = -omg[2];
		one_con.c = -omg[3];
		one_con.n = n;
		surface_self_con.push_back(one_con);
		coll_as_vert[one_con.i] = true;
		coll_as_face[one_con.j] = true;
		coll_as_face[one_con.k] = true;
		coll_as_face[one_con.l] = true;
	  }
	}else{
	  const int old_size = (int)geom_con.size();
	  boost::shared_ptr<FEMBody> b = femmesh.getBPtr(body_1);
	  boost::shared_ptr<FEMVertex> v = b->getVPtr(vert_1);
	  LinearConCollider::handle(b,v,n);
	  if (old_size < (int)geom_con.size()){
		// assert_eq(geom_con[old_size].normal,n);
		assert_ge(linear_con[vert_1].size(),1);
		const Vec3d &v1 = T2->getV0()->_pos;
		const Vec3d &v2 = T2->getV1()->_pos;
		const Vec3d &v3 = T2->getV2()->_pos;
		const double a = omg[1];
		const double b = omg[2];
		const double c = omg[3];
		const int k = linear_con[vert_1].size()-1;
		linear_con[vert_1][k][3] = n.dot(v1)*a + n.dot(v2)*b + n.dot(v3)*c;
	  }
	}
  }
}

void ContinueCollider::handle(boost::shared_ptr<ClothMesh::ClothEdge> E1,
							  boost::shared_ptr<ClothMesh::ClothEdge> E2,
							  const Vec3d normal,const Vec4d& omg,scalarD t){

  Vec3d n = E2->getNormal().normalized();
  if(n.dot(normal) > 0){
	n = normal;
  }

  if(ClothMesh::CLOTH_MESH == E1->_type){
	if( (ClothMesh::CLOTH_MESH == E2->_type) ){
	  SurfaceSelfConCache one_con;
	  one_con.i = E1->_v[0]->node_id_on_body;
	  one_con.j = E1->_v[1]->node_id_on_body;
	  one_con.k = E2->_v[0]->node_id_on_body;
	  one_con.l = E2->_v[1]->node_id_on_body;
	  const bool coll = (collided(one_con.i))|| (collided(one_con.j))
	  	|| (collided(one_con.k))|| (collided(one_con.l));
	  if( (!decouple_constraints) || (!coll) ){
		assert_gt(abs(omg[0]),1e-8);
		one_con.a = -(omg[1]/omg[0]);
		one_con.b = -(omg[2]/omg[0]);
		one_con.c = -(omg[3]/omg[0]);
		one_con.n = n;
		surface_self_con.push_back(one_con);
		coll_as_vert[one_con.i] = true;
		coll_as_face[one_con.j] = true;
		coll_as_face[one_con.k] = true;
		coll_as_face[one_con.l] = true;
	  }
	}
  }
}

void ContinueCollider::init(){

  for (int i = 0; i < deform_surface.size(); ++i){
    CCD.delMesh(deform_surface[i]);
  }
  for (int i = 0; i < scene.size(); ++i){
    CCD.delMesh(scene[i]);
  }

  deform_surface.clear();
  volNodeId_clothId_vertId.clear();
  scene.clear();
  for (int bi = 0; bi < femmesh.nrB(); ++bi){
    const FEMBody &body = femmesh.getB(bi);
	ObjMeshD obj;
	body.writeObj(obj);
	boost::shared_ptr<ClothMesh> cloth = boost::shared_ptr<ClothMesh>( new ClothMesh(obj,ClothMesh::CLOTH_MESH) );
	deform_surface.push_back(cloth);
	const vector<int>& vMap=obj.getIG();
	Vector3i one_map;
	one_map[1] = bi;
	for (int vert_id = 0, i = 0; i < vMap.size(); ++i){
	  if(vMap[i] >= 0){
		one_map[0] = i;
		one_map[2] = vMap[i];
		volNodeId_clothId_vertId.push_back(one_map);
		assert_eq(cloth->_vss[vMap[i]]->_index, vMap[i]);
		cloth->_vss[vMap[i]]->body_index = bi;
		cloth->_vss[vMap[i]]->node_id_on_body = i;
	  }
	}
  }
  
  for (int i = 0; i<geom->nrG(); ++i ){
	ObjMesh obj;
    geom->getG(i).getMesh(obj);
	scene.push_back(boost::shared_ptr<ClothMesh>(new ClothMesh(obj,ClothMesh::RIGID_MESH)));
  }

  for (int i = 0; i < deform_surface.size(); ++i){
    CCD.addMesh(deform_surface[i]);
  }
  for (int i = 0; i < scene.size(); ++i){
    CCD.addMesh(scene[i]);
  }

}

void ContinueCollider::collide(const VectorXd &last_pos, const VectorXd &cur_pos){

  for (int i = 0; i < (int)volNodeId_clothId_vertId.size(); ++i){
	const int node_id = volNodeId_clothId_vertId[i][0];
	const int cloth_id = volNodeId_clothId_vertId[i][1];
	const int vert_id = volNodeId_clothId_vertId[i][2];
	deform_surface[cloth_id]->_vss[vert_id]->_lastPos = last_pos.segment<3>(node_id*3);
	deform_surface[cloth_id]->_vss[vert_id]->_pos = cur_pos.segment<3>(node_id*3);
  }

  CCD.collide(*this);

}
