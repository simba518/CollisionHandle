#include <iomanip>
#include <assertext.h>
#include <FEMMesh.h>
#include "LinearConCollider.h"
USE_PRJ_NAMESPACE

void LinearConCollider::handle(boost::shared_ptr<FEMBody> b[5],boost::shared_ptr<FEMVertex> v[5],const Vec3d coef[5],sizeType nrV){

  /// @todo undefined function.
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
  assert_in(vert_id, 0, geom_con.size());
  assert_eq_ext(plane, plane, "n:" << n.transpose());
  addConPlane(geom_con[vert_id], plane);
  
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
