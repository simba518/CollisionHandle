#include <boost/filesystem.hpp>
#include <Eigen/Geometry> 
#include <JsonFilePaser.h>
#include <assertext.h>
#include <FEMGeom.h>
#include <FEMReducedSystem.h>
#include <FEMRigidReducedSystem.h>
#include <FEMWarppedReducedSystem.h>
#include <FEMMesh.h>
#include <ObjMesh.h>
#include <MakeMesh.h>
#include "Simulator.h"
using namespace UTILITY;
USE_PRJ_NAMESPACE

void Simulator::init(const string &json_file){

  // open init file
  init_file_name = json_file;
  JsonFilePaser jsonf;
  if ( !jsonf.open(json_file) ){
	return;
  }
  
  // load mesh
  {
	DEBUG_LOG("load mesh");
	vector<string> abq_file;
	if( jsonf.readFilePath("vol_file", abq_file) ){
	  assert_ge(abq_file.size(),1);
	  fem_solver->getMesh().reset(abq_file[0],0.0f);
	}
	vector<vector<double> > trans_rot_scale;
	if(jsonf.read("trans_rot_scale", trans_rot_scale)){
	  for(size_t i = 0; i < trans_rot_scale.size(); i++)
		fem_solver->getMesh().applyTrans(transRotScaleToMat(trans_rot_scale[i]),i,true,true);
	}
  }
  
  // set material
  {
	DEBUG_LOG("set material");
	fem_solver->getMesh().getB(0)._system.reset(new FEMSystem(fem_solver->getMesh().getB(0)));
	FEMSystem& sys=*(fem_solver->getMesh().getB(0)._system);
	sys.clearEnergy();
	vector<string> elastic_mtl_file;
	if( jsonf.readFilePath("elastic_mtl", elastic_mtl_file) ){
	  assert_ge(elastic_mtl_file.size(),1);
	  sys.readEnergy(elastic_mtl_file[0], MaterialEnergy::COROTATIONAL, true);
	}else{
	  sys.addEnergyMaterial(250000.0f, 0.45f, MaterialEnergy::COROTATIONAL, true);
	}

	double ak = 0.01, am = 0.0;
	jsonf.read("alpha_k",ak,0.01);
	jsonf.read("alpha_m",am,0.0);
	assert_ge(ak,0.0f);
	assert_ge(am,0.0f);
	sys.setDamping(ak, am);

	vector<double> gravity;
	jsonf.read("gravity",gravity);
	assert_eq(gravity.size(),3);
	sys.addEnergyMass( Vec3(gravity[0], gravity[1], gravity[2]), NULL);

	double coll_k = 1e5;
	jsonf.read("coll_pen", coll_k, 1e5);
	assert_gt(coll_k,0.0f);
	fem_solver->setCollK(coll_k);

	double mu_s = 0.3, mu_k = 0.3;
	jsonf.read("static_friction", mu_s, 0.3);
	jsonf.read("kinetic_friction", mu_k, 0.3);
	fem_solver->setFriction(mu_s, mu_k);
  }

  // set geom
  {
	DEBUG_LOG("set geom");
	boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );

	// obj mesh
  	vector<string> geom_file;
  	if( jsonf.readFilePath("scene", geom_file) ){

	  vector<vector<double> > scene_trans_scale;
	  jsonf.read("scene_trans_scale", scene_trans_scale);
	  double depth = 0.0;
	  jsonf.read("coll_dectect_depth", depth, 0.0);
	  vector<bool> revert_scene_norm;
	  jsonf.read("revert_scene_norm", revert_scene_norm);

	  for(size_t i = 0; i < geom_file.size(); i++){
		Mat4 R = Mat4::Identity();
		bool revert = false;
		if(scene_trans_scale.size() > i){
		  assert_eq(scene_trans_scale[i].size(), 4);
		  R *= scene_trans_scale[i][3];
		  R(3,3) = 1;
		  R(0,3) = scene_trans_scale[i][0];
		  R(1,3) = scene_trans_scale[i][1];
		  R(2,3) = scene_trans_scale[i][2];
		}
		if(revert_scene_norm.size() > i){
		  revert = revert_scene_norm[i];
		}
		geom->addGeomMesh( R , geom_file[i], depth, revert );
	  }
  	}
	
	// plane
	vector<vector<double> > plane;
	if( jsonf.read("plane_ori_trans_scale", plane) ){
	  Vector4d p;
	  Vector3d trans;
	  for(size_t i = 0; i < plane.size(); i++){
		assert_eq(plane[i].size(),7);
		p << plane[i][0], plane[i][1], plane[i][2], 0;
		trans << plane[i][3], plane[i][4], plane[i][5];
		addPlane(geom, p, trans, plane[i][6]);
	  }
	}

	// sphere
	vector<vector<double> > centre_rad;
	if( jsonf.read("sphere", centre_rad) ){
	  Vector3d ctr;
	  for(size_t i = 0; i < centre_rad.size(); i++){
		assert_eq(centre_rad[i].size(),4);
		ctr << centre_rad[i][0], centre_rad[i][1], centre_rad[i][2];
		geom->addGeomSphere(ctr, centre_rad[i][3]);
	  }
	}
	
	// box
	vector<vector<double> > box_ctr_ext;
	if( jsonf.read("box_ext_trans", box_ctr_ext) ){

	  Vector3d ext;
	  Matrix4d trans = Matrix4d::Identity();
	  for(size_t i = 0; i < box_ctr_ext.size(); i++){
		assert_eq(box_ctr_ext[i].size(),6);
		ext << box_ctr_ext[i][0], box_ctr_ext[i][1], box_ctr_ext[i][2];
		trans(0,3) = box_ctr_ext[i][3];
		trans(1,3) = box_ctr_ext[i][4];
		trans(2,3) = box_ctr_ext[i][5];
		geom->addGeomBox(trans, ext);
	  }
	}

	// cylinder
	vector<vector<double> > cylinder;
	if( jsonf.read("cylinder_rad_y_ori_trans", cylinder) ){
	  Vector4d orient;
	  Vector3d trans;
	  for(size_t i = 0; i < cylinder.size(); i++){
		assert_eq(cylinder[i].size(),8);
		const double rad = cylinder[i][0];
		const double y = cylinder[i][1];
		orient << cylinder[i][2], cylinder[i][3], cylinder[i][4], 0.0;
		trans << cylinder[i][5], cylinder[i][6], cylinder[i][7];
		addCylinder(geom, rad, y, orient, trans, 50, 50);
	  }
	}

	// stair
	vector<double > stair_para;
	if( jsonf.read("stair_ext_trans_yd_zd_n", stair_para) ){
	  assert_eq(stair_para.size(), 9);
	  Vector3d ext;
	  Vector3d trans;
	  ext << stair_para[0], stair_para[1], stair_para[2];
	  trans << stair_para[3], stair_para[4], stair_para[5];
	  const double y_diff = stair_para[6];
	  const double z_diff = stair_para[7];
	  const int num = (int)stair_para[8];
	  addStair(geom, ext, trans, y_diff, z_diff, num);
	}

	// assemble
	geom->assemble();
	fem_solver->_geom = geom;
  }

  // set other parameters
  {
	DEBUG_LOG("set other parameters");
	jsonf.readFilePath("save_to", save_results_to,false);
	jsonf.read("h",time_step,0.01);
	assert_gt(time_step, 0.0f);
	jsonf.read("num_frames", total_frames);
	assert_ge(total_frames, 0);

	int newton_max_iteration = 10, mprgp_max_iteration = 100;
	double newton_tolerance = 1e-4, mprgp_tolerance = 1e-4;

	jsonf.read("newton_max_iteration", newton_max_iteration, 10);
	jsonf.read("newton_tolerance", newton_tolerance, 1e-4);
	jsonf.read("mprgp_max_iteration", mprgp_max_iteration, 100);
	jsonf.read("mprgp_tolerance", mprgp_tolerance, 1e-4);
	
	fem_solver->resetImplicitEuler(newton_tolerance, newton_max_iteration);
	fem_solver->setLinearSolverParameters(mprgp_tolerance, mprgp_max_iteration);

	bool enable_self_con = false;
	jsonf.read("enable_self_con", enable_self_con, false);
	fem_solver->setSelfColl(enable_self_con);
	
	vector<vector<double> > init_vel;
	if( jsonf.read("init_vel", init_vel) ){
	  assert_ge(init_vel.size(), 1);
	  assert_eq(init_vel[0].size(), 3);
	  Vector3d vel;
	  vel << init_vel[0][0], init_vel[0][1], init_vel[0][2];
	  fem_solver->setVel(vel);
	}
  }
}

void Simulator::run(){

  boost::filesystem::create_directory(saveResultsTo());
  boost::filesystem::create_directory(saveResultsTo()+"/binary/");

  fem_solver->getMesh().writeVTK(saveResultsTo()+"/mesh.vtk");
  fem_solver->_geom->writeVTK(saveResultsTo()+"/scene.vtk");

  std::ofstream scene_file(saveResultsTo()+"/binary/scene.b", ios::binary);
  fem_solver->_geom->write(scene_file);

  for (size_t frame = 0; frame < totalFrames(); ++frame){
    
	cout << "step: " << frame << endl;

	ostringstream ossm_bin;
  	ossm_bin << saveResultsTo() << "/binary/frame_" << frame << ".b";
	std::ofstream mesh_file(ossm_bin.str(), ios::binary);
	fem_solver->getMesh().write(mesh_file);

	ostringstream ossm_vtk;
  	ossm_vtk << saveResultsTo() << "/frame_" << frame << ".vtk";
  	fem_solver->getMesh().writeVTK( ossm_vtk.str());

  	fem_solver->advance( timeStep() );
  }
}

void Simulator::print()const{
  
  fem_solver->print();
  INFO_LOG("time step: "<< timeStep());
  INFO_LOG("total frames: "<< totalFrames());
  INFO_LOG("save results to: "<< saveResultsTo());
  INFO_LOG("init file:" << init_file_name)
}

void Simulator::addStair(boost::shared_ptr<FEMGeom> geom, const Vector3d &ext, Vector3d trans,
						 const double y_diff, const double z_diff, const int num)const{
  
  assert_ge(num, 0);
  assert(geom);

  Matrix4d trans_m = Matrix4d::Identity();
  trans_m(0,3) = trans[0];
  for(int i = 0; i < num; i++){
	trans_m(1,3) = trans[1]+y_diff*(double)i;
	trans_m(2,3) = trans[2]+z_diff*(double)i;
	geom->addGeomBox(trans_m, ext);
  }
}

void Simulator::addPlane(boost::shared_ptr<FEMGeom> geom, const Vector4d &plane,
						 const Vector3d &trans, const double scale)const{
  
  Vector3d ext;
  ext << 20,2,20;
  ext *= scale;
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) = trans;
  geom->addGeomBox(T*orientToTrans(plane), ext);
}

void Simulator::addCylinder(boost::shared_ptr<FEMGeom> geom, 
							const double rad, const double y,
							const Vector4d &orient, const Vector3d &trans,
							const int slice, const int sliceY)const{
  
  ObjMesh mesh;
  MakeMesh::makeCylinder3D(mesh, rad, y, slice, sliceY, true);
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) = trans;
  geom->addGeomMesh(T*orientToTrans(orient), mesh);
}

Matrix4d Simulator::orientToTrans(const Vector4d &orient)const{
  
  const double alpha=-orient[3]/orient.block<3,1>(0,0).squaredNorm();
  const Vector3d p0=orient.block<3,1>(0,0)*alpha;
  Eigen::Quaterniond q;
  q.setFromTwoVectors(Vector3d::Unit(1),orient.block<3,1>(0,0).normalized());

  Matrix4d trans = Matrix4d::Identity();
  trans.block<3,1>(0,3)=p0-orient.block<3,1>(0,0).normalized();
  trans.block<3,3>(0,0)=q.toRotationMatrix();
  return trans;
}

Matrix4d Simulator::transRotScaleToMat(const vector<double> &trans_rot_scale)const{

  assert_eq(trans_rot_scale.size(), 7);

  // scale
  Matrix4d S = Matrix4d::Identity();
  S.block<3,3>(0,0) *= trans_rot_scale[6];

  // rotation
  Matrix4d R = Matrix4d::Identity();
  const double angle_x = trans_rot_scale[3]*M_PI/180.0;
  const double angle_y = trans_rot_scale[4]*M_PI/180.0;
  const double angle_z = trans_rot_scale[5]*M_PI/180.0;
  Matrix3d m;
  m = AngleAxisd(angle_x, Vector3d::UnitX())
	* AngleAxisd(angle_y,  Vector3d::UnitY())
	* AngleAxisd(angle_z, Vector3d::UnitZ());
  R.block<3,3>(0,0) = m;

  // translation
  Matrix4d T = Matrix4d::Identity();
  T.block<3,1>(0,3) << trans_rot_scale[0], trans_rot_scale[1], trans_rot_scale[2];

  // assemble
  return (T*R*S);
}
