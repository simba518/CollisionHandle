#include <boost/filesystem.hpp>
#include <JsonFilePaser.h>
#include <assertext.h>
#include <FEMGeom.h>
#include <FEMReducedSystem.h>
#include <FEMRigidReducedSystem.h>
#include <FEMWarppedReducedSystem.h>
#include <FEMMesh.h>
#include <ObjMesh.h>
#include "Simulator.h"
using namespace UTILITY;

void Simulator::init(const string &json_file){
  
  // open init file
  JsonFilePaser jsonf;
  if ( !jsonf.open(json_file) ){
	return;
  }
  
  // load mesh
  {
	string abq_file;
	if( jsonf.readFilePath("vol_file", abq_file) )
	  fem_solver->getMesh().reset(abq_file,0.0f);
  }
  
  // set material
  {
	double ak, am;
	jsonf.read("alpha_k",ak,0.0);
	jsonf.read("alpha_m",am,0.0);

	fem_solver->getMesh().getB(0)._system.reset(new FEMSystem(fem_solver->getMesh().getB(0)));
	FEMSystem& sys=*(fem_solver->getMesh().getB(0)._system);
	sys.clearEnergy();
	sys.addEnergyMaterial(250000.0f, 0.45f, MaterialEnergy::COROTATIONAL, true);

	vector<double> gravity;
	jsonf.read("gravity",gravity);
	assert_eq(gravity.size(),3);
	sys.addEnergyMass( Vec3(gravity[0], gravity[1], gravity[2]), NULL);

	double coll_k = 1e5;
	jsonf.read("coll_pen", coll_k);
	assert_gt(coll_k,0.0f);
	fem_solver->setCollK(coll_k);
  }

  // set geom
  {
  	ObjMesh m;
  	string geom_file;
  	if( jsonf.readFilePath("scene", geom_file) ){
	  std::ifstream mfile(geom_file);
	  m.read(mfile);
	  m.smooth();
	  m.makeUniform();
  	  boost::shared_ptr<FEMGeom> geom( new FEMGeom(3) );

	  Mat4 R = Mat4::Identity();
	  double scene_scale = 1.0f;
	  if(jsonf.read("scene_scale", scene_scale)){
		R *= scene_scale;
		R(3,3) = 1;
	  }
	  vector<double> scene_trans;
	  if(jsonf.read("scene_trans", scene_trans)){
		assert_eq(scene_trans.size(),3);
		R(0,3) = scene_trans[0];
		R(1,3) = scene_trans[1];
		R(2,3) = scene_trans[2];
	  }

	  double depth = 0.0;
	  jsonf.read("coll_dectect_depth", depth, 0.0);
	  bool revert_scene_norm = false;
	  jsonf.read("revert_scene_norm", revert_scene_norm, false);
  	  geom->addGeomMesh( R , m, depth, revert_scene_norm);
  	  geom->assemble();
  	  fem_solver->_geom = geom;

  	}
  }

  // set target directory, time step, as well as number of frames.
  {
	jsonf.readFilePath("save_to", save_results_to,false);
	jsonf.read("h",time_step,0.01);
	assert_gt(time_step, 0.0f);
	jsonf.read("num_frames", total_frames);
	assert_gt(total_frames, 1);
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
