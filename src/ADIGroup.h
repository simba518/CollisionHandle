#ifndef _ADIGROUP_H_
#define _ADIGROUP_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <tuple>
using namespace std;

typedef vector<set<int> > VVI;
typedef vector<VVI> VVVI;

/**
 * @class ADIGroup split nodes as groups for ADI preconditioning
 * 
 */
class ADIGroup{
	
public:
  ADIGroup(){
	clear();
  }
  void clear(){
	verts.clear();
	groups.clear();
	groups.resize(3);
  }
  bool loadMesh(const string &filename){

	clear();
	ifstream in(filename.c_str());
	if( !in.is_open() ){
	  cout << "ERROR: failed to open file: "<< filename << endl;
	  return false;
	}

	string line;
	while( getline(in,line) && ((line.find("NODE")==string::npos)) ){};

	char tc;
	double v[3];
	while (in.good() && !in.eof()){
	  in >> line;
	  if ( line.find("ELEMENT") != string::npos ){
		break;
	  }
	  in >> v[0] >> tc >> v[1]>> tc >> v[2];
	  verts.push_back(v[0]);
	  verts.push_back(v[1]);
	  verts.push_back(v[2]);
	}

	assert(verts.size()%3 == 0);
	cout << "num verts: " << verts.size()/3 << endl;

	const bool succ = in.good();
	in.close();
	return succ;
  }
  bool save(const string &filename)const{
	
	ofstream out(filename.c_str());
	if( !out.is_open() ){
	  cout << "ERROR: failed to open file: "<< filename << endl;
	  return false;
	}

	assert(groups.size()==3);
	const string label[3] = {string("x direction"),string("y direction"),string("z direction")};
	for (int d = 0; d < 3; d++){
	  out << label[d] << endl;
	  out << "groups " << groups[d].size() << endl;
	  for (size_t g = 0; g < groups[d].size(); ++g){
		out << "len " << groups[d][g].size() << endl;
		set<int>::const_iterator it = groups[d][g].begin();
		for (; it != groups[d][g].end(); ++it){
		  out << (*it) << " ";
		}
		out << endl;
	  }
	}

	return out.good();
	
  }
  void split(const int num_gx, const int num_gy, const int num_gz){

	groups.clear();
	groups.resize(3);
	split(groups[0], num_gx, 0);
	split(groups[1], num_gy, 1);
	split(groups[2], num_gz, 2);
	assert( checkGroups() );
  }
  double min(int direction)const{
	assert(direction >= 0 && direction <= 2);
	double v = 1e9;
	for (size_t i = direction; i < verts.size(); i+=3){
	  v = std::min(v, verts[i]);
	}
	return v;
  }
  double max(int direction)const{
	assert(direction >= 0 && direction <= 2);
	double v = -1e9;
	for (size_t i = direction; i < verts.size(); i+=3){
	  v = std::max(v, verts[i]);
	}
	return v;
  }
  bool checkGroups()const{
	
	bool valid = true;
	for (int d = 0; d < 3 && valid; d++){

	  set<int> all_verts_id;
	  for (size_t g = 0; g < groups[d].size() && valid; ++g){
		set<int>::const_iterator it = groups[d][g].begin();
		for (; it != groups[d][g].end() && valid; ++it){
		  valid = (all_verts_id.find(*it) == all_verts_id.end());
		  all_verts_id.insert(*it);
		  if(!valid){
			cout << "ERROR: nodes are duplicated\n";
		  }
		}
	  }
	  if ( all_verts_id.size()*3 != verts.size() ){
		cout << "ERROR: some nodes are missed\n";
	  }
	  valid &= (all_verts_id.size()*3 == verts.size());
	}
	return valid;
  }

  static bool loadAdiGroups(const string filename, vector<vector<set<int> > > &groups){

	ifstream in(filename.c_str());
	if( !in.is_open()){
	  return false;
	}

	groups.clear();
	groups.resize(3);

	string tempt;
	int num_g = 0;
	int len = 0;
	int vid = 0;
	for (int direction = 0; direction < 3; ++direction){

	  in >> tempt >> tempt >> tempt >> num_g;
	  assert(num_g > 0);
	  groups[direction].resize(num_g);
	  for (int g = 0; g < num_g; ++g){
		in >> tempt >> len;
		assert(len >= 0);
		for (int i = 0; i < len; ++i){
		  in >> vid;
		  groups[direction][g].insert(vid);
		}
	  }
	}

	in.close();
	return in.good();
  }

protected:
  void split(VVI &group, const int num_g, const int direction)const{
	
	assert(num_g > 0);
	assert(direction >= 0 && direction <= 2);
	
	const double min_x = min(direction)-1e-9;
	const double max_x = max(direction)+1e-9;
	assert(min_x < max_x);
	const double step = (max_x-min_x)/num_g;
	
	group.clear();
	for (int i = 0; i < num_g; ++i){
	  set<int> s;
	  addToSet( step*i+min_x, step*i+step+min_x, direction, s);
	  eraseExisted(group, s);
	  group.push_back(s);
	}	
  }
  void eraseExisted(const VVI &group, set<int> &s)const{
	
	for (size_t g = 0; g < group.size(); ++g){
	  set<int>::const_iterator it = group[g].begin();
	  for ( ;it != group[g].end(); ++it )
		s.erase(*it);
	}
  }
  void addToSet(const double begin,const double end,const int direction,set<int>&s)const{

	for (size_t i = direction; i < verts.size(); i += 3){
	  if (verts[i] >= begin && verts[i] <= end){
		s.insert(i/3);
	  }
	}
  }
	
private:
  VVVI groups;
  vector<double> verts;
};
 
#endif /* _ADIGROUP_H_ */
