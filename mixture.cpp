#include <fstream>
#include <iostream>
#include <string>
#include <random>

int main(int argc, char *argv[])
{
  if(argc<2){
    return 1;
  }
  std::string dirname = "data/";
  std::vector<std::string> filename = {"sca-0_single",
				       "sca-1_single",
				       "sca-2_single",
				       "sca-3_single",
				       "sca-4_single",
				       "sca-5_single",
				       "sca-6_single",
				       "sca-7_single",
				       "sca-8_single"
  };
  std::vector<std::vector<std::vector<double>>> flush;
  std::vector<std::vector<std::vector<double>>> idealvalue;
  std::vector<std::vector<std::vector<double>>> teachervalue;
  std::vector<std::vector<double>> flush_temp;
  std::vector<std::vector<double>> idealvalue_temp;
  std::vector<std::vector<double>> teachervalue_temp;
  std::vector<double> value_temp;
  double value;

  for(auto itr=filename.begin();itr!=filename.end();++itr){
    std::ifstream fflush(dirname+(*itr)+"_flush.dat");
    std::ifstream fideal(dirname+(*itr)+"_idealvalue.dat");
    std::ifstream fteacher(dirname+(*itr)+"_teachervalue.dat");
    flush_temp.clear();
    idealvalue_temp.clear();
    teachervalue_temp.clear();
    while(!fflush.fail()){
      value_temp.clear();
      for(int ii=0;ii<2*1024*256;ii++){
	fflush >> value;
	value_temp.push_back(value);
      }
      flush_temp.push_back(value_temp);
      value_temp.clear();
      for(int ii=0;ii<8;ii++){
	fideal >> value;
	value_temp.push_back(value);
      }
      idealvalue_temp.push_back(value_temp);
      value_temp.clear();
      for(int ii=0;ii<11;ii++){
	fteacher >> value;
	value_temp.push_back(value);
      }
      teachervalue_temp.push_back(value_temp);
    }
    flush.push_back(flush_temp);
    idealvalue.push_back(idealvalue_temp);
    teachervalue.push_back(teachervalue_temp);
    fflush.close();
    fideal.close();
    fteacher.close();
  }

  std::random_device seed_gen;
  std::mt19937 engine(seed_gen());
  std::vector<std::uniform_int_distribution<>> dist;
  std::uniform_real_distribution<> rndm(0,1);
  std::vector<bool> add(filename.size());
  for(int ii=0;ii<filename.size();ii++){
    dist.push_back(std::uniform_int_distribution<>(0,flush[ii].size()-1));
  }

  int num = 0;
  std::string outname = argv[1];
  std::ofstream oftot(dirname+outname+"_tot.dat");
  std::ofstream offlush(dirname+outname+"_flush.dat");
  std::ofstream ofideal(dirname+outname+"_idealvalue.dat");
  std::ofstream ofteach(dirname+outname+"_teachervalue.dat");
  double threshold = 1;
  std::vector<int> eve_id(filename.size());
  
  while(num<10000){
    int add_num = 0;
    for(int ii=0;ii<filename.size();ii++){
      if(rndm(engine)>0.8){
	add[ii] = true;
	add_num++;
      }else{
	add[ii] = false;
      }
    }
    if(add_num<3){
      for(int kk=0;kk<2*1024*256;kk++){
	value = 0;
	for(int jj=0;jj<filename.size();jj++){
	  eve_id[jj] = dist[jj](engine);
	}
	for(int jj=0;jj<filename.size();jj++){
	  if(add[jj]){
	    value += flush[jj][eve_id[jj]][kk];
	  }else{
	    continue;
	  }
	}
	if(value>threshold){
	  oftot << 1 << " ";
	}else{
	  oftot << 0 << " ";
	}
	offlush << value << " ";
      }
      oftot << std::endl;
      offlush << std::endl;
      for(int ii=0;ii<filename.size();ii++){
	if(add[ii]){
	  for(int jj=0;jj<idealvalue[ii][eve_id[ii]].size();jj++){
	    ofideal << idealvalue[ii][eve_id[ii]][jj] << " ";
	  }
	  for(int jj=0;jj<teachervalue[ii][eve_id[ii]].size();jj++){
	    ofteach << teachervalue[ii][eve_id[ii]][jj] << " ";
	  }
	}
      }
      ofideal << std::endl;
      ofteach << std::endl;
      num++;
    }
    std::cout << "num = " << num << "\r" << std::flush;
  }

  oftot.close();
  offlush.close();
  ofideal.close();
  ofteach.close();
    
  return 0;
}
