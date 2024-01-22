#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include <iostream> // needed for io
#include <cstdio>   // needed for io
#include <fstream>
#include <cmath>
#include <iomanip> // For setw function

using namespace std;
using namespace fastjet;

//------------------------------------------------------------------------
// the user information

class MyInfo : public PseudoJet::UserInfoBase{
public:
  // default ctor
  //  - pdg_id        the PDG id of the particle
  MyInfo(const float status_id,const float pid){
    _status_id=status_id;
    _pid=pid;
  }

  /// access to the PDG id
  float status_id() const { return _status_id;}
  float pid() const {return _pid;}
  
protected:
  float _status_id;         // the associated status id
  float _pid;               //assosciated pid
};

//--------------------------------------------------------
// the custom jet selector
class Jet_Selector: public SelectorWorker{
public:
  // default ctor
  Jet_Selector(){}

  // the selector's description
  string description() const{
    return "Jets with highest pT charged particle more than 5GeV"; 
  }

  //check if the particle is charged
  bool Check_Charged(float pid) const{
    // Charged hadrons: (pi+, K+, p+, Sigma+, Sigma-, Xi-, Omega-) 
    float charged_hadrons_pid[]={211, 321, 2212, 3222, 3112, 3312, 3334};
    for (int i=0;i<7;i++){
      if (charged_hadrons_pid[i]==abs(pid)){
        return true;
      }
    }
    return false;
  }

  // keeps the jets which have highest pt charged particle >5GeV
  bool pass(const PseudoJet &jet) const{
    const vector<PseudoJet> constituents = jet.constituents();
    for (const PseudoJet& element:constituents){
      if (Check_Charged(element.user_info<MyInfo>().pid()) && element.user_info<MyInfo>().status_id()>=0){
        if(element.pt()>=5){
          return true;
        }   
      }
    }
    return false;
  }
};

// the function that allows to write simply
Selector Selector_Criteiron(){
  return Selector(new Jet_Selector());
}
//----------------------------------------------
//Subtracts Holes
void hole_sub(std::vector<PseudoJet> jets){
  for(PseudoJet& jet:jets){   
    for(PseudoJet& p:jet.constituents()){
      if(p.user_info<MyInfo>().status_id()==-1){
        jet.reset(jet.px()-p.px(),jet.py()-p.py(),jet.pz()-p.pz(),jet.e()-p.e());
      }
    }
  }
}

//------------------------------------------------------------
/// MAIN PROGRAM
int FastJetProcedure(int lower_lim,int upper_lim,int run){
  std::string input_file="/wsu/tmp/hl9735/5020_PbPb_30-40_0.30_2.0_1_0.3/JetscapeHadronListBin"+std::to_string(lower_lim)+"_"+std::to_string(upper_lim)+"_"+"Run"+std::to_string(run)+".out";
  std::string output_file="/rs/rs_grp_majshen/tier1/Binary_output_0.3/jets/JetscapeHadronListBin"+std::to_string(lower_lim)+"_"+std::to_string(upper_lim)+"_"+"Run"+std::to_string(run)+"_jets.out";
 
  std::ofstream file1(output_file);
  std::ifstream file0(input_file);  

  if (file0.is_open()){
    std::cout<<"File successfully opened"<<std::endl;
  }
  else{
    std::cout<<"Failed to open the file "<<output_file<<std::endl;
  }

  int i,j,event;
  event=0;
  std::string str;
  float values[9];		// index >>PID>>particle_status>> E >> px >> py >> pz >> eta>> phi

  vector<fastjet::PseudoJet> input_particles;
  double R = 0.2;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  std::string s;
  while(std::getline(file0,str)){                 //std::getline(file0,str)){
    i=0;
    j=0;
    //std::cout<<str<<" "<<str.size()<<std::endl;
    if (str[0]=='#'){
      std::cout<<"New Event "<<event<<std::endl;
      file1<<"Event "<<event<<"\n";
      if(event>0){        //still there is something to be done about the last event 4000
        fastjet::ClusterSequence clust_seq(input_particles,jet_def);
        hole_sub(clust_seq.inclusive_jets());
        Selector jet_selector=Selector_Criteiron();
        vector<PseudoJet> inclusive_jets=jet_selector(clust_seq.inclusive_jets());
        file1<< std::setw(5) << "jet #" << std::setw(15) << "rapidity"<< std::setw(15) << "phi" << std::setw(15) << "pt" << std::endl;

        for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
          file1<< std::setw(5) << i << std::setw(15) << inclusive_jets[i].rap()<< std::setw(15) << inclusive_jets[i].phi() << std::setw(15) << inclusive_jets[i].perp() << std::endl;
        }
        file1<<"\n";
      }
      event+=1;
      input_particles.clear();
      continue;
    }   
    while (i!=str.size()){//(str[i]!='\0'){
      if(str[i]==' '||str[i]=='\t'||i==str.size()-1){
        //std::cout<<s<<" LOVE "<<std::endl;
        values[j]=std::atof(s.c_str());
        s.clear();
        j+=1;
      }   
      else{
        s+=str[i];
      }   
      i+=1;
    }
    //std::cout<<values[4]<<"  "<<values[5]<<"  "<<values[6]<<"  "<<values[3]<<std::endl;
    PseudoJet p(values[4],values[5],values[6],values[3]); //px,py,pz,
    p.set_user_info(new MyInfo(values[2],values[1]));   
    input_particles.push_back(p);
  }   
  
  std::cout<<"New Event "<<event<<std::endl;  
  file1<<"Event "<<event<<"\n";
  fastjet::ClusterSequence clust_seq(input_particles,jet_def);
  hole_sub(clust_seq.inclusive_jets());
  Selector jet_selector=Selector_Criteiron();
  vector<PseudoJet> inclusive_jets=jet_selector(clust_seq.inclusive_jets());
  file1<< std::setw(5) << "jet #" << std::setw(15) << "rapidity"<< std::setw(15) << "phi" << std::setw(15) << "pt" << std::endl;

  for (unsigned int i = 0; i < inclusive_jets.size(); i++) {
          file1<< std::setw(5) << i << std::setw(15) << inclusive_jets[i].rap()<< std::setw(15) << inclusive_jets[i].phi() << std::setw(15) << inclusive_jets[i].perp() << std::endl;
  }
  //file1<<"\n";
  file0.close();  
  file1.close();
  // read in input particles
  //----------------------------------------------------------
  return 0;
}

int main(){
  int pt_bin[33][2]={{1500,2510}, {800,1500}, {450, 800}, {350, 450}, {300, 350}, {280, 300}, {260, 280}, {240, 260}, {220, 240}, {210, 220}, {200, 210}, {190, 200}, {180,190}, {170, 180}, {160, 170}, {150, 160}, {140, 150}, {130, 140}, {120, 130}, {110, 120},  {100, 110}, {90,100},  {80,90}, {70,80}, {60, 70}, {55, 60}, {50, 55},  {45, 50}, {40, 45}, {35, 40}, {30, 35}, {25, 30}, {20,25}};
  for(int i=0;i<255;i++){
    for(int j=0;j<33;j++){
      FastJetProcedure(pt_bin[j][0],pt_bin[j][1],i);
    }
  }
}
