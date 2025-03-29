#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"
#include "Math/Vector4D.h"
#include "TStyle.h"
 
using namespace ROOT;
using namespace ROOT::VecOps;
using RNode = ROOT::RDF::RNode;
using namespace std;

float InvariantMass(float pt1, float eta1, float phi1, float mass1, float pt2, float eta2, float phi2, float mass2){
  float Mij;
  float px1 = pt1*cos(phi1);
  float py1 = pt1*sin(phi1);
  float pz1 = pt1*sinh(eta1);

  float E_1 = sqrt(pow(px1,2) + pow(py1,2) + pow(pz1,2));

  float px2 = pt2*cos(phi2);
  float py2 = pt2*sin(phi2);
  float pz2 = pt2*sinh(eta2);

  float E_2 = sqrt(pow(px2,2)+pow(py2,2)+pow(pz2,2));
  float E = E_1+E_2;
  
  Mij = sqrt(pow(E,2) - (pow(px1+px2,2) + pow(py1+py2,2) + pow(pz1+pz2,2)));

  return Mij;
}

float Ratio_pt(float pt1,float pt2){
	if (pt1 == 0 || pt2 == 0){
		std::cerr << "Error: Division by zero in Ratio_pt. pt1 = " << pt1 << ", pt2 = " << pt2 << std::endl;
	}

	if (pt1 > pt2){
		return pt2/pt1;
	}
	else{
		return pt1/pt2;
	}
}
