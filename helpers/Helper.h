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
  float px_1 = pt_1*cos(phi1);
  float py_1 = pt_1*sin(phi1);
  float pz_1 = pt_1*sinh(phi1);

  float E_1 = sqrt(px_1**2 + py_1**2 + pz_1**2);

  float px_2 = pt_2*cos(phi2);
  float py_2 = pt_2*sin(phi2);
  float pz_2 = pt_2*sinh(phi2);

  float E_2 = sqrt(px_2**2+py_2**2+pz_2**2);
  float E = E_1+E_2;
  
  Mij = sqrt((E)**2 - ((px_1+px_2)**2 + (py_1+py_2)**2 + (pz_1+pz_2)**2))

  return Mij;
}
