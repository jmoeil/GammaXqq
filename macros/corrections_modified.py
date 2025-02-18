import ROOT
import correctionlib
correctionlib.register_pyroot_binding()

# C++ function to get the raw pT from the corrected pT


def setupjecs(JECfile, corrfile):
    ROOT.gInterpreter.Declare('auto cset_jecs = correction::CorrectionSet::from_file("{}");'.format(JECfile))
    ROOT.gInterpreter.Declare('auto cset_l1jecs = cset_jecs->at("{}");'.format(corrfile+'L1FastJet_AK4PFPuppi'))
    ROOT.gInterpreter.Declare('auto cset_l2jecs = cset_jecs->at("{}");'.format(corrfile+'L2Relative_AK4PFPuppi'))
    ROOT.gInterpreter.Declare('auto cset_l3jecs = cset_jecs->at("{}");'.format(corrfile+'L3Absolute_AK4PFPuppi'))
    ROOT.gInterpreter.Declare('auto cset_residualjecs = cset_jecs->at("{}");'.format(corrfile+'L2L3Residual_AK4PFPuppi'))

    ROOT.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> JetRawPt(const ROOT::VecOps::RVec<float> &pt,
                                   const ROOT::VecOps::RVec<float> &rawf){

    ROOT::VecOps::RVec<float> Jet_rawPt;

    for(unsigned int i=0; i<pt.size(); i++){
       double rawPt = pt[i] * (1 - rawf[i]);
       Jet_rawPt.push_back(rawPt);
    }

    return Jet_rawPt;

    }
    ''')


    
    # C++ function to apply the Jet Energy Corrections 
    ROOT.gInterpreter.Declare('''
    ROOT::VecOps::RVec<float> JetCorPt(const ROOT::VecOps::RVec<float> &area, 
                                               const ROOT::VecOps::RVec<float> &eta,
                                               const ROOT::VecOps::RVec<float> &phi,
                                               const ROOT::VecOps::RVec<float> &pt,
                                               const ROOT::VecOps::RVec<float> &rawf,
                                               const float &rho,
    const bool &isData){

    ROOT::VecOps::RVec<float> Jet_corPt;
    bool debug = false;
    //auto cset = correction::CorrectionSet::from_file(JECfile);
//auto myelement = cset->at(corrfile + "L3Absolute_AK4PFPuppi");
//if(debug)std::cout << "element " << myelement->evaluate({0.7, 258.})  << std::endl;
  //  if(debug)std::cout << "Reading Jet Energy Corrections from files: " << JECfile<<", "<< corrfile << endl;

    for(unsigned int i=0; i<pt.size(); i++){
       if(debug)std::cout  << "Raw pt is: " << pt[i] << "\t\t";
       if(debug)std::cout  << "Eta is: " << eta[i] << "\t\t";
       if(debug)std::cout  << "Area is: " << area[i] << "\t\t";
       if(debug)std::cout  << "Rho is: " << rho << "\t\t";
/*
if(debug)std::cout  <<std::endl;
       float sf = cset->at(corrfile + "L1FastJet_AK4PFPuppi")->evaluate({area[i], eta[i], pt[i], rho});
if(debug)std::cout  << "sf is: " << sf <<std::endl; 
       sf *= cset->at(corrfile + "L2Relative_AK4PFPuppi")->evaluate({eta[i], pt[i]});
       //       sf *= cset->at(corrfile + "L2Relative_AK4PFPuppi")->evaluate({eta[i], phi[i], pt[i]}); // For phi dependent jecs (2023Bpix, 2024)
if(debug)std::cout  << "sf is: " << sf <<std::endl; 
       sf *= cset->at(corrfile + "L3Absolute_AK4PFPuppi")->evaluate({eta[i], pt[i]});
if(debug)std::cout  << "sf is: " << sf <<std::endl; 
       sf *= cset->at(corrfile + "L2L3Residual_AK4PFPuppi")->evaluate({eta[i], pt[i]});
if(debug)std::cout  << "sf is: " << sf <<std::endl; 
if(debug)std::cout  << "will now push"<<std::endl; */
    if(debug)std::cout  << "Area: " << area[i] << "\t\t" << "Eta: " << eta[i] << "\t\t" << "Pt: " << pt[i] << "\t\t" << "Rho: " << rho << std::endl;
    float sf = cset_l1jecs->evaluate({area[i], eta[i], pt[i], rho});
    if(debug)std::cout  << "L1 SF: " << cset_l1jecs->evaluate({area[i], eta[i], pt[i], rho})<< std::endl;
    if(!isData){  
          sf *= cset_l2jecs->evaluate({eta[i], pt[i]});
          if(debug)std::cout  << "L2 SF: " << cset_l2jecs->evaluate({eta[i], pt[i]})<< std::endl;
    }
    sf *= cset_l3jecs->evaluate({eta[i], pt[i]});
    if(debug)std::cout  << "L3 SF: " << cset_l3jecs->evaluate({eta[i], pt[i]})<< std::endl;
    sf *= cset_residualjecs->evaluate({eta[i], pt[i]});
    if(debug)std::cout  << "Residuals SF: " << cset_residualjecs->evaluate({eta[i], pt[i]})<< std::endl; 
    if(debug)std::cout  << "Final SF: " << sf << std::endl;       
    Jet_corPt.push_back(sf*pt[i]);

    }

    return Jet_corPt;

}
    ''')



def JECsInit(year, era, isData):
    period  = ''
    
    if year == 2022:
        JECversion = 'V2'
        if era in ['C','D']:
            JECtag     = '2022_Summer22'
            JECname    = 'Summer22_22Sep2023'
            period     = 'CD'
        else:
            JECtag     = "2022_Summer22EE"
            JECname    = "Summer22EE_22Sep2023"
            period     = era
    elif year == 2023:
        JECversion = 'V1'
        if era in ['Cv123', 'Cv4', 'C']:
            JECtag     = '2023_Summer23'
            JECname    = 'Summer23Prompt23'
            if isData:
                period = era
        else:
            JECtag     = '2023_Summer23BPix'
            JECname    = 'Summer23BPixPrompt23'
            period = era

    elif year == 2024:
        JECversion = 'V2'
        JECtag     = '2024_Winter24'
        JECname    = 'Winter24Prompt24'
        if era in ['B', 'C', 'D', 'BCD']:
            period     = 'BCD'
        elif era in ['E', 'F', 'G', 'H', 'I']:
            period     = era

    JECfile = f'JEC/{JECtag}/jet_jerc.json.gz'

    if isData:
       corrfile = f'{JECname}_Run{period}_{JECversion}_DATA_'
    else: 
       corrfile = f'{JECname}_{JECversion}_MC_'

    print(JECfile)
    print(corrfile)
    return JECfile, corrfile
