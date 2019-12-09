#include "TTree.h"
#include "TFile.h"
#include "MyDataAnalysisClass.h"
#include "TStyle.h"
#include "TROOT.h"


#include "TGFrame.h"
#include "TGFileDialog.h"
#include "TGCanvas.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGMsgBox.h"
#include "TGSlider.h"
#include "TGTab.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGTextEdit.h"
#include "TGComboBox.h"
#include "TG3DLine.h"
#include "TGClient.h"
#include "TGResourcePool.h"
#include "TStyle.h"
#include "TStyleManager.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "Riostream.h"

#include "TSystem.cxx"




//#include "TSystem.h"
//#include "/usr/local/Cellar/root/6.18.04/include/root/TSystem.h"
TString FilePath = "/Users/dengjiawei/code_project/LSMCP/20191204";
// "HV3003V_2MHz1400mV22ns_190Cyc10msPeri_00000.csv"FileCondition +
TString FileCondition = FilePath(FilePath.Last('/')+1,FilePath.Sizeof()-FilePath.Last('/')-2);

void MyMain_SPE()
{

    vector<TString> FileList;
    FileList.clear();
    GetFileList(FilePath, ".csv", FileList);

    MyDataAnalysisClass mydataanalysisclass(FileList.at(0));
    // mydataanalysisclass.CheckFileInformation();
    mydataanalysisclass.SetLEDFlag(25e-9,0.065e-6,500e-9,190);

    gStyle->SetOptStat(111111);
    gStyle->SetOptFit(1111);

    TString rootname = FileCondition + TString(".root");
    //TString rootname =  TString("test.root");

    TFile myfile(rootname, "recreate");
    TTree tree("tree", "tree");

    tree.Branch("signal", &mydataanalysisclass.vLEDArea);
    tree.Branch("Min",&mydataanalysisclass.vlowest);

    for (auto filename : FileList)
    // TString filename = FileList.at(0);
    {
        mydataanalysisclass.ReadOneFile(filename);
        mydataanalysisclass.WorkOnAFile_SimpleSum();
        tree.Fill();
    }
    tree.Write();
    myfile.Close();
}
