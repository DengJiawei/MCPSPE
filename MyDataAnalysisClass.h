
#ifndef MyDataAnalysisClass_h
#define MyDataAnalysisClass_h

#ifdef __MAKECINT__
#pragma link C++ class Signal + ;
#pragma link C++ class vector < Signal > +;
#pragma link C++ class MyDataAnalysisClass + ;
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TString.h"
#include <TStyle.h>

using namespace std;

struct Signal
{
    double baseline;
    double amplitude;
    double area;
    double starttime_CFT;
    double starttime_threshold;

    double width;
    double risetime;
    double falltime;
    double FWHM;
    bool theMaxSignal;
    Int_t startposition;
    Int_t endposition;
};

// 读取文件列表
void GetFileList(TString filePath, TString filePattern, vector<TString> &fList)
{
    char line[1000];
    fList.clear();

    FILE *fp = gSystem->OpenPipe("ls " + filePath + "/*" + filePattern, "r");
    if (!fp)
    {
        cout << "----> NO data files exists in " << filePath << "!" << endl;
        return;
    }

    while (fgets(line, sizeof(line), fp))
    {
        TString s(line);
        if (s.Index(filePattern) == -1)
            continue;
        fList.push_back(s.ReplaceAll("\n", ""));
    }
}

//开始处理数据前需要知道的文件信息有：
//  读取文件： 文件的列数、行数
//  处理数据： 每个点之间的时间间隔，第0位置：起始时间
//  单光子刻度： 方法1: 知道LED发光的时间位置（用数点来表示）

class MyDataAnalysisClass
{
private:
    // the information of raw data file
    Int_t FileQuantity;           //!
    Int_t FileMaxRow;             //!//the raw number is from 0 to max-1
    Int_t FileMaxColumn;          //! // t,1,2,3,4
    Double_t FileTimeUnitAverage; //!
    Double_t FileTimeUnit1;       //!
    Double_t FileTimeUnit2;       //!
    Double_t FileStartTime;

    Double_t *rawdata_channel_t = NULL; //!
    Double_t *rawdata_channel_1 = NULL; //!
    Double_t *rawdata_channel_2 = NULL; //!
    Double_t *rawdata_channel_3 = NULL; //!
    Double_t *rawdata_channel_4 = NULL; //!ll

public:
    MyDataAnalysisClass(TString FileName_fp);
    MyDataAnalysisClass(Int_t FileMaxRow_fp, Int_t FileMaxColumn_fp, Double_t FileTimeUnitAverage_fp, Double_t FileStartTime_fp);
    ~MyDataAnalysisClass();

    void GetFileInformation(TString FileName_fp);
    void CheckFileInformation();

    void ReadOneFile(TString FileName_fp); //输入文件名，读取文件中的数据
    void CheckResultsofOneFile();
    void CheckHistogram();

    //*******************************************************************
    //函数参数中要输出、赋值的（引用）放在前面，需要传入的放在后面

    void DrawHist_BybinWidth_GausFit(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, const Double_t binWidth_fp);
    void DrawHist_BybinWidth_GausFit(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, const Double_t binWidth_fp, bool check_fp, TString savename_fp);

    //用来大致确定信号范围的
    void FindAverageBaseline(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, TString savename_fp);
    void FindAverageBaseline(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp);

    //for SPE
    //measure informations
    vector<Int_t> UnitFlag; //单光子刻度时候分区间
    void SetUnitFlag(const Double_t TimeSignalStart, const Double_t SignalPeriod_fp, const Int_t SignalQuantity);
    void SetLEDFlag(const Double_t LEDWidth_fp, const Double_t TimeSignalStart_fp, const Double_t SignalPeriod_fp, const Int_t SignalQuantity_fp);
    vector<int> vLEDstart;
    vector<int> vLEDend;


    vector<double> vLEDBaseline;
    vector<double> vLEDBaselinesigma;
    vector<double> vLEDArea;
    void WorkOnAFile_SimpleSum();

    void CheckOverPulse(vector<int> vCheckSignal_start_fp, vector<int> vCheckSignal_end_fp, int vChecki_fp, vector<int> &vOverPulse_fp, vector<double> &vOverAmpl_fp, double baseline_fp);
    void WorkOnAFile_SPE();
    void FindTheAmplitude_min(const double *data_fp, int start_fp, int end_fp, double &A_fp);
    void FindtheMaxSignalInUnit(const Double_t *pdata_fp, const int start_fp, const int end_fp);

    void WorkOnAFile_FindMin();
    vector<int> vlowposi;
    vector<double> vlowest;

    //20191209
    void FindTheLowestPosition(const double * data, int start_fp, int end_fp, int & lowposi_fp);
};


void MyDataAnalysisClass::FindTheLowestPosition(const double * data_fp, int start_fp, int end_fp, int & lowposi_fp)
{
    double lowest = data_fp[start_fp];
    for(int i = start_fp; i < end_fp; i ++)
    {
        if(data_fp[i] < lowest)
        {
            lowest = data_fp[i];
            lowposi_fp = i;
        }
    }
    
}


//0k
MyDataAnalysisClass::MyDataAnalysisClass(TString FileName_fp)
{
    GetFileInformation(FileName_fp);
    CheckFileInformation();

    switch (FileMaxColumn)
    {
    case 4:
        cout << " creating 4 " << endl;
        rawdata_channel_4 = new Double_t[FileMaxRow]();

    case 3:
        cout << " creating 3 " << endl;
        rawdata_channel_3 = new Double_t[FileMaxRow]();

    case 2:
        cout << " creating 2 " << endl;
        rawdata_channel_2 = new Double_t[FileMaxRow]();

    case 1:
        cout << " creating 1 " << endl;
        rawdata_channel_1 = new Double_t[FileMaxRow]();

    case 0:
        cout << " creating t " << endl;
        rawdata_channel_t = new Double_t[FileMaxRow]();
        break;

    default:
        cout << " wrong FileMaxColumn " << endl;
        exit(EXIT_FAILURE);
        break;
    }
}

MyDataAnalysisClass::MyDataAnalysisClass(Int_t FileMaxRow_fp, Int_t FileMaxColumn_fp, Double_t FileTimeUnitAverage_fp, Double_t FileStartTime_fp)
{
    FileMaxRow = FileMaxRow_fp;
    FileMaxColumn = FileMaxColumn_fp;
    FileTimeUnitAverage = FileTimeUnitAverage_fp;
    FileStartTime = FileStartTime_fp;

    switch (FileMaxColumn)
    {
    case 4:
        cout << " creating 4 " << endl;
        rawdata_channel_4 = new Double_t[FileMaxRow]();

    case 3:
        cout << " creating 3 " << endl;
        rawdata_channel_3 = new Double_t[FileMaxRow]();

    case 2:
        cout << " creating 2 " << endl;
        rawdata_channel_2 = new Double_t[FileMaxRow]();

    case 1:
        cout << " creating 1 " << endl;
        rawdata_channel_1 = new Double_t[FileMaxRow]();

    case 0:
        cout << " creating t " << endl;
        rawdata_channel_t = new Double_t[FileMaxRow]();
        break;

    default:
        cout << " wrong FileMaxColumn " << endl;
        exit(EXIT_FAILURE);
        break;
    }
};

//ok
MyDataAnalysisClass::~MyDataAnalysisClass()
{
    switch (FileMaxColumn)
    {
    case 4:
        cout << " deleting 4 " << endl;
        delete[] rawdata_channel_4;
    case 3:
        cout << " deleting 3 " << endl;
        delete[] rawdata_channel_3;

    case 2:
        cout << " deleting 2 " << endl;
        delete[] rawdata_channel_2;

    case 1:
        cout << " deleting 1 " << endl;
        delete[] rawdata_channel_1;

    case 0:
        cout << " deleting 0 " << endl;
        delete[] rawdata_channel_t;
        break;

    default:
        cout << " wrong FileMaxColumn " << endl;
        exit(EXIT_FAILURE);
        break;
    }
}

//ok
//打开一个文件，通过数行中的逗号‘，’确定列数；
//再通过依次读取行数得到行数；
//读取前三行的第一列数据（时间）得到每个数据点的间隔，将个数与时间联系起来。
void MyDataAnalysisClass::GetFileInformation(TString FileName_fp)
{
    ifstream checkStream;
    checkStream.open(FileName_fp);
    if (checkStream.is_open())
    {
        cout << "open the check file successfully : " << FileName_fp << endl;
    }
    else
    {
        cout << "something wrong when get the file information " << endl;
        exit(EXIT_FAILURE);
    }

    TString firstline;
    firstline.ReadLine(checkStream);
    FileMaxRow = 1;
    FileMaxColumn = firstline.CountChar(',');

    cout << "by counting char ',', get the FileMaxColumn " << FileMaxColumn << endl;

    // Double_t time_temp[3];
    // TString rest_temp;
    // for(int i = 0; i < 3; i ++)
    // {
    // 	checkStream >> time_temp[i] >> rest_temp;
    // }
    Double_t time_temp[3];
    //	TString rest_temp;
    char rest_temp[100];
    for (int i = 0; i < 3; i++)
    {
        checkStream >> time_temp[i];
        checkStream.getline(rest_temp, 40);
        cout << time_temp[i] << rest_temp << endl;
    }

    FileMaxRow += 3;
    FileStartTime = time_temp[0];

    FileTimeUnit1 = time_temp[1] - time_temp[0];
    FileTimeUnit2 = time_temp[2] - time_temp[1];
    FileTimeUnitAverage = (time_temp[2] - time_temp[0]) / 2;

    while (firstline.ReadLine(checkStream))
    {
        FileMaxRow++;
    }
    cout << "the raw quantity in a file is " << FileMaxRow << endl;
    checkStream.close();
}

//ok
void MyDataAnalysisClass::CheckFileInformation()
{
    cout << endl
         << endl;
    cout << " Make sure the file informations following are right " << endl;
    cout << " ******************************************************************* " << endl;
    cout << " the file max row is " << FileMaxRow << endl;
    cout << " the file max column is " << FileMaxColumn << endl;
    cout << " the file time units are " << FileTimeUnit1 << "  and  " << FileTimeUnit2 << endl;
    cout << " the file time unit average is " << FileTimeUnitAverage << endl;
    cout << " ******************************************************************* " << endl;
    cout << endl
         << endl;
    cout << " input \"yes\" to continue ,input \"no\" to modify and others to quit " << endl;

    TString mks;
    if (cin >> mks)
    {
        if (mks == "yes")
        {
            cout << "finish check the file information " << endl;
        }
        else if (mks == "no")
        {
            cout << "the file information is wrong,please check it again " << endl;
            exit(EXIT_FAILURE);
        }
        else
        {
            cout << " you choose to quit the program " << endl;
            exit(EXIT_FAILURE);
        }
    }
}

//ok
//assign file data to array
void MyDataAnalysisClass::ReadOneFile(TString FileName_fp)
{
    ifstream readStream;
    readStream.open(FileName_fp);
    if (readStream.is_open())
    {
        cout << "open the file named " << FileName_fp << endl;
    }
    else
    {
        cout << "cannot open the file " << FileName_fp << "please check it again " << endl;
        exit(EXIT_FAILURE); //system("pause")
    }

    char douhao;
    switch (FileMaxColumn)
    {

    case 4:
        for (int k = 0; k < FileMaxRow; k++)
        {
            readStream >> rawdata_channel_t[k] >> douhao >> rawdata_channel_1[k] >> douhao >> rawdata_channel_2[k] >> douhao >> rawdata_channel_3[k] >> douhao >> rawdata_channel_4[k];
        }
        break;

    case 3:
        for (int k = 0; k < FileMaxRow; k++)
        {
            readStream >> rawdata_channel_t[k] >> douhao >> rawdata_channel_1[k] >> douhao >> rawdata_channel_2[k] >> douhao >> rawdata_channel_3[k];
        }
        break;

    case 2:
        for (int k = 0; k < FileMaxRow; k++)
        {
            readStream >> rawdata_channel_t[k] >> douhao >> rawdata_channel_1[k] >> douhao >> rawdata_channel_2[k];
        }
        break;

    case 1:
        for (int k = 0; k < FileMaxRow; k++)
        {
            readStream >> rawdata_channel_t[k] >> douhao >> rawdata_channel_1[k];
        }
        break;

    case 0:
        for (int k = 0; k < FileMaxRow; k++)
        {
            readStream >> rawdata_channel_t[k];
        }
        break;

    default:
        cout << " wrong FileMaxRow or wrong file " << endl;
        exit(EXIT_FAILURE);
        break;
    }

    readStream.close();
}



void MyDataAnalysisClass::FindAverageBaseline(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp)
{
    // cout << " find average baseline :  dataquantity = " << dataquantity_fp << endl;
    Double_t dataMax = 0;
    Double_t dataMin = 0;

    for (int i = 0; i < dataquantity_fp; i++)
    {
        if (data_fp[i] > dataMax)
        {
            dataMax = data_fp[i];
        }
        if (data_fp[i] < dataMin)
        {
            dataMin = data_fp[i];
        }
    }
    //    if(dataMin >= 0)
    //    {
    //        cout << "the dataMin is positive ,please check the data " << endl;
    //        exit(EXIT_FAILURE);
    //    }
    //    if(dataMax <= 0)
    //    {
    //        cout << "the dataMax is nagetive ,please check the data " << endl;
    //        exit(EXIT_FAILURE);
    //    }

    // cout << "max  = " << dataMax << endl;
    // cout << "min  = " << dataMin << endl;

    dataMin = -0.1;
    dataMax = 0.1;

    double binWidth_fp = 0.0025;
    Int_t BinNum = int(2.0 * (dataMax - dataMin) / binWidth_fp + 0.5);

    //    auto hist_temp = new TH1D("hist_temp", "hist_temp", BinNum, 2 * dataMin, 2 * dataMax);
    //    auto hist_temp = new TH1D("hist_temp", "hist_temp", 40, -0.1, 0.1);
    TH1D hist_temp("hist_temp", "hist_temp", 80, -0.1, 0.1);
    for (int i = 0; i < dataquantity_fp; i++)
    {
        if (data_fp[i] > -0.1 && data_fp[i] < 0.1)
            hist_temp.Fill(data_fp[i]);
    }

    TF1 gausf("gausf", "gaus", 2 * dataMin, 2 * dataMax);

    gausf.SetParameter(0, hist_temp.GetMaximum());
    gausf.SetParameter(1, hist_temp.GetMean());
    gausf.SetParLimits(1, -0.2, 0.2);
    gausf.SetParameter(2, hist_temp.GetStdDev());
    gausf.SetParLimits(2, -0.2, 0.2);

    //注意加上“N”不画图像
    hist_temp.Fit("gausf", "QMN");

    mean_fp = gausf.GetParameter(1);
    sigma_fp = gausf.GetParameter(2);

    // cout << "the baseline is " << mean_fp << " +- " << sigma_fp << endl;
}

void MyDataAnalysisClass::FindAverageBaseline(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, TString savename_fp)
{
    cout << " find average baseline " << endl;
    Double_t dataMax = 0;
    Double_t dataMin = 0;
    cout << " the dataquantity_fp = " << dataquantity_fp << endl;
    for (int i = 0; i < dataquantity_fp; i++)
    {
        if (data_fp[i] > dataMax)
        {
            dataMax = data_fp[i];
        }
        if (data_fp[i] < dataMin)
        {
            dataMin = data_fp[i];
        }
    }
    //    if(dataMin >= 0)
    //    {
    //        cout << "the dataMin is positive ,please check the data " << endl;
    //        exit(EXIT_FAILURE);
    //    }
    //    if(dataMax <= 0)
    //    {
    //        cout << "the dataMax is nagetive ,please check the data " << endl;
    //        exit(EXIT_FAILURE);
    //    }

    cout << "max  = " << dataMax << endl;
    cout << "min  = " << dataMin << endl;
    double binWidth_fp = 0.005;
    Int_t BinNum = int(0.5 + 2 * (dataMax - dataMin) / binWidth_fp);

    auto tc_temp = new TCanvas();
    //    auto hist_temp = new TH1D("hist_temp", "hist_temp", BinNum, 2 * dataMin, 2 * dataMax);
    auto hist_temp = new TH1D("hist_temp", "hist_temp", 120, -0.1, 0.1);
    for (int i = 0; i < dataquantity_fp; i++)
    {
        if (data_fp[i] > -0.1 && data_fp[i] < 0.1)
            hist_temp->Fill(data_fp[i]);
    }
    auto gausf = new TF1("gausf", "gaus", 2 * dataMin, 2 * dataMax);

    gausf->SetParameter(0, hist_temp->GetMaximum());
    gausf->SetParameter(1, hist_temp->GetMean());
    gausf->SetParLimits(1, -0.2, 0.2);
    gausf->SetParameter(2, hist_temp->GetStdDev());
    gausf->SetParLimits(2, -0.2, 0.2);

    hist_temp->Fit("gausf", "QM");

    mean_fp = gausf->GetParameter(1);
    sigma_fp = gausf->GetParameter(2);
    cout << "the baseline is " << mean_fp << " +- " << sigma_fp << endl;
    //hist_temp->Draw();
    tc_temp->Draw();
    tc_temp->SaveAs(savename_fp);
    delete hist_temp;
    delete gausf;
    delete tc_temp;
}

//ok
void MyDataAnalysisClass::DrawHist_BybinWidth_GausFit(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, const Double_t binWidth_fp)
{
    //    cout << " draw hist by bin width " << endl;
    Double_t dataMax = 0;
    Double_t dataMin = 0;
    //
    //    for (int i = 0; i < dataquantity_fp; i++)
    //    {
    //        if (data_fp[i] > dataMax)
    //        {
    //            dataMax = data_fp[i];
    //        }
    //        if (data_fp[i] < dataMin)
    //        {
    //            dataMin = data_fp[i];
    //        }
    //    }
    //    cout << "max  = " << dataMax << endl;
    //    cout << "min  = " << dataMin << endl;
    dataMin = -0.1;
    dataMax = 0.1;
    Int_t binNum_temp = int(0.5 + 2 * (dataMax - dataMin) / binWidth_fp);

    //    auto tc_temp = new TCanvas();
    auto hist_temp = new TH1D("hist_temp", "hist_temp", binNum_temp, 2 * dataMin, 2 * dataMax);
    for (int i = 0; i < dataquantity_fp; i++)
    {
        hist_temp->Fill(data_fp[i]);
    }
    auto gausf = new TF1("gausf", "gaus", 2 * dataMin, 2 * dataMax);
    gausf->SetParameter(0, hist_temp->GetMaximum());
    gausf->SetParLimits(0, 1, hist_temp->GetEntries());
    gausf->SetParameter(1, hist_temp->GetMean());
    gausf->SetParLimits(1, -0.2, 0.2);
    gausf->SetParameter(2, hist_temp->GetStdDev());
    gausf->SetParLimits(2, -0.2, 0.2);
    hist_temp->Fit("gausf", "QMN");

    mean_fp = gausf->GetParameter(1);
    sigma_fp = gausf->GetParameter(2);
    cout << "the baseline is " << mean_fp << " +- " << sigma_fp << endl;

    delete hist_temp;
    delete gausf;
    //    delete tc_temp;
}

void MyDataAnalysisClass::DrawHist_BybinWidth_GausFit(Double_t &mean_fp, Double_t &sigma_fp, const Double_t *data_fp, const Int_t dataquantity_fp, const Double_t binWidth_fp, bool check_fp, TString savename_fp)
{
    cout << " draw hist by bin width " << endl;
    Double_t dataMax = 0;
    Double_t dataMin = 0;

    //    for (int i = 0; i < dataquantity_fp; i++)
    //    {
    //        if (data_fp[i] > dataMax)
    //        {
    //            dataMax = data_fp[i];
    //        }
    //        if (data_fp[i] < dataMin)
    //        {
    //            dataMin = data_fp[i];
    //        }
    //    }
    //
    //    cout << "max  = " << dataMax << endl;
    //    cout << "min  = " << dataMin << endl;
    dataMin = -0.1;
    dataMax = 0.1;
    Int_t binNum_temp = int(0.5 + 2 * (dataMax - dataMin) / binWidth_fp);

    auto tc_temp = new TCanvas();
    auto hist_temp = new TH1D("hist_temp", "hist_temp", binNum_temp, 2 * dataMin, 2 * dataMax);
    for (int i = 0; i < dataquantity_fp; i++)
    {
        hist_temp->Fill(data_fp[i]);
    }
    auto gausf = new TF1("gausf", "gaus", 2 * dataMin, 2 * dataMax);
    gausf->SetParameter(0, hist_temp->GetMaximum());
    gausf->SetParLimits(0, 1, hist_temp->GetEntries());
    gausf->SetParameter(1, hist_temp->GetMean());
    gausf->SetParLimits(1, -0.2, 0.2);
    gausf->SetParameter(2, hist_temp->GetStdDev());
    gausf->SetParLimits(2, -0.2, 0.2);
    hist_temp->Fit("gausf", "QM");

    mean_fp = gausf->GetParameter(1);
    sigma_fp = gausf->GetParameter(2);
    cout << "the baseline is " << mean_fp << " +- " << sigma_fp << endl;
    if (check_fp)
    {
        tc_temp->SaveAs(savename_fp);
    }
    delete hist_temp;
    delete gausf;
    delete tc_temp;
}

void MyDataAnalysisClass::FindTheAmplitude_min(const double *data_fp, int start_fp, int end_fp, double &A_fp)
{
    double min = data_fp[start_fp];
    for (int i = start_fp; i < end_fp; i++)
    {
        // cout << "data_fp[" << i << "] = " << data_fp[i] << endl;

        if (min > data_fp[i])
        {

            min = data_fp[i];
            // cout << "                        data_fp[" << i <<"] = " << data_fp[i] << ";  min = " << min << endl;
        }
    }
    A_fp = min;
}

// void MyDataAnalysisClass::CheckResultsofOneFile()
// {
//     // for(unsigned i = 0; i < (UnitFlag.size()); i ++)
//     // {
//     //     cout << "UnitFlag.at(" << i << ") = " << UnitFlag.at(i) << endl;
//     // }
//     cout << "UnitFlag.at(" << 0 << ") = " << UnitFlag.at(0) << endl;
//     cout << "UnitFlag.at(" << UnitFlag.size() - 1 << ") = " << UnitFlag.at(UnitFlag.size() - 1) << endl;

//     auto tc_check = new TCanvas();
//     auto tg = new TGraph(FileMaxRow, rawdata_channel_t, rawdata_channel_1);
//     tg->Draw();

//     for (unsigned i = 0; i < vBaselineInOneFile.size(); i++)
//     {
//         auto bs = new TF1("bs", "[0]", rawdata_channel_t[UnitFlag.at(i)], rawdata_channel_t[UnitFlag.at(i + 1)]);
//         bs->SetParameter(0, vBaselineInOneFile.at(i));
//         bs->SetLineColor(kRed);
//         bs->Draw("same");

//         auto bs2 = new TF1("bs2", "[0]", rawdata_channel_t[UnitFlag.at(i)], rawdata_channel_t[UnitFlag.at(i + 1)]);
//         bs2->SetParameter(0, vUnitThresholdInOneFile.at(i));
//         bs2->SetLineColor(kRed);
//         bs2->Draw("same");
//     }

//     cout << "vSignalsInSingleFile.size() = " << vSignalsInSingleFile.size() << endl;

//     for (unsigned i = 0; i < vSignalsInSingleFile.size(); i++)
//     {
//         cout << "line start " << rawdata_channel_t[vSignalsInSingleFile.at(i).startposition] << endl;
//         cout << "line end " << rawdata_channel_t[vSignalsInSingleFile.at(i).endposition] << endl;
//         auto tf = new TF1("tf", "[0]", rawdata_channel_t[vSignalsInSingleFile.at(i).startposition], rawdata_channel_t[vSignalsInSingleFile.at(i).endposition]);
//         //        auto tf = new TF1("tf","[0]",rawdata_channel_t[0],rawdata_channel_t[FileMaxRow-1]);
//         tf->SetParameter(0, vSignalsInSingleFile.at(i).baseline);
//         tf->SetLineColor(kBlue);
//         tf->Draw("same");

//         auto tf2 = new TF1("tf2", "[0]", rawdata_channel_t[vSignalsInSingleFile.at(i).startposition], rawdata_channel_t[vSignalsInSingleFile.at(i).endposition]);
//         int basei_temp = (vSignalsInSingleFile.at(i).startposition - UnitFlag.at(0)) / (UnitFlag.at(1) - UnitFlag.at(0));

//         tf2->SetParameter(0, vUnitThreshold2InOneFile.at(basei_temp));
//         tf2->SetLineColor(kGreen);
//         tf2->Draw("same");

//         auto tf3 = new TF1("tf3", "[0]", rawdata_channel_t[vSignalsInSingleFile.at(i).startposition], rawdata_channel_t[vSignalsInSingleFile.at(i).endposition]);
//         tf3->SetParameter(0, 0.1 * vSignalsInSingleFile.at(i).amplitude + vSignalsInSingleFile.at(i).baseline);
//         tf3->SetLineColor(kBlack);
//         tf3->Draw("same");

//         cout << "vSignalsInSingleFile.at(" << i << ").baseline = " << vSignalsInSingleFile.at(i).baseline << endl;
//         cout << endl
//              << endl;
//     }

//     tc_check->Draw();
// }

void MyDataAnalysisClass::SetLEDFlag(const Double_t LEDWidth_fp, const Double_t TimeSignalStart_fp, const Double_t SignalPeriod_fp, const Int_t SignalQuantity_fp)
{
    vLEDstart.clear();
    vLEDend.clear();
    int ledwidth_temp = int(LEDWidth_fp / FileTimeUnitAverage + 0.5);
    int startpo_temp = int((TimeSignalStart_fp - FileStartTime) / FileTimeUnitAverage + 0.5);
    int distant_temp = int(SignalPeriod_fp / FileTimeUnitAverage + 0.5);
    for (int i = 0; i < SignalQuantity_fp; i++)
    {
        vLEDstart.push_back(startpo_temp + i * distant_temp);
        vLEDend.push_back(startpo_temp + i * distant_temp + ledwidth_temp);
    }
};

void MyDataAnalysisClass::WorkOnAFile_SimpleSum()
{
    vLEDBaseline.clear();
    vLEDBaselinesigma.clear();
    vLEDArea.clear();

    const double *time_temp = this->rawdata_channel_t;
    const double *data_temp = this->rawdata_channel_1;
    const int move_1ns = int(1e-9 / FileTimeUnitAverage + 0.5);

    for (int i = 0; i < vLEDstart.size(); i++)
    {
        double mean;
        double sigma;
        DrawHist_BybinWidth_GausFit(mean, sigma, data_temp + vLEDstart.at(i) - 20 * move_1ns, 20 * move_1ns, 0.0025);
        vLEDBaseline.push_back(mean);
        vLEDBaselinesigma.push_back(sigma);
    }

    for (int i = 0; i < vLEDstart.size(); i++)
    {
        double area_temp = 0;
        for (int j = vLEDstart.at(i); j < vLEDend.at(i); j++)
        {
            area_temp += (data_temp[j + 1] + data_temp[j]) / 2.0 * (time_temp[j + 1] - time_temp[j]);
        }
        area_temp -= vLEDBaseline.at(i) * (time_temp[vLEDend.at(i)] - time_temp[vLEDstart.at(i)]);
        vLEDArea.push_back(area_temp);
    }
}

void MyDataAnalysisClass::WorkOnAFile_FindMin()
{
    vLEDBaseline.clear();
    vLEDBaselinesigma.clear();
    vLEDArea.clear();

    vlowest.clear();
    vlowposi.clear();

    const double *time_temp = this->rawdata_channel_t;
    const double *data_temp = this->rawdata_channel_1;
    const int move_1ns = int(1e-9 / FileTimeUnitAverage + 0.5);

    for (int i = 0; i < vLEDstart.size(); i++)
    {
        double mean;
        double sigma;
        DrawHist_BybinWidth_GausFit(mean, sigma, data_temp + vLEDstart.at(i) - 20 * move_1ns, 20 * move_1ns, 0.0025);
        vLEDBaseline.push_back(mean);
        vLEDBaselinesigma.push_back(sigma);
        int l = 0;
        FindTheLowestPosition(data_temp,vLEDstart.at(i),vLEDend.at(i), l);
        vlowposi.push_back(l);
        vlowest.push_back(data_temp[l]);

    }



    for (int i = 0; i < vLEDstart.size(); i++)
    {
        double area_temp = 0;
        for (int j = vlowposi.at(i)-int(1.5*move_1ns); j < vlowposi.at(i)+int(2*move_1ns); j++)
        {
            area_temp += (data_temp[j + 1] + data_temp[j]) / 2.0 * (time_temp[j + 1] - time_temp[j]);
        }
        area_temp -= vLEDBaseline.at(i) * (time_temp[vlowposi.at(i)-int(1.5*move_1ns)] - time_temp[vlowposi.at(i)+int(2*move_1ns)]);
        vLEDArea.push_back(area_temp);
    }
}


#endif
