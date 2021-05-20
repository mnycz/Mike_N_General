#include "DataQuality.h"

using namespace ROOT;

// structure to return two values

struct Histo_Values{
  double hist_mean;
  double hist_sigma;
};
// Function for Gaussian fits
Histo_Values mean_value(TF1* fit,TH1D histo){
  histo.Fit(fit,"QR");
  Histo_Values tmp1;
  tmp1.hist_mean = fit->GetParameter("Mean");
  tmp1.hist_sigma = fit->GetParameter("Sigma");
  return tmp1; 
}


void DataQuality(){

  string Spectrometer;
  cout<<"Do you want to check HMS or SHMS?";
  cin>>Spectrometer;
  cout<<"Data quality check will be run for:"<<Spectrometer<<endl;

  TString cerenkov;
  TString event_type;
  TString spec_flag;
  TString rootfiles;
  if (Spectrometer =="HMS"){
    spec_flag ="H";
    event_type="2";
    rootfiles="HMS";
    cerenkov = "cer";
   cout<<"spec_flag"<<"\t"<<spec_flag<<endl;
  }
  else{
    cerenkov = "ngcer";
    spec_flag ="P";
    event_type="1";
    rootfiles="SHMS";
    cout <<"P"<<endl;
  }
  // cout<<"root file"<<rootfiles<<endl;
  const char *dirname;

  ofstream outFile;
   // Location of rootfiles 
  if (!Spectrometer.compare("HMS")){
    dirname = "/lustre19/expphy/volatile/hallc/c-polhe3/mnycz/ROOTfiles/HMS/";
    outFile.open("HMS_Data_Summary.csv");
  }
  else if (!Spectrometer.compare("SHMS")){
    dirname = "/lustre19/expphy/volatile/hallc/c-polhe3/mnycz/ROOTfiles/SHMS/";
    outFile.open("SHMS_Data_Summary.csv");
   }
  else{
    cout<<"Which spectrometer do you want to check????"<<endl;
  }

  const char *ext = ".root";
  vector<string> Filename;
  vector<int> Run_Number;

// We want to check all the HMS and SHMS runs (separately). This below will only check fully replayed rootfiles (_-1.root) 
  TSystemDirectory dir(dirname, dirname); 
  TList *files = dir.GetListOfFiles(); 
  if (files){
    TSystemFile *file;
    TString fname;
    TIter next(files);
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(ext)){
	Filename.push_back(fname.Data());
	cout <<fname.Data()<<endl;
      }
    }
  }

  // Get the Run# --- This assumes the file structure follows: #-1.root
  //TPRegexp filePattern("^(.*)_(\\d+)_(\\d+).root$");
  TPRegexp filePattern("(\\d+)");
    for (int i=0;i<Filename.size();i++){
      TString tmp =  Filename[i];
      int tmp_number = std::stoi(tmp(filePattern));
      Run_Number.push_back(tmp_number);
    }

  ROOT::EnableImplicitMT(1);
  TF1 *Ep_1 = new TF1("Ep_1","gaus",0.5,1.3);  // Fit for E/p
  TF1 *beta_1 = new TF1("beta_1","gaus",0.2,2);
  TF1 *NPE_1 = new TF1("NPE_1", "gaus",1,15);


  for (int j=0;j<Filename.size();j++){
 
    TString etot_cut = "";
    TString etrack_cut = "";
    TString ntrack_cut = "";
    TString combined_cut = "";
    combined_cut += Form("%s.cal.etot>0.2", spec_flag.Data());
    combined_cut += Form(" && %s.dc.ntrack>0", spec_flag.Data());
    etot_cut= Form("%s.cal.etot>0.2", spec_flag.Data());
    etrack_cut = Form("%s.cal.etracknorm>0.7",spec_flag.Data());
    ntrack_cut = Form("%s.dc.ntrack>0",spec_flag.Data());


    RDataFrame d1("T",dirname+Filename[j]);
    RDataFrame dTSP(Form("TS%s", spec_flag.Data()),dirname+Filename[j]);
    

    auto h = d1.Filter(Form(etrack_cut,spec_flag.Data())).Histo1D(
								  {"hist","E/p",100,0,2},Form("%s.cal.etracknorm",spec_flag.Data())
								  );
    auto h_beta = d1.Filter(Form(etrack_cut,spec_flag.Data())).Histo1D(
								  {"beta_hist","SHMS_beta",100,0,2},Form("%s.hod.beta",spec_flag.Data())
								  );  
    auto dc_tracks = d1.Filter(Form(ntrack_cut ,spec_flag.Data())).Histo1D(
								   {"tracks","dc_tracks",50,0,50},Form("%s.dc.ntrack",spec_flag.Data())
								   );
    auto NPE_Mean = d1.Filter(Form(combined_cut,spec_flag.Data(),spec_flag.Data())).Histo1D(
								   {"hist","NPE",40,0,2},Form("%s.%s.npeSum",spec_flag.Data(),cerenkov.Data())
							     	   );


    Histo_Values E_p1;
    Histo_Values E_p2;
    E_p1 = mean_value(Ep_1,*h);
 
    TF1 *Ep_2 = new TF1("Ep_2", "gaus",E_p1.hist_mean-E_p1.hist_sigma,E_p1.hist_mean+E_p1.hist_sigma);

    E_p2 = mean_value(Ep_2,*h);
    mean.push_back(E_p2.hist_mean);
    sigma.push_back(E_p2.hist_sigma);
 
    Histo_Values test;
    test =  mean_value(NPE_1,*NPE_Mean);
    npeSum.push_back(test.hist_mean);
    npeSigma.push_back(test.hist_sigma);

  // Beta fit                                                                                                                                                                                                         
    Histo_Values beta_p1;
    Histo_Values beta_p2;
    beta_p1 = mean_value(beta_1,*h_beta);
    TF1 *beta_2 = new TF1("beta_2", "gaus",beta_p1.hist_mean-beta_p1.hist_sigma,beta_p1.hist_mean+beta_p1.hist_sigma);

    beta_p2 = mean_value(beta_2,*h_beta);
    mean_beta.push_back(beta_p2.hist_mean);
    sigma_beta.push_back(beta_p2.hist_sigma);

  //************************************************//
  // Number of track in the drift chamber
  //auto dc_tracks = d1.Filter(Form("%s.dc.ntrack>0",spec_flag.Data())).Histo1D({"tracks","dc_tracks",50,0,50},Form("%s.dc.ntrack",spec_flag.Data()));
  track_mult.push_back(dc_tracks -> GetMean());
  track_mult_sigma.push_back(dc_tracks -> GetRMS());


  // Average drift time for each plane
  // Filter the histograms to have the same range as the plots : (300,-75,225)
  //*************************************************************************************// 
  auto h_u11= d1.Mean(Form("%s.dc.1u1.time",spec_flag.Data()));
  auto h_u12= d1.Mean(Form("%s.dc.1u2.time",spec_flag.Data()));
  auto h_v11= d1.Mean(Form("%s.dc.1v1.time",spec_flag.Data()));
  auto h_v12= d1.Mean(Form("%s.dc.1v2.time",spec_flag.Data()));
  auto h_x11= d1.Mean(Form("%s.dc.1x1.time",spec_flag.Data()));
  auto h_x12= d1.Mean(Form("%s.dc.1x2.time",spec_flag.Data()));

  auto h_u11_RMS= d1.StdDev(Form("%s.dc.1u1.time",spec_flag.Data()));
  auto h_u12_RMS= d1.StdDev(Form("%s.dc.1u2.time",spec_flag.Data()));
  auto h_v11_RMS= d1.StdDev(Form("%s.dc.1v1.time",spec_flag.Data()));
  auto h_v12_RMS= d1.StdDev(Form("%s.dc.1v2.time",spec_flag.Data()));
  auto h_x11_RMS= d1.StdDev(Form("%s.dc.1x1.time",spec_flag.Data()));
  auto h_x12_RMS= d1.StdDev(Form("%s.dc.1x2.time",spec_flag.Data()));
  //*************************************************************************************// 
  auto h_u21= d1.Mean(Form("%s.dc.1u1.time",spec_flag.Data()));
  auto h_u22= d1.Mean(Form("%s.dc.1u2.time",spec_flag.Data()));
  auto h_v21= d1.Mean(Form("%s.dc.1v1.time",spec_flag.Data()));
  auto h_v22= d1.Mean(Form("%s.dc.1v2.time",spec_flag.Data()));
  auto h_x21= d1.Mean(Form("%s.dc.1x1.time",spec_flag.Data()));
  auto h_x22= d1.Mean(Form("%s.dc.1x2.time",spec_flag.Data()));

  auto h_u21_RMS= d1.StdDev(Form("%s.dc.1u1.time",spec_flag.Data()));
  auto h_u22_RMS= d1.StdDev(Form("%s.dc.1u2.time",spec_flag.Data()));
  auto h_v21_RMS= d1.StdDev(Form("%s.dc.1v1.time",spec_flag.Data()));
  auto h_v22_RMS= d1.StdDev(Form("%s.dc.1v2.time",spec_flag.Data()));
  auto h_x21_RMS= d1.StdDev(Form("%s.dc.1x1.time",spec_flag.Data()));
  auto h_x22_RMS= d1.StdDev(Form("%s.dc.1x2.time",spec_flag.Data()));

  //*************************************************************************************// 
  u11_time_mean.push_back(*h_u11);
  u11_time_rms.push_back(*h_u11_RMS);  
  u12_time_mean.push_back(*h_u12);
  u12_time_rms.push_back(*h_u12_RMS);
  v11_time_mean.push_back(*h_v11);
  v11_time_rms.push_back(*h_v11_RMS);
  v12_time_mean.push_back(*h_v12);
  v12_time_rms.push_back(*h_v12_RMS);
  x11_time_mean.push_back(*h_x11);
  x11_time_rms.push_back(*h_x11_RMS);
  x12_time_mean.push_back(*h_x12);
  x12_time_rms.push_back(*h_x12_RMS);

  u21_time_mean.push_back(*h_u21);
  u21_time_rms.push_back(*h_u21_RMS);
  u22_time_mean.push_back(*h_u22);
  u22_time_rms.push_back(*h_u22_RMS);
  v21_time_mean.push_back(*h_v21);
  v21_time_rms.push_back(*h_v21_RMS);
  v22_time_mean.push_back(*h_v22);
  v22_time_rms.push_back(*h_v22_RMS);
  x21_time_mean.push_back(*h_x21);
  x21_time_rms.push_back(*h_x21_RMS);
  x22_time_mean.push_back(*h_x22);
  x22_time_rms.push_back(*h_x22_RMS);
  //*************************************************************************************// 


  auto Trigger = d1.Filter(Form("fEvtHdr.fEvtType ==%s",event_type.Data()));
  auto Pos_Hel = Trigger.Filter("T.helicity.hel > 0");
  auto Neg_Hel = Trigger.Filter("T.helicity.hel < 0");
  auto No_Hel = Trigger.Filter("T.helicity.hel == 0");

  auto nEntries = Trigger.Count();
  auto nPos_Hel = Pos_Hel.Count();
  auto nNeg_Hel = Neg_Hel.Count();
  auto nNo_Hel =  No_Hel.Count();
  auto Charge = dTSP.Max(Form("%s.BCM4B.scalerChargeCut",spec_flag.Data()));
  const double Total_Charge = *Charge / 1000.0;


  auto BCM_Current = dTSP.Filter(Form("%s.BCM4A.scalerCurrent>0.8",spec_flag.Data())).Histo1D(Form("%s.BCM4A.scalerCurrent",spec_flag.Data()));
  BCM_Current_Mean.push_back( BCM_Current->GetMean());


  auto BCM_Column = dTSP.Define("Current_Sq", Form("%s.BCM4A.scalerCurrent*%s.BCM4A.scalerCurrent",spec_flag.Data(),spec_flag.Data()));
  vector<double> I_sumsq2;
  BCM_Column.Foreach([&](double j){
      I_sumsq2.push_back(j);
    },{"Current_Sq"} );

  
  //Get the Median current value
  Double_t x, q;
  q=0.5;
  BCM_Current-> GetQuantiles(1,&x, &q);
  BCM_Current_Median.push_back(x);

 
  vector<double> Time_Weight;
  dTSP.Foreach([&](double i) {
      double tmp = 0.;
      tmp=i;
      Time_Weight.push_back(tmp);
    },{Form("%s.1MHz.scalerTime",spec_flag.Data())} );
  //cout<<"TW"<<"\t"<<Time_Weight[0]<<endl;

  dTSP.Foreach([&](double j){
      I_sum.push_back(j);
    },{Form("%s.BCM4A.scalerCurrent",spec_flag.Data())} );

  vector<double> Temp_Time;
  for(int i=0;i<Time_Weight.size()-1;i++){
    Temp_Time.push_back(abs(Time_Weight[i] - Time_Weight[i+1]));
  }

  for (int k;k<(Time_Weight.size())-1;k++){
    Weighted_Current += I_sum[k]*Temp_Time[k];
    Weighted_Current_Sq += I_sumsq2[k]*(Temp_Time[k]*Temp_Time[k]);
  }


  auto I_Sumsq = BCM_Column.Sum("Current_Sq");
  auto I_Sum = dTSP.Sum(Form("%s.BCM4A.scalerCurrent",spec_flag.Data()));
  auto I_entries = dTSP.Count();
  Duty_Factor.push_back((Weighted_Current*Weighted_Current)/(*I_entries*Weighted_Current_Sq) );
  // cout<<"d"<<"\t"<<Duty_Factor[0]<<endl;


  //*******Trigger Rates*******//
  std::map<std::string, std::string> triggers {
    {"TRIG1_uA", Form("%s.hTRIG1.scalerRate",spec_flag.Data())},
    {"TRIG2_uA", Form("%s.hTRIG2.scalerRate",spec_flag.Data())},
    {"TRIG3_uA", Form("%s.hTRIG3.scalerRate",spec_flag.Data())},
    {"TRIG4_uA", Form("%s.hTRIG4.scalerRate",spec_flag.Data())},
    {"TRIG5_uA", Form("%s.hTRIG5.scalerRate",spec_flag.Data())},
    {"TRIG6_uA", Form("%s.hTRIG6.scalerRate",spec_flag.Data())},
    };

  vector<double> tmp_tg;
  for (auto &it : triggers) {
    auto htrg = dTSP.Histo1D(it.second);
    auto tmp_y = htrg->GetMaximum();
    auto tmp_x = htrg->GetBinContent(tmp_y);
    TF1 *rate = new TF1("rate", "gaus",tmp_x-500,tmp_x+500);
    htrg ->Fit(rate,"QR");
    tmp_tg.push_back(rate->GetParameter("Mean"));
  }
  T1.push_back(tmp_tg[0]);
  T2.push_back(tmp_tg[1]);
  T3.push_back(tmp_tg[2]);
  T4.push_back(tmp_tg[3]);
  T5.push_back(tmp_tg[4]);
  T6.push_back(tmp_tg[5]);

  //**************************//

  }

  Int_t N_Runs = mean.size();
  Double_t RunNumber[mean.size()];
  Double_t Mean_Ep[mean.size()];
  Double_t Sigma_Ep[mean.size()];
  Double_t Current_Mean_Median[mean.size()];
  Double_t Mean_Beta[mean.size()];
  Double_t Sigma_Beta[mean.size()];
  for (int i=0;i<mean.size();i++){
    RunNumber[i]=Run_Number[i];
    Mean_Ep[i] = mean[i];
    Sigma_Ep[i] = sigma[i];
    Current_Mean_Median[i] = abs(BCM_Current_Mean[i] -  BCM_Current_Median[i]);
    Mean_Beta[i] = mean_beta[i];
    Sigma_Beta[i] = sigma_beta[i];

  }

outFile<<"Run_Number"<<","<<"E_p"<<","<<"E_p_sigma"<<","
         <<"u11_time_mean"<<","<<"u11_time_rms"<<","
	 <<"v11_time_mean"<<","<<"v11_time_rms"<<","
         <<"x11_time_mean"<<","<<"x11_time_rms"<<","
         <<"u12_time_mean"<<","<<"u12_time_rms"<<","
         <<"v12_time_mean"<<","<<"v12_time_rms"<<","
         <<"x12_time_mean"<<","<<"x12_time_rms"<<","
	 <<"u21_time_mean"<<","<<"u21_time_rms"<<","
         <<"v21_time_mean"<<","<<"v21_time_rms"<<","
         <<"x21_time_mean"<<","<<"x21_time_rms"<<","
         <<"u22_time_mean"<<","<<"u22_time_rms"<<","
         <<"v22_time_mean"<<","<<"v22_time_rms"<<","
         <<"x22_time_mean"<<","<<"x22_time_rms"<<","
	 <<"Beta_Mean"<<","<<"Beta_Sigma"<<","
	 <<"Current_Mean"<<","<<"Current_Median"<<","
	 <<"NPE_Sum"<<","<<"NPE_Sigma"<<","
	 <<"Mean_Num_Track"<<","<<"RMS_Num_Track"<<","
	 <<"T1_Rate"<<","<<"T2_Rate"<<","<<"T3_Rate"<<","
	 <<"T4_Rate"<<","<<"T5_Rate"<<","<<"T6_Rate"<<","
	 <<"Duty_Factor"<<endl;

  for (int i=0;i<mean.size();i++){
    outFile<<Run_Number[i]<<","<<Mean_Ep[i]<<","<<Sigma_Ep[i]<<","
	   <<u11_time_mean[i]<<","<<u11_time_rms[i]<<","
	   <<v11_time_mean[i]<<","<<x11_time_rms[i]<<","
	   <<v11_time_mean[i]<<","<<x11_time_rms[i]<<","
           <<u12_time_mean[i]<<","<<u12_time_rms[i]<<","
           <<v12_time_mean[i]<<","<<x12_time_rms[i]<<","
           <<v12_time_mean[i]<<","<<x12_time_rms[i]<<","      
	   <<u21_time_mean[i]<<","<<u21_time_rms[i]<<","
           <<v21_time_mean[i]<<","<<x21_time_rms[i]<<","
           <<v21_time_mean[i]<<","<<x21_time_rms[i]<<","
           <<u22_time_mean[i]<<","<<u22_time_rms[i]<<","
           <<v22_time_mean[i]<<","<<x22_time_rms[i]<<","
           <<v22_time_mean[i]<<","<<x22_time_rms[i]<<","
	   <<mean_beta[i]<<","<<sigma_beta[i]<<","
           <<BCM_Current_Mean[i]<<","<<BCM_Current_Median[i]<<","
	   <<npeSum[i]<<","<<npeSigma[i]<<","
	   <<track_mult[i]<<","<< track_mult_sigma[i]<<","
	   <<T1[i]/BCM_Current_Mean[i] <<","
           <<T2[i]/BCM_Current_Mean[i] <<","
           <<T3[i]/BCM_Current_Mean[i] <<","
	   <<T4[i]/BCM_Current_Mean[i] <<","
           <<T5[i]/BCM_Current_Mean[i] <<","
	   <<T6[i]/BCM_Current_Mean[i] <<","
	   <<Duty_Factor[i]<<endl;
  }


  outFile.close();
  return 0;
}

