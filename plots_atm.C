void read_atm(TString atmFileName, TGraph *gr_Alt_km, TGraph *gr_rho_g_cm3, TGraph *gr_thick_g_cm2, TGraph *gr_n_m, TGraph *gr_T_K, TGraph *gr_p_mbar, TGraph *gr_pw_p);  

Int_t plots_atm(){
  //
  TString atmFileName = "./atmprof_ecmwf_north_winter_fixed.dat";
  //
  TGraph *gr_Alt_km = new TGraph();
  TGraph *gr_rho_g_cm3 = new TGraph();
  TGraph *gr_thick_g_cm2 = new TGraph();
  TGraph *gr_n_m = new TGraph();
  TGraph *gr_T_K = new TGraph();
  TGraph *gr_p_mbar = new TGraph();
  TGraph *gr_pw_p = new TGraph();
  //
  read_atm( atmFileName, gr_Alt_km, gr_rho_g_cm3, gr_thick_g_cm2, gr_n_m, gr_T_K, gr_p_mbar, gr_pw_p);  
  //
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  TCanvas *c1 = new TCanvas("c1","gr_Alt_km",10,10,600,600);
  gr_Alt_km->SetLineColor(kBlack);
  gr_Alt_km->SetLineWidth(3.0);
  gr_Alt_km->SetMarkerColor(kBlack);
  gr_Alt_km->SetMarkerStyle(20);
  gr_Alt_km->Draw("APL");
  //
  TCanvas *c2 = new TCanvas("c2","gr_rho_g_cm3",10,10,600,600);
  gr_rho_g_cm3->SetLineColor(kBlack);
  gr_rho_g_cm3->SetLineWidth(3.0);
  gr_rho_g_cm3->SetMarkerColor(kBlack);
  gr_rho_g_cm3->SetMarkerStyle(20);
  gr_rho_g_cm3->Draw("APL");
  //
  TCanvas *c3 = new TCanvas("c3","gr_thick_g_cm2",10,10,600,600);
  gr_thick_g_cm2->SetLineColor(kBlack);
  gr_thick_g_cm2->SetLineWidth(3.0);
  gr_thick_g_cm2->SetMarkerColor(kBlack);
  gr_thick_g_cm2->SetMarkerStyle(20);
  gr_thick_g_cm2->Draw("APL");
  //
  TCanvas *c4 = new TCanvas("c4","gr_n_m",10,10,600,600);
  gr_n_m->SetLineColor(kBlack);
  gr_n_m->SetLineWidth(3.0);
  gr_n_m->SetMarkerColor(kBlack);
  gr_n_m->SetMarkerStyle(20);
  gr_n_m->Draw("APL");
  //
  TCanvas *c5 = new TCanvas("c5","gr_T_K",10,10,600,600);
  gr_T_K->SetLineColor(kBlack);
  gr_T_K->SetLineWidth(3.0);
  gr_T_K->SetMarkerColor(kBlack);
  gr_T_K->SetMarkerStyle(20);
  gr_T_K->Draw("APL");
  //
  TCanvas *c6 = new TCanvas("c6","gr_p_mbar",10,10,600,600);
  gr_p_mbar->SetLineColor(kBlack);
  gr_p_mbar->SetLineWidth(3.0);
  gr_p_mbar->SetMarkerColor(kBlack);
  gr_p_mbar->SetMarkerStyle(20);
  gr_p_mbar->Draw("APL");
  //
  TCanvas *c7 = new TCanvas("c7","gr_pw_p",10,10,600,600);
  gr_pw_p->SetLineColor(kBlack);
  gr_pw_p->SetLineWidth(3.0);
  gr_pw_p->SetMarkerColor(kBlack);
  gr_pw_p->SetMarkerStyle(20);
  gr_pw_p->Draw("APL");
  //
  return 0;
}

/////////////////////////////////////////////////////////////////////////////////////////
//# ----------------------------------------------------------------------------------------
//#
//# New profile generated with the following command:
//#    ./check_atmprof cfg/CTA/atmprof_ecmwf_north_winter_nexp.dat = 2.180 28.762166097 51 atmprof_ecmwf_north_winter_fixed.dat
//#
//# Primary profile: cfg/CTA/atmprof_ecmwf_north_winter_nexp.dat
//# Secondary profile: cfg/CTA/atmprof_ecmwf_north_winter_nexp.dat
//# Observation level (not part of profile): 2.180 km a.s.l.
//# Geographic latitude: 28.762166 deg (for gravity).
//# Switch-over from primary to secondary above 51.000 km a.s.l.
//#
//# ----------------------------------------------------------------------------------------
//#
//#Col. #1          #2           #3            #4        [ #5 ]     [ #6 ]       [ # 7 ]
//# Alt [km]    rho [g/cm^3] thick [g/cm^2]    n-1        T [K]    p [mbar]      pw / p
//   0.000     1.21171e-03  1.04473e+03   2.80503e-04   291.91   1.02049e+03  1.45399e-02
//   1.000     1.10668e-03  9.28757e+02   2.55945e-04   284.50   9.06892e+02  1.03318e-02
//   2.000     9.94136e-04  8.23744e+02   2.29611e-04   281.43   8.04127e+02  4.46223e-03
//    .
//    .
//    .
// 100.000     5.16989e-10  2.75865e-04   1.19328e-10   185.38   2.75100e-04  0.00000e+00
// 105.000     2.03675e-10  1.07349e-04   4.70107e-11   197.19   1.15283e-04  0.00000e+00
// 110.000     8.14116e-11  4.09604e-05   1.87908e-11   224.14   5.23782e-05  0.00000e+00
// 115.000     3.45898e-11  1.36099e-05   7.98376e-12   268.51   2.66596e-05  0.00000e+00
// 120.000     2.09799e-11  0.00000e+00   4.84242e-12   333.43   0.00000e+00  0.00000e+00
/////////////////////////////////////////////////////////////////////////////////////////

void read_atm(TString atmFileName, TGraph *gr_Alt_km, TGraph *gr_rho_g_cm3, TGraph *gr_thick_g_cm2, TGraph *gr_n_m, TGraph *gr_T_K, TGraph *gr_p_mbar, TGraph *gr_pw_p){
  //
  ifstream fileIn(atmFileName.Data());
  cout<<"atm file : "<<atmFileName<<std::endl;
  //
  Double_t Alt_km, rho_g_cm3, thick_g_cm2, n_m, T_K, p_mbar, pw_p;
  TString mot;
  Int_t counter = 0;
  //
  if (fileIn.is_open()){
    //
    while(fileIn>>mot)
      if(mot == "[g/cm^2]")
	break;
    //n-1        T [K]    p [mbar]      pw / p
    fileIn>>mot; // n-1
    fileIn>>mot; // T
    fileIn>>mot; // [K]
    fileIn>>mot; // p
    fileIn>>mot; // [mbar]
    fileIn>>mot; // pw
    fileIn>>mot; // /
    fileIn>>mot; // p
    cout<<"mot = "<<mot<<endl;
    //
    while(fileIn>>Alt_km>>rho_g_cm3>>thick_g_cm2>>n_m>>T_K>>p_mbar>>pw_p){
      gr_Alt_km->SetPoint(counter,counter,Alt_km);
      gr_rho_g_cm3->SetPoint(counter,Alt_km,rho_g_cm3);
      gr_thick_g_cm2->SetPoint(counter,Alt_km,thick_g_cm2);
      gr_n_m->SetPoint(counter,Alt_km,n_m);
      gr_T_K->SetPoint(counter,Alt_km,T_K);
      gr_p_mbar->SetPoint(counter,Alt_km,p_mbar);
      gr_pw_p->SetPoint(counter,Alt_km,pw_p);
      counter++;
    }
    fileIn.close();
    //
  }
  else {
    std::cout<<"Unable to open file"<<std::endl;
    assert(0);
  }
  //
  //for(Int_t i = 0;i<grtmp->GetN();i++){
  //grtmp->GetPoint(i,x,y);
  //_gr_wf_tmpl->SetPoint(i,x,grtmpstrached->Eval(x));
  //}
}
