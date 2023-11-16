Double_t get_angle_between_optical_axis_particle( Double_t tel_theta, Double_t tel_phi, Double_t part_theta, Double_t part_phi);

Int_t angle_between_optical_axis_particle(){
  //
  Double_t tel_theta = 20.0/180.0*TMath::Pi();
  Double_t tel_phi   = 180.0/180.0*TMath::Pi();
  //
  Double_t part_theta = 20.0/180.0*TMath::Pi();
  Double_t part_phi   = 175.0/180.0*TMath::Pi();
  Double_t part_theta_max = tel_theta + 5.0/180.0*TMath::Pi();
  Double_t part_theta_min = tel_theta - 5.0/180.0*TMath::Pi();
  Double_t part_phi_max   = tel_phi   + 10.0/180.0*TMath::Pi();
  Double_t part_phi_min   = tel_phi   - 10.0/180.0*TMath::Pi();
  //
  Int_t n_bin_phi   = 1000;
  Int_t n_bin_theta = 1000;
  Double_t bin_width_phi = (part_phi_max - part_phi_min)/n_bin_phi;
  Double_t bin_width_theta = (part_theta_max - part_theta_min)/n_bin_theta;
  //
  Double_t val;
  //
  TH2D *h2 = new TH2D("h2","h2",
		      n_bin_phi, part_phi_min*180.0/TMath::Pi(), part_phi_max*180.0/TMath::Pi(),
		      n_bin_theta, 90-part_theta_max*180.0/TMath::Pi(), 90-part_theta_min*180.0/TMath::Pi());
  for(Int_t i = 1;i<=n_bin_phi;i++){
    part_theta = part_theta_min + bin_width_theta/2.0 + bin_width_theta*(i-1);
    for(Int_t j = 1;j<=n_bin_theta;j++){
      part_phi = part_phi_min + bin_width_phi/2.0 + bin_width_phi*(j-1);
      //val = (i-1)*n_bin_theta+j;
      val = get_angle_between_optical_axis_particle( tel_theta, tel_phi, part_theta, part_phi);
      //val = part_phi*180.0/TMath::Pi();
      //val = part_theta*180.0/TMath::Pi();
      h2->SetBinContent(j,i,val);
    }
  }
  //
  //cout<<get_angle_between_optical_axis_particle( tel_theta, tel_phi, part_theta, part_phi)<<endl;
  //
  //h2->Draw("TEXT");
  //
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  //gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  //
  h2->SetTitle("Angle between optical axis and particle");
  h2->SetMinimum(0);
  h2->SetMaximum(3);
  //
  h2->Draw("ZCOLOR");
 
  h2->GetXaxis()->SetTitle("Azimuth, deg");
  h2->GetYaxis()->SetTitle("Altitude, deg");

  //  
  return 0;
}

Double_t get_angle_between_optical_axis_particle( Double_t tel_theta, Double_t tel_phi, Double_t part_theta, Double_t part_phi){
  TVector3 tel;
  TVector3 part;  
  tel.SetMagThetaPhi(1.0,tel_theta,tel_phi);
  part.SetMagThetaPhi(1.0,part_theta,part_phi);
  return TMath::ACos(tel.Dot(part))*180.0/TMath::Pi();
}
