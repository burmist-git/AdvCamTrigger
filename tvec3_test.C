Int_t tvec3_test(){
  //
  Double_t tel_theta = 20.0/180.0*TMath::Pi();
  Double_t tel_phi   = 180.0/180.0*TMath::Pi();
  //
  Double_t part_theta = 20.0/180.0*TMath::Pi();
  Double_t part_phi   = 180.0/180.0*TMath::Pi();
  //
  TVector3 vx_tel(1.0,0.0,0.0);
  TVector3 vy_tel(0.0,1.0,0.0);
  TVector3 vz_tel(0.0,0.0,1.0);
  //////////////////////////
  vx_tel.RotateZ(tel_phi);
  vx_tel.RotateY(-tel_theta);
  cout<<"vx_tel.x() "<<vx_tel.x()<<endl
      <<"vx_tel.y() "<<vx_tel.y()<<endl
      <<"vx_tel.z() "<<vx_tel.z()<<endl;    
  //
  vy_tel.RotateZ(tel_phi);
  vy_tel.RotateY(-tel_theta);
  cout<<"vy_tel.x() "<<vy_tel.x()<<endl
      <<"vy_tel.y() "<<vy_tel.y()<<endl
      <<"vy_tel.z() "<<vy_tel.z()<<endl;
  //
  vz_tel.RotateZ(tel_phi);
  vz_tel.RotateY(-tel_theta);
  cout<<"vz_tel.x() "<<vz_tel.x()<<endl
      <<"vz_tel.y() "<<vz_tel.y()<<endl
      <<"vz_tel.z() "<<vz_tel.z()<<endl;    
  //
  TVector3 part;
  part.SetMagThetaPhi(1.0,part_theta,part_phi);
  cout<<"part.x() "<<part.x()<<endl
      <<"part.y() "<<part.y()<<endl
      <<"part.z() "<<part.z()<<endl;    
  //
  TVector3 v_in_tel( vx_tel.Dot(part), vy_tel.Dot(part), vz_tel.Dot(part));
  //
  //
  //
  TVector3 vr(1.0,1.0,1.0);
  cout<<"vr.x() "<<vr.x()<<endl
      <<"vr.y() "<<vr.y()<<endl
      <<"vr.z() "<<vr.z()<<endl;      
  //
  vr.RotateY(TMath::Pi());
  cout<<"vr.x() "<<vr.x()<<endl
      <<"vr.y() "<<vr.y()<<endl
      <<"vr.z() "<<vr.z()<<endl;
  //
  vr.RotateZ(TMath::Pi());
  cout<<"vr.x() "<<vr.x()<<endl
      <<"vr.y() "<<vr.y()<<endl
      <<"vr.z() "<<vr.z()<<endl;
  //
  //
  //
  /*
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
  // 
  h2->GetXaxis()->SetTitle("Azimuth, deg");
  h2->GetYaxis()->SetTitle("Altitude, deg");
  //
  tel_theta = 20.0/180.0*TMath::Pi();
  tel_phi   = 180.0/180.0*TMath::Pi();
  //
  part_theta = 20.0/180.0*TMath::Pi();
  part_phi   = 177.5/180.0*TMath::Pi();
  cout<<get_angle_between_optical_axis_particle( tel_theta, tel_phi, part_theta, part_phi)/180.0*TMath::Pi()<<endl;
  //
  part_theta = 21.0/180.0*TMath::Pi();
  part_phi   = 180.0/180.0*TMath::Pi();
  cout<<get_angle_between_optical_axis_particle( tel_theta, tel_phi, part_theta, part_phi)/180.0*TMath::Pi()<<endl;
  //
  part_theta = 21.0/180.0*TMath::Pi();
  part_phi   = 177.5/180.0*TMath::Pi();
  cout<<get_angle_between_optical_axis_particle( tel_theta, tel_phi, part_theta, part_phi)/180.0*TMath::Pi()<<endl;  
  //
  //
  Double_t part_theta_tel;
  Double_t part_phi_tel;
  //
  //
  tel_theta = 10.0/180.0*TMath::Pi();
  tel_phi   = 180.0/180.0*TMath::Pi();
  //
  //
  part_theta = 20.0/180.0*TMath::Pi();
  part_phi   = 180.0/180.0*TMath::Pi();
  //
  //
  Double_t tel_alpha;
  Double_t tel_betta;
  Double_t tel_gamma;
  get_Euler_alpha_betta_gamma( tel_theta, tel_phi, tel_alpha, tel_betta, tel_gamma);
  cout<<"tel_alpha "<<tel_alpha*180.0/TMath::Pi()<<endl
      <<"tel_betta "<<tel_betta*180.0/TMath::Pi()<<endl
      <<"tel_gamma "<<tel_gamma*180.0/TMath::Pi()<<endl;  
  //
  //
  get_theta_phi_in_telescope_frame(tel_alpha, tel_betta, tel_gamma, part_theta, part_phi, part_theta_tel, part_phi_tel);
  cout<<"part_theta_tel "<<part_theta_tel*180.0/TMath::Pi()<<endl
      <<"part_phi_tel   "<<part_phi_tel*180.0/TMath::Pi()<<endl;
  //
  //
  Double_t theta_in_tel;
  Double_t phi_in_tel;
  part_theta = 90.0/180.0*TMath::Pi();
  part_phi   = 180.0/180.0*TMath::Pi();
  get_part_coordinates_in_tel_frame( part_theta, part_phi, theta_in_tel, phi_in_tel);
  //
  cout<<"theta_in_tel "<<theta_in_tel*180.0/TMath::Pi()<<endl
      <<"phi_in_tel   "<<phi_in_tel*180.0/TMath::Pi()<<endl;
  //
  part_theta = 90.0/180.0*TMath::Pi();
  part_phi   = 0.0/180.0*TMath::Pi();
  test( part_theta, part_phi);
  //

  */

  return 0;
}

/*
Double_t get_angle_between_optical_axis_particle( Double_t tel_theta, Double_t tel_phi, Double_t part_theta, Double_t part_phi){
  TVector3 tel;
  TVector3 part;  
  tel.SetMagThetaPhi(1.0,tel_theta,tel_phi);
  part.SetMagThetaPhi(1.0,part_theta,part_phi);
  return TMath::ACos(tel.Dot(part))*180.0/TMath::Pi();
}

//void get_theta_phi_in_telescope_frame( Double_t tel_theta, Double_t tel_phi, Double_t part_theta, Double_t part_phi, Double_t &part_theta_tel, Double_t &part_phi_tel){
void get_theta_phi_in_telescope_frame( Double_t tel_alpha, Double_t tel_betta, Double_t tel_gamma,
				       Double_t part_theta, Double_t part_phi, Double_t &part_theta_tel, Double_t &part_phi_tel){
  //
  part_theta_tel = 0.0;
  part_phi_tel = 0.0;
  TVector3 part;
  part.SetMagThetaPhi(1.0,part_theta,part_phi);
  part.RotateZ(tel_alpha);
  part.RotateX(tel_betta);
  part.RotateZ(tel_gamma);
  part_theta_tel = part.Theta();
  part_phi_tel = part.Phi();
  //return TMath::ACos(tel.Dot(part))*180.0/TMath::Pi();
}

Bool_t get_Euler_alpha_betta_gamma( Double_t theta, Double_t phi, Double_t &alpha, Double_t &betta, Double_t &gamma){
  TVector3 vx(1.0,0.0,0.0);
  TVector3 vy(0.0,1.0,0.0);
  TVector3 vz(0.0,0.0,1.0);
  //
  vx.RotateZ(phi);
  vx.RotateY(theta);
  //cout<<"vx.x() "<<vx.x()<<endl
  //   <<"vx.y() "<<vx.y()<<endl
  //   <<"vx.z() "<<vx.z()<<endl;    
  //
  vy.RotateZ(phi);
  vy.RotateY(theta);
  //cout<<"vy.x() "<<vy.x()<<endl
  //  <<"vy.y() "<<vy.y()<<endl
  //  <<"vy.z() "<<vy.z()<<endl;
  //
  vz.RotateZ(phi);
  vz.RotateY(theta);
  //cout<<"vz.x() "<<vz.x()<<endl
  //  <<"vz.y() "<<vz.y()<<endl
  //  <<"vz.z() "<<vz.z()<<endl;      
  //
  Double_t Z3 = vz.z();
  Double_t Z2 = vz.y();
  Double_t Y3 = vy.z();
  //
  if(Z3*Z3 != 1.0){
    alpha = TMath::ACos(-Z2/TMath::Sqrt(1.0-Z3*Z3));
    betta = TMath::ACos(Z3);
    gamma = TMath::ACos(Y3/TMath::Sqrt(1.0-Z3*Z3));
    return true;
  }
  alpha = phi/2.0;
  betta = 0.0;
  gamma = phi/2.0;
  return false;
}

void get_part_coordinates_in_tel_frame( Double_t theta, Double_t phi, Double_t &theta_in_tel, Double_t &phi_in_tel){
  //////////////////////////
  Double_t tel_theta = 20.0/180.0*TMath::Pi();
  Double_t tel_phi   = 180.0/180.0*TMath::Pi();
  //////////////////////////
  TVector3 vx_tel(1.0,0.0,0.0);
  TVector3 vy_tel(0.0,1.0,0.0);
  TVector3 vz_tel(0.0,0.0,1.0);
  //////////////////////////
  vx_tel.RotateZ(tel_phi);
  vx_tel.RotateY(tel_theta);
  cout<<"vx_tel.x() "<<vx_tel.x()<<endl
      <<"vx_tel.y() "<<vx_tel.y()<<endl
      <<"vx_tel.z() "<<vx_tel.z()<<endl;    
  //
  vy_tel.RotateZ(tel_phi);
  vy_tel.RotateY(tel_theta);
  vz_tel.RotateZ(tel_phi);
  vz_tel.RotateY(tel_theta);
  //////////////////////////
  TVector3 part;
  part.SetMagThetaPhi(1.0,theta,phi);
  //
  TVector3 v_in_tel(vx_tel.Dot(part),vy_tel.Dot(part),vz_tel.Dot(part));
  theta_in_tel = v_in_tel.Theta();
  phi_in_tel = v_in_tel.Phi();
}

void test( Double_t theta, Double_t phi){
    TVector3 part;
    part.SetMagThetaPhi(1.0,theta,phi);
    cout<<"part.x() "<<part.x()<<endl
	<<"part.y() "<<part.y()<<endl
	<<"part.z() "<<part.z()<<endl;
}
*/
