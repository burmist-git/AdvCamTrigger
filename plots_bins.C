void copyHist(TH1D *h1, TH1D *h1_to_cp);

Int_t plots_bins(){

  TH1D *h1 = new TH1D();
  h1->SetNameTitle("h1","h1");
  Double_t gg[5]={0,1,2,3,4};
  h1->SetBins(4,gg);

  h1->Fill(0.5);
  h1->Fill(0.5);
  h1->Fill(1.5);
  h1->Fill(0.5);

  h1->Fill(3.5);

  TH1D *h12 = new TH1D();  
  copyHist(h12, h1);
  h12->Draw();
  
  return 0;
}

void copyHist(TH1D *h1, TH1D *h1_to_cp){
  int nBins = h1_to_cp->GetNbinsX();
  double *bins_low_edge= new double[nBins+1];
  for(int i = 1;i<=nBins;i++)
    bins_low_edge[i-1] = h1_to_cp->GetBinLowEdge(i);
  bins_low_edge[nBins] = h1_to_cp->GetBinLowEdge(nBins) + h1_to_cp->GetBinWidth(nBins);
  h1->SetBins(nBins,bins_low_edge);
  //  
  for(Int_t i = 1;i<=h1_to_cp->GetNbinsX();i++)
    h1->SetBinContent(i,h1_to_cp->GetBinContent(i));
}
