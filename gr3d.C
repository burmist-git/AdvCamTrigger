#include <TPolyLine3D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLegend.h>

void gr3d ( Double_t zmin=-10., Double_t zmax=+10. )
{
   Int_t nPoints=5 ;
   Double_t x1[] = {1., 3., 5., 7., 9.} ;
   Double_t x2[] = {1., 3., 5., 7., 9.} ;
   Double_t y1[] = {0., 0., 0., 0., 0.} ; 
   Double_t y2[] = {3., 3., 3., 3., 3.} ;
   Double_t z1[] = {0., 3.,-2., 3., 4.} ; 
   Double_t z2[] = {3.,-1., 2., 2., 0.} ;

   TPolyLine3D *line3D_1 = new TPolyLine3D (nPoints,x1,y1,z1) ; 
   line3D_1->SetLineColor(kRed ) ;
   TPolyLine3D *line3D_2 = new TPolyLine3D (nPoints,x2,y2,z2) ; 
   line3D_2->SetLineColor(kBlue) ;

   Double_t xmin=0., xmax=10., resolution=0.5 ; 
   Int_t nBins = Int_t((xmax-xmin)/resolution) ;
   /// The binning should be adjusted to the interactive zoom resolution desired

   TH2F* histo = new TH2F("histo","", nBins,xmin,xmax, 1,y1[0],y2[0]) ;
   histo->SetStats(kFALSE) ;
   histo->SetMinimum(zmin) ; 
   histo->SetMaximum(zmax) ;
   histo->SetXTitle("#lambda  [nm]"    ) ; 
   histo->GetXaxis()->CenterTitle() ;
   histo->SetYTitle("time  [min]"      ) ; 
   histo->GetYaxis()->CenterTitle() ;
   histo->SetZTitle("intensity  [a.u.]") ; 
   histo->GetZaxis()->CenterTitle() ;

   TCanvas *canvas = new TCanvas("canvas","canvas",0,0,1000,618) ;
   canvas->SetTheta(23.) ; 
   canvas->SetPhi(-23.) ;
   gPad->SetLeftMargin(0.18) ; 
   gPad->SetRightMargin(0.10) ; 
   gPad->SetTopMargin(0.20) ; 
   gPad->SetBottomMargin(0.20) ;
   histo->Draw("lego0,fb") ; 
   line3D_1->Draw() ; 
   line3D_2->Draw() ;
   TLegend *leg = new TLegend(0.81,0.86,0.99,0.94) ; 
   leg->SetTextSize(0.04) ; 
   leg->SetBorderSize(0.) ;
   leg->AddEntry (line3D_1,"data 1","L") ; 
   leg->AddEntry (line3D_2,"data 2","L") ;
   leg->Draw() ;
   //canvas->SaveAs("test.png") ;

   return ;
}
