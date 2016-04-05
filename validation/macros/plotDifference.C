#include "RiceStyle.h"

using namespace std;

void plotDifference(){
	
	TFile* file1 = new TFile("../files/800_pre1.root");
	TFile* file2 = new TFile("../files/800_pre6.root");

	TH1D* Ntrk1 = (TH1D*) file1->Get("ana/Ntrk");
	TH1D* nVtx1 = (TH1D*) file1->Get("ana/nVtx");
	TH1D* vtxTracksSize1 = (TH1D*) file1->Get("ana/vtxTracksSize");
	TH1D* vtxZ1 = (TH1D*) file1->Get("ana/vtxZ");
	TH1D* vtxX1 = (TH1D*) file1->Get("ana/vtxX");
	TH1D* vtxY1 = (TH1D*) file1->Get("ana/vtxY");
	TH1D* pt1 = (TH1D*) file1->Get("ana/pt");
	TH1D* ptError1 = (TH1D*) file1->Get("ana/ptError");
	TH1D* DCAz1 = (TH1D*) file1->Get("ana/DCAz");
	TH1D* DCAxy1 = (TH1D*) file1->Get("ana/DCAxy");
	TH1D* eta1 = (TH1D*) file1->Get("ana/eta");
	TH1D* phi1 = (TH1D*) file1->Get("ana/phi");
	TH1D* Chi2n1 = (TH1D*) file1->Get("ana/Chi2n");
	TH1D* numberOfHits1 = (TH1D*) file1->Get("ana/numberOfHits");
	TH1D* Algo1 = (TH1D*) file1->Get("ana/Algo");
	TH2D* caloVsCbin1 = (TH2D*) file1->Get("ana/caloVsCbin");

	Ntrk1->SetMarkerStyle(20);
	nVtx1->SetMarkerStyle(20);
	vtxTracksSize1->SetMarkerStyle(20);
	vtxZ1->SetMarkerStyle(20);
	vtxX1->SetMarkerStyle(20);
	vtxY1->SetMarkerStyle(20);
	pt1->SetMarkerStyle(20);
	ptError1->SetMarkerStyle(20);
	DCAz1->SetMarkerStyle(20);
	DCAxy1->SetMarkerStyle(20);
	eta1->SetMarkerStyle(20);
	phi1->SetMarkerStyle(20);
	Chi2n1->SetMarkerStyle(20);
    numberOfHits1->SetMarkerStyle(20);
	Algo1->SetMarkerStyle(20);

	Ntrk1->SetMarkerColor(kBlack);
	nVtx1->SetMarkerColor(kBlack);
	vtxTracksSize1->SetMarkerColor(kBlack);
	vtxZ1->SetMarkerColor(kBlack);
	vtxX1->SetMarkerColor(kBlack);
	vtxY1->SetMarkerColor(kBlack);
	pt1->SetMarkerColor(kBlack);
	ptError1->SetMarkerColor(kBlack);
	DCAz1->SetMarkerColor(kBlack);
	DCAxy1->SetMarkerColor(kBlack);
	eta1->SetMarkerColor(kBlack);
	phi1->SetMarkerColor(kBlack);
	Chi2n1->SetMarkerColor(kBlack);
    numberOfHits1->SetMarkerColor(kBlack);
	Algo1->SetMarkerColor(kBlack);

	Ntrk1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	nVtx1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	vtxTracksSize1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	vtxZ1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	vtxX1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	vtxY1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	pt1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	ptError1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	DCAz1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	DCAxy1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	eta1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	phi1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	Chi2n1->Scale( 1.0/ (vtxZ1->GetEntries()) );
    numberOfHits1->Scale( 1.0/ (vtxZ1->GetEntries()) );
	Algo1->Scale( 1.0/ (vtxZ1->GetEntries()) );


	TH1D* Ntrk2 = (TH1D*) file2->Get("ana/Ntrk");
	TH1D* nVtx2 = (TH1D*) file2->Get("ana/nVtx");
	TH1D* vtxTracksSize2 = (TH1D*) file2->Get("ana/vtxTracksSize");
	TH1D* vtxZ2 = (TH1D*) file2->Get("ana/vtxZ");
	TH1D* vtxX2 = (TH1D*) file2->Get("ana/vtxX");
	TH1D* vtxY2 = (TH1D*) file2->Get("ana/vtxY");
	TH1D* pt2 = (TH1D*) file2->Get("ana/pt");
	TH1D* ptError2 = (TH1D*) file2->Get("ana/ptError");
	TH1D* DCAz2 = (TH1D*) file2->Get("ana/DCAz");
	TH1D* DCAxy2 = (TH1D*) file2->Get("ana/DCAxy");
	TH1D* eta2 = (TH1D*) file2->Get("ana/eta");
	TH1D* phi2 = (TH1D*) file2->Get("ana/phi");
	TH1D* Chi2n2 = (TH1D*) file2->Get("ana/Chi2n");
	TH1D* numberOfHits2 = (TH1D*) file2->Get("ana/numberOfHits");
	TH1D* Algo2 = (TH1D*) file2->Get("ana/Algo");
	TH2D* caloVsCbin2 = (TH2D*) file2->Get("ana/caloVsCbin");

	Ntrk2->SetMarkerStyle(20);
	nVtx2->SetMarkerStyle(20);
	vtxTracksSize2->SetMarkerStyle(20);
	vtxZ2->SetMarkerStyle(20);
	vtxX2->SetMarkerStyle(20);
	vtxY2->SetMarkerStyle(20);
	pt2->SetMarkerStyle(20);
	ptError2->SetMarkerStyle(20);
	DCAz2->SetMarkerStyle(20);
	DCAxy2->SetMarkerStyle(20);
	eta2->SetMarkerStyle(20);
	phi2->SetMarkerStyle(20);
	Chi2n2->SetMarkerStyle(20);
    numberOfHits2->SetMarkerStyle(20);
	Algo2->SetMarkerStyle(20);

	Ntrk2->SetMarkerColor(kRed);
	nVtx2->SetMarkerColor(kRed);
	vtxTracksSize2->SetMarkerColor(kRed);
	vtxZ2->SetMarkerColor(kRed);
	vtxX2->SetMarkerColor(kRed);
	vtxY2->SetMarkerColor(kRed);
	pt2->SetMarkerColor(kRed);
	ptError2->SetMarkerColor(kRed);
	DCAz2->SetMarkerColor(kRed);
	DCAxy2->SetMarkerColor(kRed);
	eta2->SetMarkerColor(kRed);
	phi2->SetMarkerColor(kRed);
	Chi2n2->SetMarkerColor(kRed);
    numberOfHits2->SetMarkerColor(kRed);
	Algo2->SetMarkerColor(kRed);

	Ntrk2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	nVtx2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	vtxTracksSize2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	vtxZ2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	vtxX2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	vtxY2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	pt2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	ptError2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	DCAz2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	DCAxy2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	eta2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	phi2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	Chi2n2->Scale( 1.0/ (vtxZ2->GetEntries()) );
    numberOfHits2->Scale( 1.0/ (vtxZ2->GetEntries()) );
	Algo2->Scale( 1.0/ (vtxZ2->GetEntries()) );

	eta1->Rebin(10);
	phi1->Rebin(10);

	eta2->Rebin(10);
	phi2->Rebin(10);

	TLegend *w1 = new TLegend(0.20,0.65,0.5,0.80);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    w1->SetTextSize(10);
    w1->SetTextFont(43);
    w1->AddEntry(Ntrk1,"800_pre1,76X_dataRun1_v10","P");
    w1->AddEntry(Ntrk2,"800_pre6,80X_dataRun2_v8","P");


	TCanvas* c1 = makeMultiCanvas("c1","c1",2,1);
	c1->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetLogz(1);
	caloVsCbin1->SetStats(kFALSE);
	caloVsCbin1->GetYaxis()->SetRangeUser(0,2);
	caloVsCbin1->GetXaxis()->SetRangeUser(0,40);
	caloVsCbin1->Draw("colz");

	c1->cd(2);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetLogz(1);
	caloVsCbin2->SetStats(kFALSE);
	caloVsCbin2->GetYaxis()->SetRangeUser(0,2);	
	caloVsCbin2->GetXaxis()->SetRangeUser(0,40);	
	caloVsCbin2->Draw("colz");

	TCanvas* c2 = makeMultiCanvas("c2", "c2", 3,3);
	c2->cd(1);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	nVtx1->Draw("P");
	nVtx2->Draw("Psame");
	w1->Draw("same");

	c2->cd(2);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLogx(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	vtxTracksSize1->Draw("P");
	vtxTracksSize2->Draw("Psame");

	c2->cd(3);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	vtxZ1->Draw("P");
	vtxZ2->Draw("Psame");

	c2->cd(4);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	vtxX1->Draw("P");
	vtxX2->Draw("Psame");

	c2->cd(5);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	vtxY1->Draw("P");
	vtxY2->Draw("Psame");

	c2->cd(6);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLogx(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	pt1->Draw("P");
	pt2->Draw("Psame");

	c2->cd(7);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	eta1->Draw("P");
	eta2->Draw("Psame");

	c2->cd(8);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	phi1->Draw("P");
	phi2->Draw("Psame");

	c2->cd(9);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLogx(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	ptError1->Draw("P");
	ptError2->Draw("Psame");

	TCanvas* c3 = makeMultiCanvas("c3", "c3", 3,2);
	c3->cd(1);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	DCAz1->GetXaxis()->SetRangeUser(0,20);
	DCAz1->Draw("P");
	DCAz2->Draw("Psame");
	w1->Draw("same");

	c3->cd(2);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	DCAxy1->GetXaxis()->SetRangeUser(0,20);	
	DCAxy1->Draw("P");
	DCAxy2->Draw("Psame");

	c3->cd(3);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	Chi2n1->Draw("P");
	Chi2n2->Draw("Psame");

	c3->cd(4);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	Algo1->Draw("P");
	Algo2->Draw("Psame");

	c3->cd(5);
	gPad->SetTicks();
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	numberOfHits1->Draw("P");
	numberOfHits2->Draw("Psame");

	c3->cd(6);
	gPad->SetTicks();
	gPad->SetLogy(1);
	gPad->SetLogx(1);	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	ptError1->Draw("P");
	ptError2->Draw("Psame");



}