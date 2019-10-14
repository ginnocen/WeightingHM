#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TDirectoryFile.h"
#include "TH1D.h"
#include "TCanvas.h"
//#include "AliCFTaskVertexingHF.h"
//#include "AliCFContainer.h"
#include "TLegend.h"
#include <Riostream.h>


void MakeD0AccEff(int bunch=16,int nameCase=0){
    

    TString suffix = "EvWithD";
    TFile *_file16 = TFile::Open("../../../../../../../../newStat_July/MC_july/D0efficiencies/../rootfiles/output_EvWithD_16merged.root");//i pesi sono gli stessi per D0 e Lc
    TFile *_file17 = TFile::Open("../../../../../../../../newStat_July/MC_july/D0efficiencies/../rootfiles/output_EvWithD_17merged.root");
    TFile *_file18 = TFile::Open("../../../../../../../../newStat_July/MC_july/D0efficiencies/../rootfiles/output_EvWithD_18merged.root");
    TH1F *hWeights;
    if(bunch==16) hWeights = (TH1F*)_file16->Get("Weights0");
    else if(bunch==17) hWeights = (TH1F*)_file17->Get("Weights0");
    else if(bunch==18) hWeights = (TH1F*)_file18->Get("Weights0");
    else {
        Printf("Wrong bunch number, check!");
        return;
    }
TString nameFadd[8]={"_central","_CrossPt","_crossRowFind","_l1","_l2","_nclus","_t1","_t2"};
TString namePadd[8]={"","_CrossPt","_crossRowFind","_l1","_l2","_nclus","_t1","_t2"};        

    TString nameP=Form("_prompt_2016%s",namePadd[nameCase].Data());;
    TString nameF=Form("_fromB_13t%s",nameFadd[nameCase].Data());//PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly_fromB_13t_central
//AnalysisResults_13TeV_2016_TOTMC.root
//AnalysisResults_1574_18MC.root
//AnalysisResults_1575_17MC.root
   
    TString fileNameEff;
    if(bunch==16)fileNameEff= Form("../AnalysisResults_13TeV_2016_TOTMC.root",bunch);
    if(bunch==18) fileNameEff = Form("../AnalysisResults_1678_2018.root",bunch);
    if(bunch==17) fileNameEff = Form("../AnalysisResults_1677_2017.root",bunch);

    TString baseDirFilNameEff=Form("PWG3_D2H_CFtaskD0toKpi%s",nameP.Data());//PR o FD
    TString baseContainerNameEff=Form("CFHFccontainer0%s",nameP.Data());
    TString baseDirFilNameBEff=Form("PWG3_D2H_CFtaskD0toKpiKeepDfromBOnly%s",nameF.Data());
    TString baseContainerNameBEff=Form("CFHFccontainer0DfromB%s",nameF.Data());
    const Int_t nPtBins=6;
    Double_t ptLimits[nPtBins+1];
    
    ptLimits[0]=1.;
    ptLimits[1]=2.;
//    ptLimits[1]=3.;
    ptLimits[2]=4.;
//    ptLimit[2]=5.;
    ptLimits[3]=6.;
    ptLimits[4]=8.;
    ptLimits[5]=12.;
    ptLimits[6]=24.;

    const int nMultBins = 4;
    Double_t multLims[nMultBins+1] = {1.,10.,30.,60.,1000};//160 con HMSPD

    Int_t middle=AliCFTaskVertexingHF::kStepAcceptance;//mesoni Ds con tutte le figlie in accettanza
    Int_t numer=AliCFTaskVertexingHF::kStepRecoPID;//mesoni Ds selezionati dopo tagli+PID
    
    
    TFile* fileff = new TFile(fileNameEff.Data());
    cout<<baseDirFilNameEff.Data()<<"  "<<baseDirFilNameBEff.Data()<<endl;
    TDirectoryFile* dfeff=(TDirectoryFile*)fileff->Get(Form("%s",baseDirFilNameEff.Data()));
    if(!dfeff) Printf("dfeff not exists");
    TDirectoryFile* dfeffB=(TDirectoryFile*)fileff->Get(Form("%s",baseDirFilNameBEff.Data()));
    
    AliCFContainer* conteff=(AliCFContainer*)dfeff->Get(Form("%s",baseContainerNameEff.Data()));
    AliCFContainer* conteffB=(AliCFContainer*)dfeffB->Get(Form("%s",baseContainerNameBEff.Data()));
    cout<< " conteffB " << conteffB << " dfeffB " << dfeffB << endl;
    
    TH1D* hpt;
    TH1D* hrec;
    TH1D* hrecB;
    TH1D* hgeneff;
    TH1D* hgeneffB;
    
    TH1D* hnum[nMultBins+1];
    TH1D* hnumB[nMultBins+1];
    TH1D* hden[nMultBins+1];
    TH1D* hdenB[nMultBins+1];

    for(int iMultBin = 0; iMultBin <= nMultBins; iMultBin++) {
        if(iMultBin<nMultBins) {
            hnum[iMultBin]  = new TH1D(Form("hnumMultBin%.0f%.0f",multLims[iMultBin],multLims[iMultBin+1]-1),"",nPtBins,ptLimits);
            hnumB[iMultBin] = new TH1D(Form("hnumBMultBin%.0f%.0f",multLims[iMultBin],multLims[iMultBin+1]-1),"",nPtBins,ptLimits);
            hden[iMultBin]  = new TH1D(Form("hdenMultBin%.0f%.0f",multLims[iMultBin],multLims[iMultBin+1]-1),"",nPtBins,ptLimits);
            hdenB[iMultBin] = new TH1D(Form("hdenBMultBin%.0f%.0f",multLims[iMultBin],multLims[iMultBin+1]-1),"",nPtBins,ptLimits);
        }
        else {//minimum bias
            hnum[nMultBins]  = new TH1D(Form("hnumMultBin%.0f%.0f",multLims[0],multLims[nMultBins]-1),"",nPtBins,ptLimits);
            hnumB[nMultBins] = new TH1D(Form("hnumBMultBin%.0f%.0f",multLims[0],multLims[nMultBins]-1),"",nPtBins,ptLimits);
            hden[nMultBins]  = new TH1D(Form("hdenMultBin%.0f%.0f",multLims[0],multLims[nMultBins]-1),"",nPtBins,ptLimits);
            hdenB[nMultBins] = new TH1D(Form("hdenBMultBin%.0f%.0f",multLims[0],multLims[nMultBins]-1),"",nPtBins,ptLimits);
        }
    }
    
    for(int iPtBin = 0; iPtBin < nPtBins; iPtBin++) {
cout<< " for iPtBin " << iPtBin <<endl;

        //reset range
        conteff->SetRangeUser(0,0.,24.);
        conteffB->SetRangeUser(0,0.,24.);
        //set new range
        cout<< " e poi set new range " << endl;
        conteff->SetRangeUser(0,ptLimits[iPtBin],ptLimits[iPtBin+1]-0.001);
        conteffB->SetRangeUser(0,ptLimits[iPtBin],ptLimits[iPtBin+1]-0.001);
      
    
        for(int iMultBin = 0; iMultBin <= nMultBins; iMultBin++) {
          cout<< " iMultBin "<< iMultBin <<endl;
            
            //reset range multiplicity
            conteff->SetRangeUser(3,0.,1000.);
            conteffB->SetRangeUser(3,0.,1000.);
            //set new range
            
            if(iMultBin<nMultBins) {
                conteff->SetRangeUser(3,multLims[iMultBin],multLims[iMultBin+1]-0.001);
                conteffB->SetRangeUser(3,multLims[iMultBin],multLims[iMultBin+1]-0.001);
            }
            else {
              conteff->SetRangeUser(3,multLims[0],multLims[nMultBins]-0.001);
                conteffB->SetRangeUser(3,multLims[0],multLims[nMultBins]-0.001);
            }
            
            hpt     = (TH1D*)conteff->Project(numer,0);
            Printf("pt lims %f-%f",hpt->GetBinLowEdge(1),hpt->GetBinLowEdge(hpt->GetNbinsX()+1));

            hrec     = (TH1D*)conteff->Project(numer,3);
            hrecB    = (TH1D*)conteffB->Project(numer,3);
            hgeneff  = (TH1D*)conteff->Project(middle,3);
            hgeneffB = (TH1D*)conteffB->Project(middle,3);
            
            hrec->Sumw2();
            hrecB->Sumw2();
            hgeneff->Sumw2();
            hgeneffB->Sumw2();
            Printf("mult %f-%f",hrec->GetBinLowEdge(1),hrec->GetBinLowEdge(hrec->GetNbinsX()+1));

            for(int j = 1; j <= hrec->GetNbinsX(); j++) {
                double mult = hrec->GetBinCenter(j)-0.5;
                cout<< " mult " << mult << " original " << hrec->GetBinCenter(j) <<endl;
                double w = hWeights->GetBinContent(hWeights->FindBin(mult));
                double y = hrec->GetBinContent(j);
                
                hrec->SetBinContent(j,y*w);
                hrec->SetBinError(j,TMath::Abs(y*w));
                y = hrecB->GetBinContent(j);
                hrecB->SetBinContent(j,y*w);
                hrecB->SetBinError(j,TMath::Abs(y*w));
                y = hgeneff->GetBinContent(j);
                hgeneff->SetBinContent(j,y*w);
                hgeneff->SetBinError(j,TMath::Abs(y*w));
                y = hgeneffB->GetBinContent(j);
                hgeneffB->SetBinContent(j,y*w);
                hgeneffB->SetBinError(j,TMath::Abs(y*w));
            }
            hnum[iMultBin]->SetBinContent(iPtBin+1,hrec->Integral());
            hnumB[iMultBin]->SetBinContent(iPtBin+1,hrecB->Integral());
            hden[iMultBin]->SetBinContent(iPtBin+1,hgeneff->Integral());
            hdenB[iMultBin]->SetBinContent(iPtBin+1,hgeneffB->Integral());

        }
    }
    

    TFile* outfil=new TFile(Form("D0Eff_year%d_1_1000_Ntrkl_%s_%s.root",bunch,suffix.Data(),nameFadd[nameCase].Data()),"recreate");
    for(int iMultBin = 0; iMultBin <= nMultBins; iMultBin++) {
        hnum[iMultBin]->Write();
        hnumB[iMultBin]->Write();
        hden[iMultBin]->Write();
        hdenB[iMultBin]->Write();
    }
    outfil->Close();
 
    
}


