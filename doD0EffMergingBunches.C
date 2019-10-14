#include "TFile.h"
#include "TString.h"
#include "TH1D.h"

void doDsEffMergingBunches(int nameCase=0) {
    
    TString suffix = "EvWithD";

    
    TString fileNameAcc = "D0Acc_FONLL5mio.root";
    TFile *fileacc = TFile::Open(fileNameAcc.Data());
    TH1D *hAcc = (TH1D*)fileacc->Get("hAcc");

    TFile *_file0;
    const int nmult = 5;
    TString mult[nmult] = {"19","1029","3059","60999","1999"};
    TH1D *hnum;
    TH1D *hnumB;
    TH1D *hden;
    TH1D *hdenB;
    
    TH1D *hNumer[nmult];
    TH1D *hNumerB[nmult];
    TH1D *hDenom[nmult];
    TH1D *hDenomB[nmult];
    TH1D *hEff[nmult];
    TH1D *hEffB[nmult];
TString nameFadd[8]={"_central","_CrossPt","_crossRowFind","_l1","_l2","_nclus","_t1","_t2"};
TString namePadd[8]={"","_CrossPt","_crossRowFind","_l1","_l2","_nclus","_t1","_t2"};         
    
    for(int ibunch=16; ibunch<=18; ibunch++) {
        _file0 = TFile::Open(Form("D0Eff_year%d_1_1000_Ntrkl_%s_%s.root",ibunch,suffix.Data(),nameFadd[nameCase].Data()));
//        D0Eff_year16_1_100_Ntrkl_EvWithD.root
       _file0->ls(); 
        for(int jmult=0; jmult<nmult; jmult++) {
          cout<<  " nome mult " <<mult[jmult] <<endl;
            hnum  = (TH1D*)_file0->Get(Form("hnumMultBin%s",mult[jmult].Data()));
            hnumB = (TH1D*)_file0->Get(Form("hnumBMultBin%s",mult[jmult].Data()));
            hden  = (TH1D*)_file0->Get(Form("hdenMultBin%s",mult[jmult].Data()));
            hdenB = (TH1D*)_file0->Get(Form("hdenBMultBin%s",mult[jmult].Data()));
            
            hnum->Sumw2();
            hnumB->Sumw2();
            hden->Sumw2();
            hdenB->Sumw2();
            
            if(ibunch==16) {
                hNumer[jmult] = (TH1D*)hnum->Clone(Form("hNumMult%s",mult[jmult].Data()));
                hNumerB[jmult] = (TH1D*)hnumB->Clone(Form("hNumB%s",mult[jmult].Data()));
                hDenom[jmult] = (TH1D*)hden->Clone(Form("hDen%s",mult[jmult].Data()));
                hDenomB[jmult] = (TH1D*)hdenB->Clone(Form("hDenB%s",mult[jmult].Data()));
                hNumer[jmult]->Sumw2();
                hNumerB[jmult]->Sumw2();
                hDenom[jmult]->Sumw2();
                hDenomB[jmult]->Sumw2();
            }
            else {
                hNumer[jmult]->Add(hnum);
                hNumerB[jmult]->Add(hnumB);
                hDenom[jmult]->Add(hden);
                hDenomB[jmult]->Add(hdenB);
            }
            
        }
        
    }
    
    
    for(int jmult=0; jmult<nmult; jmult++) {
        hEff[jmult] = 0x0;
        hEffB[jmult] = 0x0;
        
        hEff[jmult]  = (TH1D*)hNumer[jmult]->Clone(Form("hEff_C_Mult%s",mult[jmult].Data()));
        hEffB[jmult] = (TH1D*)hNumerB[jmult]->Clone(Form("hEff_B_Mult%s",mult[jmult].Data()));
   
        hEff[jmult]->Divide(hNumer[jmult],hDenom[jmult],1,1,"B");
        hEffB[jmult]->Divide(hNumerB[jmult],hDenomB[jmult],1,1,"B");
        hEff[jmult]->SetName("hEff_C");//Form("hEff_C_Mult%s",mult[jmult].Data()));
        hEffB[jmult]->SetName("hEff_B");//Form("hEff_B_Mult%s",mult[jmult].Data()));

        for(int ibin=1;ibin<=hEff[jmult]->GetNbinsX();ibin++){
            double ptcenter = hEff[jmult]->GetBinCenter(ibin);
            int accbin = hAcc->FindBin(ptcenter);
            hEff[jmult]->SetBinContent(ibin,hEff[jmult]->GetBinContent(ibin)*hAcc->GetBinContent(accbin));
            hEffB[jmult]->SetBinContent(ibin,hEffB[jmult]->GetBinContent(ibin)*hAcc->GetBinContent(accbin));
        }
        
        TFile* outfil=new TFile(Form("D0AccEff_161718_wNtrklWeights_%s_%s_%s.root",mult[jmult].Data(),suffix.Data(),nameFadd[nameCase].Data()),"recreate");
        hEff[jmult]->Write();
        hEffB[jmult]->Write();
        outfil->Close();

    }
    
    
}
