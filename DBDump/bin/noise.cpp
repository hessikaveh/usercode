#include "usercode/DBDump/interface/NoiseDumper.h"
#include "usercode/DBDump/interface/NoiseLaserDumper.h"

#include "CondFormats/BeamSpotObjects/interface/BeamSpotObjects.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/EcalObjects/interface/EcalClusterLocalContCorrParameters.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/EcalObjects/interface/EcalLaserAlphas.h"
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalPulseShapes.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLinearizationConst.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLutGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGLutIdMap.h"
#include "CondFormats/EcalObjects/interface/EcalTPGPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightGroup.h"
#include "CondFormats/EcalObjects/interface/EcalTPGWeightIdMap.h"
#include "CondFormats/EcalObjects/interface/EcalTPGSlidingWindow.h"
#include "CondFormats/EcalObjects/interface/EcalTPGSpike.h"
#include "CondFormats/ESObjects/interface/ESEEIntercalibConstants.h"
#include "CondFormats/ESObjects/interface/ESGain.h"
#include "CondFormats/ESObjects/interface/ESIntercalibConstants.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include <TROOT.h>
#include <TStyle.h>
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TF1.h"

//TH2F* dum_pedEB;
//TH2F* dum_pedEE;
//TH2F* dum_icEB;
//TH2F* dum_icEE;
//TH2F* dum_alphaEB;
//TH2F* dum_alphaEE;
//TH2F* dum_laserEB;
//TH2F* dum_laserEE;
TH2F* histoEB_IC = new TH2F ("histoEB_IC", "IC" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
TH2F* histoEE_IC = new TH2F ("histoEE_IC", "IC" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
TH2F* histoEB_alpha = new TH2F ("histoEB_alpha", "#alpha" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
TH2F* histoEE_alpha = new TH2F ("histoEE_alpha", "#alpha" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
TH2F* histoEE_p1 = new TH2F ("histoEE_p1", "Relative response to laser" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
TH2F* histoEB_p1 = new TH2F ("histoEB_p1", "Relative response to laser" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
TH2F* histoEB_G12rms = new TH2F ("histoEB_G12rms", "G12 rms"       ,  360, 0.5, 360.5,  171, -85.5, 85.5);
TH2F* histoEE_G12rms = new TH2F ("histoEE_G12rms", "G12 rms"       ,  200, 0.5, 200.5,  100, 0.5, 100.5);


void plotIC(cond::NoiseDumper<EcalIntercalibConstants> ICmain) {

    gStyle->SetOptStat(0);

    TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);

    TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);


    //
    // input file format:
    //
    // ix/ieta iy/iphi   iz/0       IC        rawId
    //
    // -85    2    0        0.983292 838904322
    //



    std::vector<int> ix_ieta;
    std::vector<int> iy_iphi;
    std::vector<int> iz;

    std::vector<float> IC;



    //std::cout <<ICmain.NoiseVariables.size() << std::endl;
    cond::NoiseDumper<EcalIntercalibConstants>::NoiseStruct a;

    for(int i =0 ; i < (int) ICmain.NoiseVariables.size(); ++i)
    {
        a = ICmain.NoiseVariables.at(i);

        std::cout << a.IOVBegin  << " IOV "<< a.IOVEnd  << " IC: "  <<a.noise_ic.size() << std::endl;
        for(int i =0 ; i < (int) a.noise_ic.size(); ++i)
        {
            ix_ieta.push_back(a.noise_ic.at(i).ix);
            iy_iphi.push_back(a.noise_ic.at(i).iy);
            iz.push_back(a.noise_ic.at(i).iz);
            IC.push_back(a.noise_ic.at(i).ic);


            //            std::cout << "time: " << a.time << "   " <<a.noise_ic.at(i).ix  << " Det "<< a.noise_ic.at(i).iy   << " IC: "  <<a.noise_ic.at(i).ic << std::endl;

        }

    }





    std::cout << " ix_ieta.size() = " << ix_ieta.size() << std::endl;

    if (ix_ieta.size() > 75848) {
        std::cout << " Attention: you appended the tag twice!" << std::endl;
    }


    //---- EB ----

    histoEB_IC->GetXaxis()->SetTitle("i#phi");
    histoEB_IC->GetYaxis()->SetTitle("i#eta");


    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            histoEB_IC->Fill(iy_iphi.at(iter), ix_ieta.at(iter), IC.at(iter) );

        }
    }

    ccEB->cd();
    histoEB_IC->Draw("colz");
    ccEB->SaveAs("IC_EB.png");
    ccEB->SaveAs("IC_EB.root");






    //---- EE ----



    histoEE_IC->GetXaxis()->SetTitle("x");
    histoEE_IC->GetYaxis()->SetTitle("y");

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) != 0) {

            histoEE_IC->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), IC.at(iter) );

        }
    }

    ccEE->cd();
    histoEE_IC->Draw("colz");
    ccEE->SaveAs("IC_EE.png");
    ccEE->SaveAs("IC_EE.root");


    // ---- read iring definition from file
    //        (ix, iy, iz) -> ring
    //


    std::map < std::pair<int, int> , int > iring_map_plus;
    std::map < std::pair<int, int> , int > iring_map_minus;

    std::ifstream fileEEring ("../data/eerings.dat");

    std::string buffer;
//    int num;
    if (!fileEEring.is_open()) {
        std::cerr << "** ERROR: Can't open for input" << std::endl;
        //return false;
    }

    while(!fileEEring.eof()) {
        getline(fileEEring,buffer);
        if (buffer != ""){ ///---> save from empty line at the end!
            int ix;
            int iy;
            int iz;
            int iring;

            //       std::cout << " buffer = " << buffer << std::endl;

            std::stringstream line( buffer );
            line >> ix;
            line >> iy;
            line >> iz;
            line >> iring;

            std::pair<int, int> ixiy (ix, iy);
            //       if (iz>0) iring_map_plus  [ixiy] = 38 - iring;
            //       else      iring_map_minus [ixiy] = 38 - iring;

            if (iz>0) iring_map_plus  [ixiy] = iring;
            else      iring_map_minus [ixiy] = iring;

        }
    }


    TH2F *iring_map = (TH2F*) histoEE_p1->Clone("iring_map");

    for(int ix=1; ix<=histoEE_p1->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_p1->GetNbinsY(); iy++){
            int iring = -99;
            std::pair<int, int> ixiy (ix%100, iy%100);
            if ((1*(ix>=100)-1*(ix<100)) > 0) {
                //         iring = 38 - iring_map_plus [ixiy];
                iring = iring_map_plus [ixiy];
            }
            else {
                //         iring = 38 - iring_map_minus [ixiy];
                iring = iring_map_minus [ixiy];
            }

//            std::cout << " ix,iy,ring = " << ix%100 << "  " << iy%100 << "  " << iring << std::endl;

            iring_map->SetBinContent(ix,iy,iring);

        }
    }


    // TCanvas* ccRingMap = new TCanvas ("ccRingMap","",800,600);
    iring_map->Draw("colz");

    std::vector<int> d_ix_ieta;
    std::vector<int> d_iy_iphi;
    std::vector<int> d_iz;

    std::vector<float> v_Product;



    std::vector<float> ringPlus_max_Product;      std::vector<float> ringMinus_max_Product;
    std::vector<float> ringPlus_min_Product;      std::vector<float> ringMinus_min_Product;

    std::vector<float> ringPlus_Product;          std::vector<float> ringMinus_Product;
    std::vector<float> ringPlus_ProductSpread;    std::vector<float> ringMinus_ProductSpread;
    std::vector<int> ringPlus_count;              std::vector<int> ringMinus_count;


    for (int iter = 0; iter < (50-11+3); iter++) {
        ringPlus_max_Product.push_back(0);        ringMinus_max_Product.push_back(0);
        ringPlus_min_Product.push_back(99);       ringMinus_min_Product.push_back(99);
        ringPlus_Product.push_back(0);            ringMinus_Product.push_back(0);
        ringPlus_ProductSpread.push_back(0);      ringMinus_ProductSpread.push_back(0);
        ringPlus_count.push_back(0);              ringMinus_count.push_back(0);
    }


    for(int ix=1; ix<=histoEE_IC->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_IC->GetNbinsY(); iy++){

            d_ix_ieta.push_back(ix%100);
            d_iy_iphi.push_back(iy%100);
            d_iz.push_back     (1*(ix>=100)-1*(ix<100));

            if (histoEE_IC->GetBinContent(ix,iy) > 0 ) {
                v_Product.push_back( histoEE_IC->GetBinContent(ix,iy) );
            }
            else {
                v_Product.push_back( -99 );
            }
        }
    }

    std::cout << d_ix_ieta.size() << std::endl;
    for (int iter = 0; iter < (int) d_ix_ieta.size(); iter++) {
        if (d_iz.at(iter) != 0) {

            float dx = d_ix_ieta.at(iter) - 50;
            float dy = d_iy_iphi.at(iter) - 50;

            float ring = sqrt( dx*dx + dy*dy );

            int iring = round(ring) - 11;  //---- 12 [ = (62 - 50 - 1) from the 2D plot] is the first ring

            iring = -99;
            std::pair<int, int> ixiy (d_ix_ieta.at(iter), d_iy_iphi.at(iter));

            //       std::map<std::pair<int, int>,int>::const_iterator it = iring_map_plus.find(ixiy);
            //       if (it==iring_map_plus.end() ) {
            //         continue;
            //       }
            //

            if (d_iz.at(iter) > 0) {
                iring = iring_map_plus [ixiy];
            }
            else {
                iring = iring_map_minus [ixiy];
            }

            //       std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;

            //       if (iring < 0 ) continue;

            //       if (iring > 38 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;
            //       if (iring > 37 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << std::endl;

            if (v_Product.at(iter) > 0 ) if (iring == 0 )  std::cout << " ix, iy, iring = " << d_ix_ieta.at(iter) << "  " <<  d_iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << " tot = " << ringPlus_count.at(iring) << std::endl;


            if (v_Product.at(iter) > 0 ) {
                if (iring > (50-11+2) || iring < 0) std::cout << " what ?!?   iring = " << iring << " dx = " << dx << " dy = " << dy << " :::: ix = " << d_ix_ieta.at(iter) << "  iy = " << d_iy_iphi.at(iter) << " prod = " << v_Product.at(iter) << std::endl;

                if (d_iz.at(iter) > 0) {

                    if (ringPlus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringPlus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringPlus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringPlus_min_Product.at(iring)   =  v_Product.at(iter);
                    }

                    ringPlus_Product.at(iring)       = ringPlus_Product.at(iring) + v_Product.at(iter);
                    ringPlus_ProductSpread.at(iring) = ringPlus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringPlus_count.at(iring)         = ringPlus_count.at(iring) + 1 ;

                }
                else {

                    if (ringMinus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringMinus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringMinus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringMinus_min_Product.at(iring)   =  v_Product.at(iter);
                    }

                    ringMinus_Product.at(iring)       = ringMinus_Product.at(iring) + v_Product.at(iter);
                    ringMinus_ProductSpread.at(iring) = ringMinus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringMinus_count.at(iring)         = ringMinus_count.at(iring) + 1 ;
                }
            }

        }
    }
//    std::cout << "Is it here now?!??!" << std::endl;

    for (int iring = 0; iring < (50-11+3); iring++) {
        //     std::cout << " ringPlus_ProductSpread.at(" << iring << ") = " << ringPlus_ProductSpread.at(iring) << " ---> " << ringPlus_ProductSpread.at(iring)  << " - " <<  ringPlus_Product.at(iring)  * ringPlus_Product.at(iring) << std::endl;
        ringMinus_ProductSpread.at(iring) = sqrt(ringMinus_ProductSpread.at(iring) -  ringMinus_Product.at(iring) / ringMinus_count.at(iring) * ringMinus_Product.at(iring)) / ringMinus_count.at(iring);
        ringPlus_ProductSpread.at(iring)  = sqrt(ringPlus_ProductSpread.at(iring)  -  ringPlus_Product.at(iring)  / ringPlus_count.at(iring)  * ringPlus_Product.at(iring))  / ringPlus_count.at(iring);
    }




    TGraph *gr_EEPlus_IC = new TGraph();     TGraph *gr_EEMinus_IC = new TGraph();


    for (int iter = 0; iter < (39); iter++) {

        gr_EEPlus_IC-> SetPoint (iter, iter,   ringPlus_count.at(iter)  ?  ringPlus_Product.at(iter)  / ringPlus_count.at(iter)  : 0 ) ;
        gr_EEMinus_IC-> SetPoint (iter,  iter,  ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;

    }


    //---- style ----

    gr_EEPlus_IC->SetMarkerSize  (1);                           gr_EEMinus_IC->SetMarkerSize  (1);
    gr_EEPlus_IC->SetMarkerStyle (24);                          gr_EEMinus_IC->SetMarkerStyle (22);
    gr_EEPlus_IC->SetMarkerColor (kRed);                        gr_EEMinus_IC->SetMarkerColor (kRed);
    gr_EEPlus_IC->SetLineWidth (1);                             gr_EEMinus_IC->SetLineWidth (1);
    gr_EEPlus_IC->SetLineColor (kRed);                          gr_EEMinus_IC->SetLineColor (kRed);


    //---- style (end) ----


    TCanvas* ccRing = new TCanvas ("ccRing","",800,600);

    gr_EEPlus_IC->Draw("APL");
    gr_EEMinus_IC->Draw("PL");

    gr_EEPlus_IC->GetYaxis()->SetTitle("<IC>");
    gr_EEPlus_IC->GetXaxis()->SetTitle("iRing");

    ccRing->SaveAs("IC_EE_perRing.png");
    ccRing->SaveAs("IC_EE_perRing.root");



    //
    //   iring = ieta in EB
    //

    std::vector<float> EBringPlus_IC;          std::vector<float> EBringMinus_IC;
    std::vector<int> EBringPlus_ICcount;       std::vector<int> EBringMinus_ICcount;

    for (int iter = 0; iter < 85; iter++) {
        EBringPlus_IC.push_back(0);            EBringMinus_IC.push_back(0);
        EBringPlus_ICcount.push_back(0);       EBringMinus_ICcount.push_back(0);
    }

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            int iEBring = abs(ix_ieta.at(iter)) - 1 ;


            if (iEBring > 84 || iEBring < 0) std::cout << " what ?!?   iEBring = " << iEBring << std::endl;

            if (ix_ieta.at(iter) > 0) {

                EBringPlus_IC.at(iEBring) = EBringPlus_IC.at(iEBring) + IC.at(iter);
                EBringPlus_ICcount.at(iEBring) =  EBringPlus_ICcount.at(iEBring) + 1 ;

            }
            else {

                EBringMinus_IC.at(iEBring) = EBringMinus_IC.at(iEBring) + IC.at(iter);
                EBringMinus_ICcount.at(iEBring) =  EBringMinus_ICcount.at(iEBring) + 1 ;

            }

        }
    }


    TGraph *gr_EBPlus_IC = new TGraph();     TGraph *gr_EBMinus_IC = new TGraph();

    for (int iter = 0; iter < 85; iter++) {

        gr_EBPlus_IC-> SetPoint (iter,    iter+1,   EBringPlus_ICcount.at(iter) ? EBringPlus_IC.at(iter) / EBringPlus_ICcount.at(iter) : 0 ) ;
        gr_EBMinus_IC-> SetPoint (iter,    iter+1,   EBringMinus_ICcount.at(iter) ?  EBringMinus_IC.at(iter) / EBringMinus_ICcount.at(iter) : 0 ) ;

    }


    //---- style ----

    gr_EBPlus_IC->SetMarkerSize  (1);                           gr_EBMinus_IC->SetMarkerSize  (1);
    gr_EBPlus_IC->SetMarkerStyle (24);                          gr_EBMinus_IC->SetMarkerStyle (22);
    gr_EBPlus_IC->SetMarkerColor (kRed);                        gr_EBMinus_IC->SetMarkerColor (kRed);
    gr_EBPlus_IC->SetLineWidth (1);                             gr_EBMinus_IC->SetLineWidth (1);
    gr_EBPlus_IC->SetLineColor (kRed);                          gr_EBMinus_IC->SetLineColor (kRed);


    //---- style (end) ----


    TCanvas* ccRingEB = new TCanvas ("ccRingEB","",800,600);

    gr_EBPlus_IC->Draw("APL");
    gr_EBMinus_IC->Draw("PL");

    gr_EBPlus_IC->GetYaxis()->SetTitle("<IC>");
    gr_EBPlus_IC->GetXaxis()->SetTitle("i#eta");

    ccRingEB->SaveAs("IC_EB_perRing.png");
    ccRingEB->SaveAs("IC_EB_perRing.root");


    TFile fileout("fileOutIC.root","RECREATE");
    fileout.cd();
    histoEB_IC->Write();
    histoEE_IC->Write();
}

void plotAlpha(cond::NoiseDumper<EcalLaserAlphas> Alpha) {

    gStyle->SetOptStat(0);

    TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);

    TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);


    //
    // input file format:
    //
    // ix/ieta iy/iphi   iz/0       alpha        rawId
    //
    // -85    2    0       1.520000 838904322
    //



    std::vector<int> ix_ieta;
    std::vector<int> iy_iphi;
    std::vector<int> iz;

    std::vector<float> alpha;
    cond::NoiseDumper<EcalLaserAlphas>::NoiseStruct a;
    for(int i =0 ; i < (int) Alpha.NoiseVariables.size(); ++i)
    {
        a = Alpha.NoiseVariables.at(i);

        std::cout << a.IOVBegin  << " IOV "<< a.IOVEnd  << " alpha " <<a.noise_ic.size() << std::endl;
        for(int i =0 ; i < (int) a.noise_ic.size(); ++i)
        {

            ix_ieta.push_back(a.noise_ic.at(i).ix);
            iy_iphi.push_back(a.noise_ic.at(i).iy);
            iz.push_back(a.noise_ic.at(i).iz);
            alpha.push_back(a.noise_ic.at(i).ic);

            // std::cout << "time: " << a.time << "   "<< a.noise_ic.at(i).ix  << " Det "<< a.noise_ic.at(i).iy   << " alpha: "  <<a.noise_ic.at(i).ic << std::endl;

        }


    }


    std::cout << " ix_ieta.size() = " << ix_ieta.size() << std::endl;

    if (ix_ieta.size() > 75848) {
        std::cout << " Attention: you appended the tag twice!" << std::endl;
    }


    //---- EB ----

    histoEB_alpha->GetXaxis()->SetTitle("i#phi");
    histoEB_alpha->GetYaxis()->SetTitle("i#eta");


    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            histoEB_alpha->Fill(iy_iphi.at(iter), ix_ieta.at(iter), alpha.at(iter) );

        }
    }

    ccEB->cd();
    histoEB_alpha->Draw("colz");
    ccEB->SaveAs("Alpha_EB.png");
    ccEB->SaveAs("Alpha_EB.root");






    //---- EE ----



    histoEE_alpha->GetXaxis()->SetTitle("x");
    histoEE_alpha->GetYaxis()->SetTitle("y");

    for (int iter = 0; iter <(int) ix_ieta.size(); iter++) {
        if (iz.at(iter) != 0) {

            histoEE_alpha->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), alpha.at(iter) );

        }
    }

    ccEE->cd();
    histoEE_alpha->Draw("colz");
    ccEE->SaveAs("Alphas_EE.png");
    ccEE->SaveAs("Alphas_EE.root");
    //  dum_alphaEE = (TH2F*) histoEE_alpha->Clone();
    //  dum_alphaEB = (TH2F*) histoEB_alpha->Clone();

    TFile fileout("fileOutAlpha.root","RECREATE");
    fileout.cd();
    histoEB_alpha->Write();
    histoEE_alpha->Write();
}


void plotPedestals(cond::NoiseDumper<EcalPedestals> Ped) {

    gStyle->SetOptStat(0);

    TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);
    ccEB->Divide(3,2);

    TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);
    ccEE->Divide(3,2);

    TCanvas* ccEB2 = new TCanvas ("ccEB2","",1600,600);
    TCanvas* ccEE2 = new TCanvas ("ccEE2","",1600,600);


    //
    //    ------------------------------------
    //        G12-ped  |  G6-ped  |  G1-ped
    //    ------------------------------------
    //        G12-rms  |  G6-rms  |  G1-rms
    //    ------------------------------------
    //

    // input file format:
    //
    // ix/ieta iy/iphi   iz/0       G12-ped  G12-rms     G6-ped  G6-rms      G1-ped  G1-rms      rawId
    //



    std::vector<int> ix_ieta;
    std::vector<int> iy_iphi;
    std::vector<int> iz;

    std::vector<float> G12ped;
    std::vector<float> G6ped;
    std::vector<float> G1ped;

    std::vector<float> G12rms;
    std::vector<float> G6rms;
    std::vector<float> G1rms;

    cond::NoiseDumper<EcalPedestals>::NoiseStruct c;
    for(int i =0 ; i < (int) Ped.NoiseVariables.size(); ++i)
    {
        c = Ped.NoiseVariables.at(i);

        std::cout << c.IOVBegin  << " IOV "<< c.IOVEnd  << " alpha " <<c.noise_ped.size() << std::endl;
        for(int i =0 ; i < (int) c.noise_ped.size(); ++i)
        {
            ix_ieta.push_back(c.noise_ped.at(i).ix);
            iy_iphi.push_back(c.noise_ped.at(i).iy);
            iz.push_back(c.noise_ped.at(i).iz);

            G12ped.push_back(c.noise_ped.at(i).mean1);
            G6ped.push_back(c.noise_ped.at(i).mean2);
            G1ped.push_back(c.noise_ped.at(i).mean3);

            G12rms.push_back(c.noise_ped.at(i).rms1);
            G6rms.push_back(c.noise_ped.at(i).rms2);
            G1rms.push_back(c.noise_ped.at(i).rms3);

            // std::cout << "time: " << c.time << "   " << c.noise_ped.at(i).ix  << " Det "<< c.noise_ped.at(i).iy   << " mean1: "  <<c.noise_ped.at(i).mean1 << " rms1: " << c.noise_ped.at(i).rms1 << std::endl;

        }

    }


    std::cout << " ix_ieta.size() = " << ix_ieta.size() << std::endl;

    if (ix_ieta.size() > 75848) {
        std::cout << " Attention: you appended the tag twice!" << std::endl;
    }


    //---- EB ----

    TH2F* histoEB_G12ped = new TH2F ("histoEB_G12ped", "G12 pedestals" ,  360, 0.5, 360.5,  171, -85.5, 85.5);
    TH2F* histoEB_G6ped  = new TH2F ("histoEB_G6ped",  "G6 pedestals"  ,  360, 0.5, 360.5,  171, -85.5, 85.5);
    TH2F* histoEB_G1ped  = new TH2F ("histoEB_G1ped",  "G1 pedestals"  ,  360, 0.5, 360.5,  171, -85.5, 85.5);

    TH2F* histoEB_G6rms  = new TH2F ("histoEB_G6rms",  "G6 rms"        ,  360, 0.5, 360.5,  171, -85.5, 85.5);
    TH2F* histoEB_G1rms  = new TH2F ("histoEB_G1rms",  "G1 rms"        ,  360, 0.5, 360.5,  171, -85.5, 85.5);


    histoEB_G12ped->GetXaxis()->SetTitle("i#phi");
    histoEB_G12ped->GetYaxis()->SetTitle("i#eta");
    histoEB_G1ped->GetXaxis()->SetTitle("i#phi");
    histoEB_G1ped->GetYaxis()->SetTitle("i#eta");
    histoEB_G6ped->GetXaxis()->SetTitle("i#phi");
    histoEB_G6ped->GetYaxis()->SetTitle("i#eta");

    histoEB_G12rms->GetXaxis()->SetTitle("i#phi");
    histoEB_G12rms->GetYaxis()->SetTitle("i#eta");
    histoEB_G1rms->GetXaxis()->SetTitle("i#phi");
    histoEB_G1rms->GetYaxis()->SetTitle("i#eta");
    histoEB_G6rms->GetXaxis()->SetTitle("i#phi");
    histoEB_G6rms->GetYaxis()->SetTitle("i#eta");


    //   histoEB_G12ped->GetZaxis()->SetRangeUser( 0.5, 3.0 );
    //   histoEB_G12rms->GetZaxis()->SetRangeUser( 0.5, 3.0 );
    //   histoEB_G1ped->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEB_G1rms->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEB_G6ped->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEB_G6rms->GetZaxis() ->SetRangeUser( 0.5, 3.0 );


    for (int iter = 0; iter <(int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            histoEB_G12ped->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G12ped.at(iter) );
            histoEB_G12rms->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G12rms.at(iter) );

            histoEB_G6ped ->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G6ped.at(iter)  );
            histoEB_G6rms ->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G6rms.at(iter)  );

            histoEB_G1ped ->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G1ped.at(iter)  );
            histoEB_G1rms ->Fill(iy_iphi.at(iter), ix_ieta.at(iter), G1rms.at(iter)  );

        }
    }


    ccEB->cd(1);
    histoEB_G12ped->Draw("colz");
    ccEB->cd(2);
    histoEB_G6ped->Draw("colz");
    ccEB->cd(3);
    histoEB_G1ped->Draw("colz");


    ccEB->cd(4);
    histoEB_G12rms->Draw("colz");
    ccEB->cd(5);
    histoEB_G6rms->Draw("colz");
    ccEB->cd(6);
    histoEB_G1rms->Draw("colz");

    ccEB->SaveAs("pedestalsEB.png");
    ccEB->SaveAs("pedestalsEB.root");

    ccEB2->cd();
    //histoEB_G12rms->SetMinimum(1);
    //histoEB_G12rms->SetMaximum(2.5);
    histoEB_G12rms->Draw("colz");
    ccEB2->SaveAs("pedestalsRMS_G12_EB.png");
    ccEB2->SaveAs("pedestalsRMS_G12_EB.root");







    //---- EE ----

    TH2F* histoEE_G12ped = new TH2F ("histoEE_G12ped", "G12 pedestals" ,  200, 0.5, 200.5,  100, 0.5, 100.5);
    TH2F* histoEE_G6ped  = new TH2F ("histoEE_G6ped",  "G6 pedestals"  ,  200, 0.5, 200.5,  100, 0.5, 100.5);
    TH2F* histoEE_G1ped  = new TH2F ("histoEE_G1ped",  "G1 pedestals"  ,  200, 0.5, 200.5,  100, 0.5, 100.5);

    TH2F* histoEE_G6rms  = new TH2F ("histoEE_G6rms",  "G6 rms"        ,  200, 0.5, 200.5,  100, 0.5, 100.5);
    TH2F* histoEE_G1rms  = new TH2F ("histoEE_G1rms",  "G1 rms"        ,  200, 0.5, 200.5,  100, 0.5, 100.5);


    histoEE_G12ped->GetXaxis()->SetTitle("x");
    histoEE_G12ped->GetYaxis()->SetTitle("y");
    histoEE_G1ped->GetXaxis()->SetTitle("x");
    histoEE_G1ped->GetYaxis()->SetTitle("y");
    histoEE_G6ped->GetXaxis()->SetTitle("x");
    histoEE_G6ped->GetYaxis()->SetTitle("y");

    histoEE_G12rms->GetXaxis()->SetTitle("x");
    histoEE_G12rms->GetYaxis()->SetTitle("y");
    histoEE_G1rms->GetXaxis()->SetTitle("x");
    histoEE_G1rms->GetYaxis()->SetTitle("y");
    histoEE_G6rms->GetXaxis()->SetTitle("x");
    histoEE_G6rms->GetYaxis()->SetTitle("y");


    //   histoEE_G12ped->GetZaxis()->SetRangeUser( 0.5, 3.0 );
    //   histoEE_G12rms->GetZaxis()->SetRangeUser( 0.5, 3.0 );
    //   histoEE_G1ped->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEE_G1rms->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEE_G6ped->GetZaxis() ->SetRangeUser( 0.5, 3.0 );
    //   histoEE_G6rms->GetZaxis() ->SetRangeUser( 0.5, 3.0 );


    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) != 0) {

            histoEE_G12ped->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G12ped.at(iter) );
            histoEE_G12rms->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G12rms.at(iter) );

            histoEE_G6ped ->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G6ped.at(iter)  );
            histoEE_G6rms ->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G6rms.at(iter)  );

            histoEE_G1ped ->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G1ped.at(iter)  );
            histoEE_G1rms ->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), G1rms.at(iter)  );

        }
    }


    ccEE->cd(1);
    histoEE_G12ped->Draw("colz");
    ccEE->cd(2);
    histoEE_G6ped->Draw("colz");
    ccEE->cd(3);
    histoEE_G1ped->Draw("colz");


    ccEE->cd(4);
    histoEE_G12rms->Draw("colz");
    ccEE->cd(5);
    histoEE_G6rms->Draw("colz");
    ccEE->cd(6);
    histoEE_G1rms->Draw("colz");

    ccEE->SaveAs("pedestalsEE.png");
    ccEE->SaveAs("pedestalsEE.root");

    ccEE2->cd();
    //histoEE_G12rms->SetMinimum(1.5);
    //histoEE_G12rms->SetMaximum(2.5);
    histoEE_G12rms->Draw("colz");
    ccEE2->SaveAs("pedestalsRMS_G12_EE.png");
    ccEE2->SaveAs("pedestalsRMS_G12_EE.root");

    // ---- read iring definition from file
    //        (ix, iy, iz) -> ring
    //


    std::map < std::pair<int, int> , int > iring_map_plus;
    std::map < std::pair<int, int> , int > iring_map_minus;

    std::ifstream fileEEring ("../usercode/data/eerings.dat");

    std::string buffer;
//    int num;
    if (!fileEEring.is_open()) {
        std::cerr << "** ERROR: Can't open for input" << std::endl;
        //return false;
    }

    while(!fileEEring.eof()) {
        getline(fileEEring,buffer);
        if (buffer != ""){ ///---> save from empty line at the end!
            int ix;
            int iy;
            int iz;
            int iring;

            //       std::cout << " buffer = " << buffer << std::endl;

            std::stringstream line( buffer );
            line >> ix;
            line >> iy;
            line >> iz;
            line >> iring;

            std::pair<int, int> ixiy (ix, iy);
            //       if (iz>0) iring_map_plus  [ixiy] = 38 - iring;
            //       else      iring_map_minus [ixiy] = 38 - iring;

            if (iz>0) iring_map_plus  [ixiy] = iring;
            else      iring_map_minus [ixiy] = iring;

        }
    }

    TH2F *iring_map = (TH2F*) histoEE_G12ped->Clone("iring_map");

    for(int ix=1; ix<=histoEE_G12ped->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_G12ped->GetNbinsY(); iy++){
            int iring = -99;
            std::pair<int, int> ixiy (ix%100, iy%100);
            if ((1*(ix>=100)-1*(ix<100)) > 0) {
                //         iring = 38 - iring_map_plus [ixiy];
                iring = iring_map_plus [ixiy];
            }
            else {
                //         iring = 38 - iring_map_minus [ixiy];
                iring = iring_map_minus [ixiy];
            }

            //std::cout << " ix,iy,ring = " << ix%100 << "  " << iy%100 << "  " << iring << std::endl;

            iring_map->SetBinContent(ix,iy,iring);

        }
    }


    // TCanvas* ccRingMap = new TCanvas ("ccRingMap","",800,600);
    iring_map->Draw("colz");
    std::vector<int> d_ix_ieta;
    std::vector<int> d_iy_iphi;
    std::vector<int> d_iz;

    std::vector<float> v_ProductG12ped;
    std::vector<float> v_ProductG6ped;
    std::vector<float> v_ProductG1ped;
    std::vector<float> v_ProductG12rms;
    std::vector<float> v_ProductG6rms;
    std::vector<float> v_ProductG1rms;

    for(int ix=1; ix<=histoEE_G12ped->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_G12ped->GetNbinsY(); iy++){

            d_ix_ieta.push_back(ix%100);
            d_iy_iphi.push_back(iy%100);
            d_iz.push_back     (1*(ix>=100)-1*(ix<100));

            if (histoEE_G12ped->GetBinContent(ix,iy) > 0 ) {
                v_ProductG12ped.push_back( histoEE_G12ped->GetBinContent(ix,iy) );
            }
            //                  else {
            //                    v_ProductG12ped.push_back( -99 );
            //                  }
            if (histoEE_G6ped->GetBinContent(ix,iy) > 0 ) {
                v_ProductG6ped.push_back( histoEE_G6ped->GetBinContent(ix,iy) );
            }
            //      else {
            //        v_ProductG6ped.push_back( -99 );
            //      }
            if (histoEE_G1ped->GetBinContent(ix,iy) > 0 ) {
                v_ProductG1ped.push_back( histoEE_G1ped->GetBinContent(ix,iy) );
            }
            //      else {
            //        v_ProductG1ped.push_back( -99 );
            //      }
            if (histoEE_G12rms->GetBinContent(ix,iy) > 0 ) {
                v_ProductG12rms.push_back( histoEE_G12rms->GetBinContent(ix,iy) );
            }
            //      else {
            //        v_ProductG12rms.push_back( -99 );
            //      }
            if (histoEE_G6rms->GetBinContent(ix,iy) > 0 ) {
                v_ProductG6rms.push_back( histoEE_G6rms->GetBinContent(ix,iy) );
            }
            //      else {
            //        v_ProductG6rms.push_back( -99 );
            //      }
            if (histoEE_G1rms->GetBinContent(ix,iy) > 0 ) {
                v_ProductG1rms.push_back( histoEE_G1rms->GetBinContent(ix,iy) );
            }
            //      else {
            //        v_ProductG1rms.push_back( -99 );
            //      }
        }
    }

    //
    //   iring = sqrt ((ix-100)^2 + (iy-100)^2)
    //

    std::vector<float> ringPlus_G12ped;          std::vector<float> ringMinus_G12ped;
    std::vector<float> ringPlus_G6ped;           std::vector<float> ringMinus_G6ped;
    std::vector<float> ringPlus_G1ped;           std::vector<float> ringMinus_G1ped;

    std::vector<float> ringPlus_G12rms;          std::vector<float> ringMinus_G12rms;
    std::vector<float> ringPlus_G6rms;           std::vector<float> ringMinus_G6rms;
    std::vector<float> ringPlus_G1rms;           std::vector<float> ringMinus_G1rms;

    std::vector<int> ringPlus_G12count;          std::vector<int> ringMinus_G12count;
    std::vector<int> ringPlus_G6count;           std::vector<int> ringMinus_G6count;
    std::vector<int> ringPlus_G1count;           std::vector<int> ringMinus_G1count;


    for (int iter = 0; iter < (50-11+3); iter++) {
        ringPlus_G12ped.push_back(0);            ringMinus_G12ped.push_back(0);
        ringPlus_G6ped.push_back(0);             ringMinus_G6ped.push_back(0);
        ringPlus_G1ped.push_back(0);             ringMinus_G1ped.push_back(0);

        ringPlus_G12rms.push_back(0);            ringMinus_G12rms.push_back(0);
        ringPlus_G6rms.push_back(0);             ringMinus_G6rms.push_back(0);
        ringPlus_G1rms.push_back(0);             ringMinus_G1rms.push_back(0);

        ringPlus_G12count.push_back(0);            ringMinus_G12count.push_back(0);
        ringPlus_G6count.push_back(0);             ringMinus_G6count.push_back(0);
        ringPlus_G1count.push_back(0);             ringMinus_G1count.push_back(0);

    }
    std::cout << d_ix_ieta.size() << std::endl;
    for (int iter = 0; iter < (int) d_ix_ieta.size(); iter++) {
        if (d_iz.at(iter) != 0) {

            float dx = d_ix_ieta.at(iter) - 50;
            float dy = d_iy_iphi.at(iter) - 50;

            float ring = sqrt( dx*dx + dy*dy );

            int iring = round(ring) - 11;  //---- 12 [ = (62 - 50 - 1) from the 2D plot] is the first ring

            iring = -99;
            std::pair<int, int> ixiy (d_ix_ieta.at(iter), d_iy_iphi.at(iter));

            if (d_iz.at(iter) > 0) {
                iring = iring_map_plus [ixiy];
            }
            else {
                iring = iring_map_minus [ixiy];
            }

            if (d_iz.at(iter) > 0) {
                ringPlus_G12ped.at(iring) = ringPlus_G12ped.at(iring) + G12ped.at(iter);
                ringPlus_G6ped.at(iring) = ringPlus_G6ped.at(iring) + G6ped.at(iter);
                ringPlus_G1ped.at(iring) = ringPlus_G1ped.at(iring) + G1ped.at(iter);

                ringPlus_G12rms.at(iring) = ringPlus_G12rms.at(iring) + G12rms.at(iter);
                ringPlus_G6rms.at(iring) = ringPlus_G6rms.at(iring) + G6rms.at(iter);
                ringPlus_G1rms.at(iring) = ringPlus_G1rms.at(iring) + G1rms.at(iter);

                ringPlus_G12count.at(iring) =  ringPlus_G12count.at(iring) + 1 ;
                ringPlus_G6count.at(iring) =   ringPlus_G6count. at(iring) + 1 ;
                ringPlus_G1count.at(iring) =   ringPlus_G1count. at(iring) + 1 ;

            }
            else {

                ringMinus_G12ped.at(iring) = ringMinus_G12ped.at(iring) + G12ped.at(iter);
                ringMinus_G6ped.at(iring) = ringMinus_G6ped.at(iring) + G6ped.at(iter);
                ringMinus_G1ped.at(iring) = ringMinus_G1ped.at(iring) + G1ped.at(iter);

                ringMinus_G12rms.at(iring) = ringMinus_G12rms.at(iring) + G12rms.at(iter);
                ringMinus_G6rms.at(iring) = ringMinus_G6rms.at(iring) + G6rms.at(iter);
                ringMinus_G1rms.at(iring) = ringMinus_G1rms.at(iring) + G1rms.at(iter);

                ringMinus_G12count.at(iring) =  ringMinus_G12count.at(iring) + 1 ;
                ringMinus_G6count.at(iring) =   ringMinus_G6count. at(iring) + 1 ;
                ringMinus_G1count.at(iring) =   ringMinus_G1count. at(iring) + 1 ;
            }

        }
    }


    TGraph *gr_EEPlus_G12ped = new TGraph();     TGraph *gr_EEMinus_G12ped = new TGraph();
    TGraph *gr_EEPlus_G6ped  = new TGraph();     TGraph *gr_EEMinus_G6ped  = new TGraph();
    TGraph *gr_EEPlus_G1ped  = new TGraph();     TGraph *gr_EEMinus_G1ped  = new TGraph();

    TGraph *gr_EEPlus_G12rms = new TGraph();     TGraph *gr_EEMinus_G12rms = new TGraph();
    TGraph *gr_EEPlus_G6rms  = new TGraph();     TGraph *gr_EEMinus_G6rms  = new TGraph();
    TGraph *gr_EEPlus_G1rms  = new TGraph();     TGraph *gr_EEMinus_G1rms  = new TGraph();



    for (int iter = 0; iter < (39); iter++) {


        gr_EEPlus_G12ped-> SetPoint (iter, iter,   ringPlus_G12count.at(iter) ? ringPlus_G12ped.at(iter) / ringPlus_G12count.at(iter) : 0 ) ;
        gr_EEPlus_G6ped-> SetPoint (iter,   iter,    ringPlus_G6count.at(iter)  ? ringPlus_G6ped.at(iter)  / ringPlus_G6count.at(iter)  : 0 ) ;
        gr_EEPlus_G1ped-> SetPoint (iter,   iter,    ringPlus_G1count.at(iter)  ? ringPlus_G1ped.at(iter)  / ringPlus_G1count.at(iter)  : 0 ) ;

        gr_EEMinus_G12ped-> SetPoint (iter,   iter,   ringMinus_G12count.at(iter) ?  ringMinus_G12ped.at(iter) / ringMinus_G12count.at(iter) : 0 ) ;
        gr_EEMinus_G6ped-> SetPoint (iter,  iter,    ringMinus_G6count.at(iter)  ?  ringMinus_G6ped.at(iter)  / ringMinus_G6count.at(iter)  : 0 ) ;
        gr_EEMinus_G1ped-> SetPoint (iter,  iter,    ringMinus_G1count.at(iter)  ?  ringMinus_G1ped.at(iter)  / ringMinus_G1count.at(iter)  : 0 ) ;

        gr_EEPlus_G12rms-> SetPoint (iter,  iter,   ringPlus_G12count.at(iter) ?  ringPlus_G12rms.at(iter) / ringPlus_G12count.at(iter) : 0 ) ;
        gr_EEPlus_G6rms-> SetPoint (iter,   iter,    ringPlus_G6count.at(iter)  ?  ringPlus_G6rms.at(iter)  / ringPlus_G6count.at(iter)  : 0 ) ;
        gr_EEPlus_G1rms-> SetPoint (iter,    iter,    ringPlus_G1count.at(iter)  ?  ringPlus_G1rms.at(iter)  / ringPlus_G1count.at(iter)  : 0 ) ;

        gr_EEMinus_G12rms-> SetPoint (iter,  iter,   ringMinus_G12count.at(iter) ?  ringMinus_G12rms.at(iter) / ringMinus_G12count.at(iter) : 0 ) ;
        gr_EEMinus_G6rms-> SetPoint (iter,   iter,    ringMinus_G6count.at(iter)  ?  ringMinus_G6rms.at(iter)  / ringMinus_G6count.at(iter)  : 0 ) ;
        gr_EEMinus_G1rms-> SetPoint (iter,   iter,    ringMinus_G1count.at(iter)  ?  ringMinus_G1rms.at(iter)  / ringMinus_G1count.at(iter)  : 0 ) ;


    }


    //---- style ----

    gr_EEPlus_G12ped->SetMarkerSize  (1);                           gr_EEMinus_G12ped->SetMarkerSize  (1);
    gr_EEPlus_G12ped->SetMarkerStyle (24);                          gr_EEMinus_G12ped->SetMarkerStyle (22);
    gr_EEPlus_G12ped->SetMarkerColor (kRed);                        gr_EEMinus_G12ped->SetMarkerColor (kRed);
    gr_EEPlus_G12ped->SetLineWidth (1);                             gr_EEMinus_G12ped->SetLineWidth (1);
    gr_EEPlus_G12ped->SetLineColor (kRed);                          gr_EEMinus_G12ped->SetLineColor (kRed);

    gr_EEPlus_G6ped->SetMarkerSize  (1);                            gr_EEMinus_G6ped->SetMarkerSize  (1);
    gr_EEPlus_G6ped->SetMarkerStyle (24);                           gr_EEMinus_G6ped->SetMarkerStyle (22);
    gr_EEPlus_G6ped->SetMarkerColor (kBlue);                        gr_EEMinus_G6ped->SetMarkerColor (kBlue);
    gr_EEPlus_G6ped->SetLineWidth (1);                              gr_EEMinus_G6ped->SetLineWidth (1);
    gr_EEPlus_G6ped->SetLineColor (kBlue);                          gr_EEMinus_G6ped->SetLineColor (kBlue);

    gr_EEPlus_G1ped->SetMarkerSize  (1);                            gr_EEMinus_G1ped->SetMarkerSize  (1);
    gr_EEPlus_G1ped->SetMarkerStyle (24);                           gr_EEMinus_G1ped->SetMarkerStyle (22);
    gr_EEPlus_G1ped->SetMarkerColor (kAzure);                       gr_EEMinus_G1ped->SetMarkerColor (kAzure);
    gr_EEPlus_G1ped->SetLineWidth (1);                              gr_EEMinus_G1ped->SetLineWidth (1);
    gr_EEPlus_G1ped->SetLineColor (kAzure);                         gr_EEMinus_G1ped->SetLineColor (kAzure);


    gr_EEPlus_G12rms->SetMarkerSize  (1);                           gr_EEMinus_G12rms->SetMarkerSize  (1);
    gr_EEPlus_G12rms->SetMarkerStyle (24);                          gr_EEMinus_G12rms->SetMarkerStyle (22);
    gr_EEPlus_G12rms->SetMarkerColor (kRed);                        gr_EEMinus_G12rms->SetMarkerColor (kRed);
    gr_EEPlus_G12rms->SetLineWidth (1);                             gr_EEMinus_G12rms->SetLineWidth (1);
    gr_EEPlus_G12rms->SetLineColor (kRed);                          gr_EEMinus_G12rms->SetLineColor (kRed);

    gr_EEPlus_G6rms->SetMarkerSize  (1);                            gr_EEMinus_G6rms->SetMarkerSize  (1);
    gr_EEPlus_G6rms->SetMarkerStyle (24);                           gr_EEMinus_G6rms->SetMarkerStyle (22);
    gr_EEPlus_G6rms->SetMarkerColor (kBlue);                        gr_EEMinus_G6rms->SetMarkerColor (kBlue);
    gr_EEPlus_G6rms->SetLineWidth (1);                              gr_EEMinus_G6rms->SetLineWidth (1);
    gr_EEPlus_G6rms->SetLineColor (kBlue);                          gr_EEMinus_G6rms->SetLineColor (kBlue);

    gr_EEPlus_G1rms->SetMarkerSize  (1);                            gr_EEMinus_G1rms->SetMarkerSize  (1);
    gr_EEPlus_G1rms->SetMarkerStyle (24);                           gr_EEMinus_G1rms->SetMarkerStyle (22);
    gr_EEPlus_G1rms->SetMarkerColor (kAzure);                       gr_EEMinus_G1rms->SetMarkerColor (kAzure);
    gr_EEPlus_G1rms->SetLineWidth (1);                              gr_EEMinus_G1rms->SetLineWidth (1);
    gr_EEPlus_G1rms->SetLineColor (kAzure);                         gr_EEMinus_G1rms->SetLineColor (kAzure);


    //---- style (end) ----


    TCanvas* ccRing = new TCanvas ("ccRing","",800,600);
    ccRing->Divide(2,1);

    ccRing->cd(1);

    gr_EEPlus_G12ped->Draw("APL");
    gr_EEPlus_G6ped->Draw("PL");
    gr_EEPlus_G1ped->Draw("PL");

    gr_EEMinus_G12ped->Draw("PL");
    gr_EEMinus_G6ped->Draw("PL");
    gr_EEMinus_G1ped->Draw("PL");


    ccRing->cd(2);

    gr_EEPlus_G12rms->Draw("APL");
    gr_EEPlus_G6rms->Draw("PL");
    gr_EEPlus_G1rms->Draw("PL");

    gr_EEMinus_G12rms->Draw("PL");
    gr_EEMinus_G6rms->Draw("PL");
    gr_EEMinus_G1rms->Draw("PL");


    gr_EEPlus_G12ped->GetYaxis()->SetTitle("ped ADC");
    gr_EEPlus_G12ped->GetXaxis()->SetTitle("iRing");

    gr_EEPlus_G12rms->GetYaxis()->SetTitle("rms ADC");
    gr_EEPlus_G12rms->GetXaxis()->SetTitle("iRing");

    ccRing->SaveAs("PedestalRMSVsRingAllGains_EE.png");
    ccRing->SaveAs("PedestalRMSVsRingAllGains_EE.root");


    TCanvas* ccRing2 = new TCanvas ("ccRing2","",800,600);
    gr_EEPlus_G12rms->Draw("APL");
    gr_EEMinus_G12rms->Draw("PL");

    TLegend* leg_EE = new TLegend(0.20,0.70,0.50,0.90);
    leg_EE->AddEntry(gr_EEPlus_G12rms, "EE+","pl");
    leg_EE->AddEntry(gr_EEMinus_G12rms,"EE-","pl");
    leg_EE->Draw();


    ccRing2->SaveAs("PedestalRMSVsRingGain12_EE.png");
    ccRing2->SaveAs("PedestalRMSVsRingGain12_EE.root");







































    //
    //   iring = ieta in EB
    //

    std::vector<float> EBringPlus_G12ped;          std::vector<float> EBringMinus_G12ped;
    std::vector<float> EBringPlus_G6ped;           std::vector<float> EBringMinus_G6ped;
    std::vector<float> EBringPlus_G1ped;           std::vector<float> EBringMinus_G1ped;

    std::vector<float> EBringPlus_G12rms;          std::vector<float> EBringMinus_G12rms;
    std::vector<float> EBringPlus_G6rms;           std::vector<float> EBringMinus_G6rms;
    std::vector<float> EBringPlus_G1rms;           std::vector<float> EBringMinus_G1rms;

    std::vector<int> EBringPlus_G12count;          std::vector<int> EBringMinus_G12count;
    std::vector<int> EBringPlus_G6count;           std::vector<int> EBringMinus_G6count;
    std::vector<int> EBringPlus_G1count;           std::vector<int> EBringMinus_G1count;


    for (int iter = 0; iter < 85; iter++) {
        EBringPlus_G12ped.push_back(0);            EBringMinus_G12ped.push_back(0);
        EBringPlus_G6ped.push_back(0);             EBringMinus_G6ped.push_back(0);
        EBringPlus_G1ped.push_back(0);             EBringMinus_G1ped.push_back(0);

        EBringPlus_G12rms.push_back(0);            EBringMinus_G12rms.push_back(0);
        EBringPlus_G6rms.push_back(0);             EBringMinus_G6rms.push_back(0);
        EBringPlus_G1rms.push_back(0);             EBringMinus_G1rms.push_back(0);

        EBringPlus_G12count.push_back(0);            EBringMinus_G12count.push_back(0);
        EBringPlus_G6count.push_back(0);             EBringMinus_G6count.push_back(0);
        EBringPlus_G1count.push_back(0);             EBringMinus_G1count.push_back(0);

    }

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            int iEBring = abs(ix_ieta.at(iter)) - 1 ;


            if (iEBring > 84 || iEBring < 0) std::cout << " what ?!?   iEBring = " << iEBring << std::endl;

            if (ix_ieta.at(iter) > 0) {

                EBringPlus_G12ped.at(iEBring) = EBringPlus_G12ped.at(iEBring) + G12ped.at(iter);
                EBringPlus_G6ped.at(iEBring)  = EBringPlus_G6ped.at(iEBring) + G6ped.at(iter);
                EBringPlus_G1ped.at(iEBring)  = EBringPlus_G1ped.at(iEBring) + G1ped.at(iter);

                EBringPlus_G12rms.at(iEBring) = EBringPlus_G12rms.at(iEBring) + G12rms.at(iter);
                EBringPlus_G6rms.at(iEBring)  = EBringPlus_G6rms.at(iEBring) + G6rms.at(iter);
                EBringPlus_G1rms.at(iEBring)  = EBringPlus_G1rms.at(iEBring) + G1rms.at(iter);

                EBringPlus_G12count.at(iEBring) =  EBringPlus_G12count.at(iEBring) + 1 ;
                EBringPlus_G6count.at(iEBring)  =  EBringPlus_G6count. at(iEBring) + 1 ;
                EBringPlus_G1count.at(iEBring)  =  EBringPlus_G1count. at(iEBring) + 1 ;

            }
            else {

                EBringMinus_G12ped.at(iEBring) = EBringMinus_G12ped.at(iEBring) + G12ped.at(iter);
                EBringMinus_G6ped.at(iEBring)  = EBringMinus_G6ped.at(iEBring) + G6ped.at(iter);
                EBringMinus_G1ped.at(iEBring)  = EBringMinus_G1ped.at(iEBring) + G1ped.at(iter);

                EBringMinus_G12rms.at(iEBring) = EBringMinus_G12rms.at(iEBring) + G12rms.at(iter);
                EBringMinus_G6rms.at(iEBring)  = EBringMinus_G6rms.at(iEBring) + G6rms.at(iter);
                EBringMinus_G1rms.at(iEBring)  = EBringMinus_G1rms.at(iEBring) + G1rms.at(iter);

                EBringMinus_G12count.at(iEBring) =  EBringMinus_G12count.at(iEBring) + 1 ;
                EBringMinus_G6count.at(iEBring)  =  EBringMinus_G6count. at(iEBring) + 1 ;
                EBringMinus_G1count.at(iEBring)  =  EBringMinus_G1count. at(iEBring) + 1 ;
            }

        }
    }


    TGraph *gr_EBPlus_G12ped = new TGraph();     TGraph *gr_EBMinus_G12ped = new TGraph();
    TGraph *gr_EBPlus_G6ped  = new TGraph();     TGraph *gr_EBMinus_G6ped  = new TGraph();
    TGraph *gr_EBPlus_G1ped  = new TGraph();     TGraph *gr_EBMinus_G1ped  = new TGraph();

    TGraph *gr_EBPlus_G12rms = new TGraph();     TGraph *gr_EBMinus_G12rms = new TGraph();
    TGraph *gr_EBPlus_G6rms  = new TGraph();     TGraph *gr_EBMinus_G6rms  = new TGraph();
    TGraph *gr_EBPlus_G1rms  = new TGraph();     TGraph *gr_EBMinus_G1rms  = new TGraph();



    for (int iter = 0; iter < 85; iter++) {


        gr_EBPlus_G12ped-> SetPoint (iter,    iter+1,   EBringPlus_G12count.at(iter) ? EBringPlus_G12ped.at(iter) / EBringPlus_G12count.at(iter) : 0 ) ;
        gr_EBPlus_G6ped-> SetPoint (iter,    iter+1,    EBringPlus_G6count.at(iter)  ? EBringPlus_G6ped.at(iter)  / EBringPlus_G6count.at(iter)  : 0 ) ;
        gr_EBPlus_G1ped-> SetPoint (iter,    iter+1,    EBringPlus_G1count.at(iter)  ? EBringPlus_G1ped.at(iter)  / EBringPlus_G1count.at(iter)  : 0 ) ;

        gr_EBMinus_G12ped-> SetPoint (iter,    iter+1,   EBringMinus_G12count.at(iter) ?  EBringMinus_G12ped.at(iter) / EBringMinus_G12count.at(iter) : 0 ) ;
        gr_EBMinus_G6ped-> SetPoint (iter,    iter+1,    EBringMinus_G6count.at(iter)  ?  EBringMinus_G6ped.at(iter)  / EBringMinus_G6count.at(iter)  : 0 ) ;
        gr_EBMinus_G1ped-> SetPoint (iter,    iter+1,    EBringMinus_G1count.at(iter)  ?  EBringMinus_G1ped.at(iter)  / EBringMinus_G1count.at(iter)  : 0 ) ;

        gr_EBPlus_G12rms-> SetPoint (iter,    iter+1,   EBringPlus_G12count.at(iter) ?  EBringPlus_G12rms.at(iter) / EBringPlus_G12count.at(iter) : 0 ) ;
        gr_EBPlus_G6rms-> SetPoint (iter,    iter+1,    EBringPlus_G6count.at(iter)  ?  EBringPlus_G6rms.at(iter)  / EBringPlus_G6count.at(iter)  : 0 ) ;
        gr_EBPlus_G1rms-> SetPoint (iter,    iter+1,    EBringPlus_G1count.at(iter)  ?  EBringPlus_G1rms.at(iter)  / EBringPlus_G1count.at(iter)  : 0 ) ;

        gr_EBMinus_G12rms-> SetPoint (iter,    iter+1,   EBringMinus_G12count.at(iter) ?  EBringMinus_G12rms.at(iter) / EBringMinus_G12count.at(iter) : 0 ) ;
        gr_EBMinus_G6rms-> SetPoint (iter,    iter+1,    EBringMinus_G6count.at(iter)  ?  EBringMinus_G6rms.at(iter)  / EBringMinus_G6count.at(iter)  : 0 ) ;
        gr_EBMinus_G1rms-> SetPoint (iter,    iter+1,    EBringMinus_G1count.at(iter)  ?  EBringMinus_G1rms.at(iter)  / EBringMinus_G1count.at(iter)  : 0 ) ;


    }


    //---- style ----

    gr_EBPlus_G12ped->SetMarkerSize  (1);                           gr_EBMinus_G12ped->SetMarkerSize  (1);
    gr_EBPlus_G12ped->SetMarkerStyle (24);                          gr_EBMinus_G12ped->SetMarkerStyle (22);
    gr_EBPlus_G12ped->SetMarkerColor (kRed);                        gr_EBMinus_G12ped->SetMarkerColor (kRed);
    gr_EBPlus_G12ped->SetLineWidth (1);                             gr_EBMinus_G12ped->SetLineWidth (1);
    gr_EBPlus_G12ped->SetLineColor (kRed);                          gr_EBMinus_G12ped->SetLineColor (kRed);

    gr_EBPlus_G6ped->SetMarkerSize  (1);                            gr_EBMinus_G6ped->SetMarkerSize  (1);
    gr_EBPlus_G6ped->SetMarkerStyle (24);                           gr_EBMinus_G6ped->SetMarkerStyle (22);
    gr_EBPlus_G6ped->SetMarkerColor (kBlue);                        gr_EBMinus_G6ped->SetMarkerColor (kBlue);
    gr_EBPlus_G6ped->SetLineWidth (1);                              gr_EBMinus_G6ped->SetLineWidth (1);
    gr_EBPlus_G6ped->SetLineColor (kBlue);                          gr_EBMinus_G6ped->SetLineColor (kBlue);

    gr_EBPlus_G1ped->SetMarkerSize  (1);                            gr_EBMinus_G1ped->SetMarkerSize  (1);
    gr_EBPlus_G1ped->SetMarkerStyle (24);                           gr_EBMinus_G1ped->SetMarkerStyle (22);
    gr_EBPlus_G1ped->SetMarkerColor (kAzure);                       gr_EBMinus_G1ped->SetMarkerColor (kAzure);
    gr_EBPlus_G1ped->SetLineWidth (1);                              gr_EBMinus_G1ped->SetLineWidth (1);
    gr_EBPlus_G1ped->SetLineColor (kAzure);                         gr_EBMinus_G1ped->SetLineColor (kAzure);


    gr_EBPlus_G12rms->SetMarkerSize  (1);                           gr_EBMinus_G12rms->SetMarkerSize  (1);
    gr_EBPlus_G12rms->SetMarkerStyle (24);                          gr_EBMinus_G12rms->SetMarkerStyle (22);
    gr_EBPlus_G12rms->SetMarkerColor (kRed);                        gr_EBMinus_G12rms->SetMarkerColor (kRed);
    gr_EBPlus_G12rms->SetLineWidth (1);                             gr_EBMinus_G12rms->SetLineWidth (1);
    gr_EBPlus_G12rms->SetLineColor (kRed);                          gr_EBMinus_G12rms->SetLineColor (kRed);

    gr_EBPlus_G6rms->SetMarkerSize  (1);                            gr_EBMinus_G6rms->SetMarkerSize  (1);
    gr_EBPlus_G6rms->SetMarkerStyle (24);                           gr_EBMinus_G6rms->SetMarkerStyle (22);
    gr_EBPlus_G6rms->SetMarkerColor (kBlue);                        gr_EBMinus_G6rms->SetMarkerColor (kBlue);
    gr_EBPlus_G6rms->SetLineWidth (1);                              gr_EBMinus_G6rms->SetLineWidth (1);
    gr_EBPlus_G6rms->SetLineColor (kBlue);                          gr_EBMinus_G6rms->SetLineColor (kBlue);

    gr_EBPlus_G1rms->SetMarkerSize  (1);                            gr_EBMinus_G1rms->SetMarkerSize  (1);
    gr_EBPlus_G1rms->SetMarkerStyle (24);                           gr_EBMinus_G1rms->SetMarkerStyle (22);
    gr_EBPlus_G1rms->SetMarkerColor (kAzure);                       gr_EBMinus_G1rms->SetMarkerColor (kAzure);
    gr_EBPlus_G1rms->SetLineWidth (1);                              gr_EBMinus_G1rms->SetLineWidth (1);
    gr_EBPlus_G1rms->SetLineColor (kAzure);                         gr_EBMinus_G1rms->SetLineColor (kAzure);


    //---- style (end) ----

    TLegend* leg_EB = new TLegend(0.20,0.70,0.50,0.90);
    leg_EB->AddEntry(gr_EBPlus_G12ped, "EB+","pl");
    leg_EB->AddEntry(gr_EBMinus_G12ped,"EB-","pl");




    TCanvas* ccRingEB = new TCanvas ("ccRingEB","",800,600);
    ccRingEB->Divide(2,1);

    ccRingEB->cd(1);

    gr_EBPlus_G12ped->Draw("APL");
    gr_EBPlus_G6ped->Draw("PL");
    gr_EBPlus_G1ped->Draw("PL");

    gr_EBMinus_G12ped->Draw("PL");
    gr_EBMinus_G6ped->Draw("PL");
    gr_EBMinus_G1ped->Draw("PL");


    ccRingEB->cd(2);

    gr_EBPlus_G12rms->Draw("APL");
    gr_EBPlus_G6rms->Draw("PL");
    gr_EBPlus_G1rms->Draw("PL");

    gr_EBMinus_G12rms->Draw("PL");
    gr_EBMinus_G6rms->Draw("PL");
    gr_EBMinus_G1rms->Draw("PL");

    gr_EBPlus_G12ped->GetYaxis()->SetTitle("ped ADC");
    gr_EBPlus_G12ped->GetXaxis()->SetTitle("i#eta");

    gr_EBPlus_G12rms->GetYaxis()->SetTitle("rms ADC");
    gr_EBPlus_G12rms->GetXaxis()->SetTitle("i#eta");
    ccRingEB->SaveAs("PedestalRMSVsRingAllGains_EB.png");
    ccRingEB->SaveAs("PedestalRMSVsRingAllGains_EB.root");


    TCanvas* ccRingEB2 = new TCanvas ("ccRingEB2","",800,600);
    gr_EBPlus_G12rms->Draw("APL");
    gr_EBMinus_G12rms->Draw("PL");
    ccRingEB2->SaveAs("PedestalRMSVsRingGain12_EB.png");
    ccRingEB2->SaveAs("PedestalRMSVsRingGain12_EB.root");

    leg_EB->Draw();
    //  dum_pedEB = (TH2F*) histoEB_G12rms->Clone();
    //  dum_pedEE = (TH2F*) histoEE_G12rms->Clone();
    TFile myOutfile("myOutfile.root","RECREATE");
    myOutfile.cd();
    histoEB_G12rms->Write();
    histoEE_G12rms->Write();
    gr_EEPlus_G12rms->Write("gr_EEPlus_G12rms");
    gr_EEMinus_G12rms->Write("gr_EEMinus_G12rms");
    gr_EBPlus_G12rms->Write("gr_EBPlus_G12rms");
    gr_EBMinus_G12rms->Write("gr_EBMinus_G12rms");
}


void plotLaser(cond::NoiseLaserDumper laser) {
    //void plotLaser(std::string nameInputFile = "/afs/cern.ch/work/c/crovelli/dpgTests/dbDumper_902/src/usercode/DBDump/bin/Run304292__6472469664432652288") {

    gStyle->SetOptStat(0);

    TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);
    TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);

    // input file format:
    // rawId     p1    p2   p3   ix/ieta iy/iphi   iz/0

    std::vector<float> p1;
    std::vector<float> p2;
    std::vector<float> p3;

    std::vector<int> ix_ieta;
    std::vector<int> iy_iphi;
    std::vector<int> iz;

    std::vector<int> rawId;

    cond::NoiseLaserDumper::NoiseStruct strct_laser;
    for(int i =0 ; i < (int) laser.NoiseVariables.size(); ++i)
    {
        strct_laser = laser.NoiseVariables.at(i);

        std::cout << strct_laser.IOVBegin  << " IOV "<< strct_laser.IOVEnd  << " apdpn size " <<strct_laser.noise_apdpn.size() << std::endl;
        for(int i =0 ; i < (int) strct_laser.noise_apdpn.size(); ++i)
        {
            ix_ieta.push_back(strct_laser.noise_apdpn.at(i).ix);
            iy_iphi.push_back(strct_laser.noise_apdpn.at(i).iy);
            iz.push_back(strct_laser.noise_apdpn.at(i).iz);

            rawId.push_back(strct_laser.noise_apdpn.at(i).id);

            p1.push_back(strct_laser.noise_apdpn.at(i).p1);
            p2.push_back(strct_laser.noise_apdpn.at(i).p2);
            p3.push_back(strct_laser.noise_apdpn.at(i).p3);

            //   std::cout << "p2: "<< strct_laser.noise_apdpn.at(i).p2 << std::endl;

        }

    }
    // checks
    std::cout << " ix_ieta.size() = " << ix_ieta.size() << std::endl;
    if (ix_ieta.size() > 75848) {
        std::cout << " Attention: you appended the tag twice or you are running on more IOVs!" << std::endl;
    }


    //---- EB ----
    histoEB_p1->GetXaxis()->SetTitle("i#phi");
    histoEB_p1->GetYaxis()->SetTitle("i#eta");

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {
            histoEB_p1->Fill(iy_iphi.at(iter), ix_ieta.at(iter), p1.at(iter) );
            // cout << "eta = " << ix_ieta.at(iter) << ", phi = " << iy_iphi.at(iter) << ", z = " << iz.at(iter) << ", p1 = " << p1.at(iter) << endl;
        }
    }

    ccEB->cd();
    histoEB_p1->Draw("colz");
    ccEB->SaveAs("laserEB.root");
    ccEB->SaveAs("laserEB.png");


    //---- EE ----
    histoEE_p1->GetXaxis()->SetTitle("x");
    histoEE_p1->GetYaxis()->SetTitle("y");

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) != 0) {
            histoEE_p1->Fill(ix_ieta.at(iter) + 100*(iz.at(iter)>0), iy_iphi.at(iter), p1.at(iter) );
            //cout << "x = " << ix_ieta.at(iter) << ", y = " << iy_iphi.at(iter) << ", z = " << iz.at(iter) << ", p1 = " << p1.at(iter) << endl;
        }
    }

    ccEE->cd();
    histoEE_p1->Draw("colz");
    ccEE->SaveAs("laserEE.root");
    ccEE->SaveAs("laserEE.png");


    // ---- read iring definition from file
    //        (ix, iy, iz) -> ring
    //


    std::map < std::pair<int, int> , int > iring_map_plus;
    std::map < std::pair<int, int> , int > iring_map_minus;

    std::ifstream fileEEring ("../usercode/data/eerings.dat");

    std::string buffer;
//    int num;
    if (!fileEEring.is_open()) {
        std::cerr << "** ERROR: Can't open for input" << std::endl;
        //return false;
    }

    while(!fileEEring.eof()) {
        getline(fileEEring,buffer);
        if (buffer != ""){ ///---> save from empty line at the end!
            int ix;
            int iy;
            int iz;
            int iring;

            //       std::cout << " buffer = " << buffer << std::endl;

            std::stringstream line( buffer );
            line >> ix;
            line >> iy;
            line >> iz;
            line >> iring;

            std::pair<int, int> ixiy (ix, iy);
            //       if (iz>0) iring_map_plus  [ixiy] = 38 - iring;
            //       else      iring_map_minus [ixiy] = 38 - iring;

            if (iz>0) iring_map_plus  [ixiy] = iring;
            else      iring_map_minus [ixiy] = iring;

        }
    }


    TH2F *iring_map = (TH2F*) histoEE_p1->Clone("iring_map");

    for(int ix=1; ix<=histoEE_p1->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_p1->GetNbinsY(); iy++){
            int iring = -99;
            std::pair<int, int> ixiy (ix%100, iy%100);
            if ((1*(ix>=100)-1*(ix<100)) > 0) {
                //         iring = 38 - iring_map_plus [ixiy];
                iring = iring_map_plus [ixiy];
            }
            else {
                //         iring = 38 - iring_map_minus [ixiy];
                iring = iring_map_minus [ixiy];
            }

           // std::cout << " ix,iy,ring = " << ix%100 << "  " << iy%100 << "  " << iring << std::endl;

            iring_map->SetBinContent(ix,iy,iring);

        }
    }


    // TCanvas* ccRingMap = new TCanvas ("ccRingMap","",800,600);
    iring_map->Draw("colz");

    std::vector<int> d_ix_ieta;
    std::vector<int> d_iy_iphi;
    std::vector<int> d_iz;

    std::vector<float> v_Product;



    std::vector<float> ringPlus_max_Product;      std::vector<float> ringMinus_max_Product;
    std::vector<float> ringPlus_min_Product;      std::vector<float> ringMinus_min_Product;

    std::vector<float> ringPlus_Product;          std::vector<float> ringMinus_Product;
    std::vector<float> ringPlus_ProductSpread;    std::vector<float> ringMinus_ProductSpread;
    std::vector<int> ringPlus_count;              std::vector<int> ringMinus_count;


    for (int iter = 0; iter < (50-11+3); iter++) {
        ringPlus_max_Product.push_back(0);        ringMinus_max_Product.push_back(0);
        ringPlus_min_Product.push_back(99);       ringMinus_min_Product.push_back(99);
        ringPlus_Product.push_back(0);            ringMinus_Product.push_back(0);
        ringPlus_ProductSpread.push_back(0);      ringMinus_ProductSpread.push_back(0);
        ringPlus_count.push_back(0);              ringMinus_count.push_back(0);
    }


    for(int ix=1; ix<=histoEE_p1->GetNbinsX(); ix++){
        for(int iy=1; iy<=histoEE_p1->GetNbinsY(); iy++){

            d_ix_ieta.push_back(ix%100);
            d_iy_iphi.push_back(iy%100);
            d_iz.push_back     (1*(ix>=100)-1*(ix<100));

            if (histoEE_p1->GetBinContent(ix,iy) > 0 ) {
                v_Product.push_back( histoEE_p1->GetBinContent(ix,iy) );
            }
            else {
                v_Product.push_back( -99 );
            }
        }
    }

    std::cout << d_ix_ieta.size() << std::endl;
    for (int iter = 0; iter < (int) d_ix_ieta.size(); iter++) {
        if (d_iz.at(iter) != 0) {

            float dx = d_ix_ieta.at(iter) - 50;
            float dy = d_iy_iphi.at(iter) - 50;

            float ring = sqrt( dx*dx + dy*dy );

            int iring = round(ring) - 11;  //---- 12 [ = (62 - 50 - 1) from the 2D plot] is the first ring

            iring = -99;
            std::pair<int, int> ixiy (d_ix_ieta.at(iter), d_iy_iphi.at(iter));

            //       std::map<std::pair<int, int>,int>::const_iterator it = iring_map_plus.find(ixiy);
            //       if (it==iring_map_plus.end() ) {
            //         continue;
            //       }
            //

            if (d_iz.at(iter) > 0) {
                iring = iring_map_plus [ixiy];
            }
            else {
                iring = iring_map_minus [ixiy];
            }

            //       std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;

            //       if (iring < 0 ) continue;

            //       if (iring > 38 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;
            //       if (iring > 37 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << std::endl;

            if (v_Product.at(iter) > 0 ) if (iring == 0 )  std::cout << " ix, iy, iring = " << d_ix_ieta.at(iter) << "  " <<  d_iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << " tot = " << ringPlus_count.at(iring) << std::endl;


            if (v_Product.at(iter) > 0 ) {
                if (iring > (50-11+2) || iring < 0) std::cout << " what ?!?   iring = " << iring << " dx = " << dx << " dy = " << dy << " :::: ix = " << d_ix_ieta.at(iter) << "  iy = " << d_iy_iphi.at(iter) << " prod = " << v_Product.at(iter) << std::endl;

                if (d_iz.at(iter) > 0) {

                    if (ringPlus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringPlus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringPlus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringPlus_min_Product.at(iring)   =  v_Product.at(iter);
                    }

                    ringPlus_Product.at(iring)       = ringPlus_Product.at(iring) + v_Product.at(iter);
                    ringPlus_ProductSpread.at(iring) = ringPlus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringPlus_count.at(iring)         = ringPlus_count.at(iring) + 1 ;

                }
                else {

                    if (ringMinus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringMinus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringMinus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringMinus_min_Product.at(iring)   =  v_Product.at(iter);
                    }

                    ringMinus_Product.at(iring)       = ringMinus_Product.at(iring) + v_Product.at(iter);
                    ringMinus_ProductSpread.at(iring) = ringMinus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringMinus_count.at(iring)         = ringMinus_count.at(iring) + 1 ;
                }
            }

        }
    }
  //  std::cout << "Is it here now?!??!" << std::endl;

    for (int iring = 0; iring < (50-11+3); iring++) {
        //     std::cout << " ringPlus_ProductSpread.at(" << iring << ") = " << ringPlus_ProductSpread.at(iring) << " ---> " << ringPlus_ProductSpread.at(iring)  << " - " <<  ringPlus_Product.at(iring)  * ringPlus_Product.at(iring) << std::endl;
        ringMinus_ProductSpread.at(iring) = sqrt(ringMinus_ProductSpread.at(iring) -  ringMinus_Product.at(iring) / ringMinus_count.at(iring) * ringMinus_Product.at(iring)) / ringMinus_count.at(iring);
        ringPlus_ProductSpread.at(iring)  = sqrt(ringPlus_ProductSpread.at(iring)  -  ringPlus_Product.at(iring)  / ringPlus_count.at(iring)  * ringPlus_Product.at(iring))  / ringPlus_count.at(iring);
    }



    //  // Average over rings
    //  std::vector<float> ringPlus_p1;   std::vector<float> ringMinus_p1;
    //  std::vector<int> ringPlus_count;  std::vector<int> ringMinus_count;

    //  for (int iter = 0; iter < (50-11+3); iter++) {
    //    ringPlus_p1.push_back(0);            ringMinus_p1.push_back(0);
    //    ringPlus_count.push_back(0);         ringMinus_count.push_back(0);
    //  }

    //  for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
    //    if (iz.at(iter) != 0) {

    //      float dx = ix_ieta.at(iter) - 50;
    //      float dy = iy_iphi.at(iter) - 50;

    //      float ring = sqrt( dx*dx + dy*dy );

    //      int iring = round(ring) - 11;  //---- 12 [ = (62 - 50 - 1) from the 2D plot] is the first ring

    //      if (iring > (50-11+2) || iring < 0) std::cout << " what ?!?   iring = " << iring << " dx = " << dx << " dy = " << dy << " :::: ix = " << ix_ieta.at(iter) << "  iy = " << iy_iphi.at(iter) << std::endl;

    //      if (iz.at(iter) > 0) {
    //        ringPlus_p1.at(iring)    = ringPlus_p1.at(iring) + p1.at(iter);
    //        ringPlus_count.at(iring) = ringPlus_count.at(iring) + 1 ;
    //      }
    //      else {
    //        ringMinus_p1.at(iring)    = ringMinus_p1.at(iring) + p1.at(iter);
    //        ringMinus_count.at(iring) = ringMinus_count.at(iring) + 1 ;
    //      }
    //    }
    //  }

    TGraph *gr_EEPlus_p1 = new TGraph();
    TGraph *gr_EEMinus_p1 = new TGraph();
    for (int iter = 0; iter < (39); iter++) {
        //std::cout << "Hi! ringPlus_Product.at(" << iter << ") = " << ringPlus_Product.at(iter) << " N = " << ringPlus_count.at(iter) <<  " --> " << ringPlus_Product.at(iter) / ringPlus_count.at(iter) << std::endl;
        //gr_EEPlus_Product-> SetPoint (iter,  iter,   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_EEPlus_p1->  SetPoint (iter,  iter,  ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_EEMinus_p1-> SetPoint (iter, iter,  ringMinus_count.at(iter) ? ringPlus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;
    }


    //---- style ----
    gr_EEPlus_p1->SetMarkerSize  (1);                            gr_EEMinus_p1->SetMarkerSize  (1);
    gr_EEPlus_p1->SetMarkerStyle (24);                           gr_EEMinus_p1->SetMarkerStyle (22);
    gr_EEPlus_p1->SetMarkerColor (kBlue);                        gr_EEMinus_p1->SetMarkerColor (kBlue);
    gr_EEPlus_p1->SetLineWidth (1);                              gr_EEMinus_p1->SetLineWidth (1);
    gr_EEPlus_p1->SetLineColor (kBlue);                          gr_EEMinus_p1->SetLineColor (kBlue);

    // plots
    TCanvas* ccRing = new TCanvas ("ccRing","",800,600);
    ccRing->cd();
    gr_EEPlus_p1->Draw("APL");
    gr_EEMinus_p1->Draw("PL");
    gr_EEPlus_p1->GetYaxis()->SetTitle("Relative response to laser");
    gr_EEPlus_p1->GetXaxis()->SetTitle("iRing");
    ccRing->SaveAs("laserRingEE.root");
    ccRing->SaveAs("laserRingEE.png");


    //
    //   iring = ieta in EB
    //
    std::vector<float> EBringPlus_p1;      std::vector<float> EBringMinus_p1;
    std::vector<int> EBringPlus_count;     std::vector<int> EBringMinus_count;
    for (int iter = 0; iter < 85; iter++) {
        EBringPlus_p1.push_back(0);         EBringMinus_p1.push_back(0);
        EBringPlus_count.push_back(0);      EBringMinus_count.push_back(0);
    }

    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) == 0) {

            int iEBring = abs(ix_ieta.at(iter)) - 1 ;

            if (iEBring > 84 || iEBring < 0) std::cout << " what ?!?   iEBring = " << iEBring << std::endl;

            if (ix_ieta.at(iter) > 0) {
                EBringPlus_p1.at(iEBring)    = EBringPlus_p1.at(iEBring) + p1.at(iter);
                EBringPlus_count.at(iEBring) = EBringPlus_count.at(iEBring) + 1 ;
            }
            else {
                EBringMinus_p1.at(iEBring)    = EBringMinus_p1.at(iEBring) + p1.at(iter);
                EBringMinus_count.at(iEBring) = EBringMinus_count.at(iEBring) + 1 ;
            }
        }
    }

    TGraph *gr_EBPlus_p1 = new TGraph();     TGraph *gr_EBMinus_p1 = new TGraph();
    for (int iter = 0; iter < 85; iter++) {
        gr_EBPlus_p1 -> SetPoint (iter,    iter+1,   EBringPlus_count.at(iter) ? EBringPlus_p1.at(iter) / EBringPlus_count.at(iter) : 0 ) ;
        gr_EBMinus_p1-> SetPoint (iter,    iter+1,   EBringMinus_count.at(iter) ?  EBringMinus_p1.at(iter) / EBringMinus_count.at(iter) : 0 ) ;
    }


    //---- style ----
    gr_EBPlus_p1->SetMarkerSize  (1);                           gr_EBMinus_p1->SetMarkerSize  (1);
    gr_EBPlus_p1->SetMarkerStyle (24);                          gr_EBMinus_p1->SetMarkerStyle (22);
    gr_EBPlus_p1->SetMarkerColor (kRed);                        gr_EBMinus_p1->SetMarkerColor (kRed);
    gr_EBPlus_p1->SetLineWidth (1);                             gr_EBMinus_p1->SetLineWidth (1);
    gr_EBPlus_p1->SetLineColor (kRed);                          gr_EBMinus_p1->SetLineColor (kRed);

    // plots
    TCanvas* ccRingEB = new TCanvas ("ccRingEB","",800,600);
    ccRingEB->Divide(2,1);
    gr_EBPlus_p1->Draw("APL");
    gr_EBMinus_p1->Draw("PL");
    gr_EBPlus_p1->GetYaxis()->SetTitle("Relative response to laser");
    gr_EBPlus_p1->GetXaxis()->SetTitle("i#eta");
    ccRingEB->SaveAs("laserRingEB.root");
    ccRingEB->SaveAs("laserRingEB.png");
    //  dum_laserEB = (TH2F*) histoEB_p1->Clone();
    //  dum_laserEE = (TH2F*) histoEE_p1->Clone();
    TFile outFile("fileOutLaser.root","RECREATE");
    outFile.cd();
    histoEE_p1->Write();
    histoEB_p1->Write();
    gr_EBPlus_p1->Write("gr_EBPlus_p1");
    gr_EBMinus_p1->Write("gr_EBMinus_p1");
    gr_EEPlus_p1->Write("gr_EEPlus_p1");
    gr_EEMinus_p1->Write("gr_EEMinus_p1");
}

void productOfTag(cond::NoiseDumper<EcalADCToGeVConstant> ADCtoGEV)
{
    float adcToGeVEB = ADCtoGEV.a.adcToGeVEB;
    float adcToGeVEE = ADCtoGEV.a.adcToGeVEE;
    TH2F *pedEB = (TH2F*) histoEB_G12rms->Clone();
    TH2F *pedEE = (TH2F*) histoEE_G12rms->Clone();

    TH2F *icEB = (TH2F*) histoEB_IC->Clone();
    TH2F *icEE = (TH2F*) histoEE_IC->Clone();

    TH2F *alphaEB = (TH2F*) histoEB_alpha->Clone();
    TH2F *alphaEE = (TH2F*) histoEE_alpha->Clone();

    TH2F *laserEB = (TH2F*) histoEB_p1->Clone();
    TH2F *laserEE = (TH2F*) histoEE_p1->Clone();
    std::cout << adcToGeVEB << " ADC " << adcToGeVEE << std::endl;

    //---- output
    TFile fileOut("product.root","RECREATE");
    fileOut.cd();

    // Product
    TH2F *productEB = (TH2F*) icEB->Clone("productEB");
    TH2F *productEBChinese = (TH2F*) icEB->Clone("productEBChinese");
    TH2F *productEBRussian = (TH2F*) icEB->Clone("productEBRussian");

    TH2F *productEE = (TH2F*) icEE->Clone("productEE");
    TH2F *productEEChinese = (TH2F*) icEE->Clone("productEEChinese");
    TH2F *productEERussian = (TH2F*) icEE->Clone("productEERussian");




    TGraphErrors *gr_ring_EBPlus_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EBMinus_eta_Product = new TGraphErrors();
    TGraphErrors *gr_ring_EBPlus_max_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EBMinus_max_eta_Product = new TGraphErrors();
    TGraphErrors *gr_ring_EBPlus_min_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EBMinus_min_eta_Product = new TGraphErrors();



    std::vector<float> ringPlus_EB_max_Product;      std::vector<float> ringMinus_EB_max_Product;
    std::vector<float> ringPlus_EB_min_Product;      std::vector<float> ringMinus_EB_min_Product;

    std::vector<float> ringPlus_EB_Product;          std::vector<float> ringMinus_EB_Product;
    std::vector<int> ringPlus_EB_count;              std::vector<int> ringMinus_EB_count;


    for (int iter = 0; iter < (85); iter++) {
        ringPlus_EB_max_Product.push_back(0);        ringMinus_EB_max_Product.push_back(0);
        ringPlus_EB_min_Product.push_back(99);       ringMinus_EB_min_Product.push_back(99);
        ringPlus_EB_Product.push_back(0);            ringMinus_EB_Product.push_back(0);
        ringPlus_EB_count.push_back(0);              ringMinus_EB_count.push_back(0);
    }



    for(int ix=1; ix<=productEB->GetNbinsX(); ix++){ //---- phi
        for(int iy=1; iy<=productEB->GetNbinsY(); iy++){ //---- eta

            float VpedestalEB = pedEB     -> GetBinContent(ix,iy);
            float ValphaEB    = alphaEB   -> GetBinContent(ix,iy);
            float VlaserEB    = laserEB   -> GetBinContent(ix,iy);
            float VicEB       = productEB -> GetBinContent(ix,iy);
            float myProdEB = VpedestalEB*adcToGeVEB*pow((1./VlaserEB),ValphaEB)*VicEB;

            //       std::cout << VpedestalEB << " " << adcToGeVEB << " " << VlaserEB << " " << ValphaEB << " " << VicEB << " ==> " << myProdEB << std::endl;

            productEB->SetBinContent(ix,iy,myProdEB);
            if((ValphaEB - 1.52) < 0.0001 )
            {
                productEBRussian->SetBinContent(ix,iy,myProdEB);
                productEBChinese->SetBinContent(ix,iy,-0.05);

            }
            if((ValphaEB - 1.50) < 0.0001)
            {


                productEBChinese->SetBinContent(ix,iy,myProdEB);
                productEBRussian->SetBinContent(ix,iy,-0.05);

            }
//            else
//            {
//                std::cout << ValphaEB << std::endl;

//                productEBRussian->SetBinContent(ix,iy,0);
//                productEBChinese->SetBinContent(ix,iy,0);

//            }



            //       std::cout << " iy = " << iy << std::endl;

            if (iy<=85) {
                //---- EB-
                ringMinus_EB_count.at(84 - (iy -1) ) += 1;
                ringMinus_EB_Product.at(84 - (iy -1) ) += myProdEB;

                if (ringMinus_EB_max_Product.at(84 - (iy -1) ) < myProdEB ) {
                    ringMinus_EB_max_Product.at(84 - (iy -1) )   =  myProdEB;
                }
                if (ringMinus_EB_min_Product.at(84 - (iy -1) ) > myProdEB ) {
                    ringMinus_EB_min_Product.at(84 - (iy -1) )   =  myProdEB;
                }
            }
            else if (iy>86) {
                //---- EB+
                ringPlus_EB_count.at(iy - 85 - 2) += 1;
                ringPlus_EB_Product.at(iy - 85 - 2) += myProdEB;

                if (ringPlus_EB_max_Product.at(iy - 85 - 2) < myProdEB ) {
                    ringPlus_EB_max_Product.at(iy - 85 - 2)   =  myProdEB;
                }
                if (ringPlus_EB_min_Product.at(iy - 85 - 2) > myProdEB ) {
                    ringPlus_EB_min_Product.at(iy - 85 - 2)   =  myProdEB;
                }

            }

        }
    }


    for(int ix=1; ix<=productEE->GetNbinsX(); ix++){
        for(int iy=1; iy<=productEE->GetNbinsY(); iy++){
            float VpedestalEE = pedEE->GetBinContent(ix,iy);
            float ValphaEE    = alphaEE->GetBinContent(ix,iy);
            float VlaserEE    = laserEE->GetBinContent(ix,iy);
            float VicEE       = productEE->GetBinContent(ix,iy);
            float myProdEE = VpedestalEE*adcToGeVEE*pow((1./VlaserEE),ValphaEE)*VicEE;

            productEE->SetBinContent(ix,iy,myProdEE);
            //std::cout << ValphaEE << std::endl;

            if((ValphaEE - 1.16) < 0.0001)
            {
                productEERussian->SetBinContent(ix,iy,myProdEE);
                productEEChinese->SetBinContent(ix,iy,-0.05);
            }
            if((ValphaEE - 1.0) < 0.0001 )
            {
                productEERussian->SetBinContent(ix,iy,-0.05);

                productEEChinese->SetBinContent(ix,iy,myProdEE);
            }
//            else
//            {
//                std::cout << ValphaEE << std::endl;
//                productEERussian->SetBinContent(ix,iy,0);
//                productEEChinese->SetBinContent(ix,iy,0);

//            }
        }
    }



    TCanvas* ccEB = new TCanvas ("ccEB","",1600,600);
    productEB->Draw("colz");

    TCanvas* ccEE = new TCanvas ("ccEE","",1600,600);
    productEE->Draw("colz");



    productEB->Write("resultEB");
    //productEB->SaveAs("resultEB.png");
    productEE->Write("resultEE");
    //productEE->SaveAs("resultEE.png");
    productEBChinese->Write("resultEBChinese");
    productEEChinese->Write("resultEEChinese");
    productEBRussian->Write("resultEBRussian");
    productEERussian->Write("resultEERussian");
    ccEB->Write();
    ccEE->Write();







    for (int iter = 0; iter < (85); iter++) {


        gr_ring_EBPlus_eta_Product-> SetPoint (iter,  iter,   ringPlus_EB_count.at(iter) ? ringPlus_EB_Product.at(iter) / ringPlus_EB_count.at(iter) : 0 ) ;
        gr_ring_EBMinus_eta_Product-> SetPoint (iter,  iter,   ringMinus_EB_count.at(iter) ? ringMinus_EB_Product.at(iter) / ringMinus_EB_count.at(iter) : 0 ) ;

        gr_ring_EBPlus_max_eta_Product -> SetPoint (iter, iter,   ringPlus_EB_max_Product.at(iter)) ;
        gr_ring_EBPlus_min_eta_Product -> SetPoint (iter, iter,   ringPlus_EB_min_Product.at(iter)) ;

        gr_ring_EBMinus_max_eta_Product -> SetPoint (iter, iter,   ringMinus_EB_max_Product.at(iter)) ;
        gr_ring_EBMinus_min_eta_Product -> SetPoint (iter, iter,   ringMinus_EB_min_Product.at(iter)) ;

    }

    gr_ring_EBPlus_eta_Product->SetMarkerSize  (1);                           gr_ring_EBMinus_eta_Product->SetMarkerSize  (1);
    gr_ring_EBPlus_eta_Product->SetMarkerStyle (24);                          gr_ring_EBMinus_eta_Product->SetMarkerStyle (22);
    gr_ring_EBPlus_eta_Product->SetMarkerColor (kRed);                        gr_ring_EBMinus_eta_Product->SetMarkerColor (kRed);
    gr_ring_EBPlus_eta_Product->SetLineWidth (1);                             gr_ring_EBMinus_eta_Product->SetLineWidth (1);
    gr_ring_EBPlus_eta_Product->SetLineColor (kRed);                          gr_ring_EBMinus_eta_Product->SetLineColor (kRed);

    gr_ring_EBPlus_max_eta_Product->SetMarkerSize  (1);                           gr_ring_EBMinus_max_eta_Product->SetMarkerSize  (1);
    gr_ring_EBPlus_max_eta_Product->SetMarkerStyle (24);                          gr_ring_EBMinus_max_eta_Product->SetMarkerStyle (22);
    gr_ring_EBPlus_max_eta_Product->SetMarkerColor (kBlue);                       gr_ring_EBMinus_max_eta_Product->SetMarkerColor (kBlue);
    gr_ring_EBPlus_max_eta_Product->SetLineWidth (1);                             gr_ring_EBMinus_max_eta_Product->SetLineWidth (1);
    gr_ring_EBPlus_max_eta_Product->SetLineColor (kBlue);                         gr_ring_EBMinus_max_eta_Product->SetLineColor (kBlue);

    gr_ring_EBPlus_min_eta_Product->SetMarkerSize  (1);                           gr_ring_EBMinus_min_eta_Product->SetMarkerSize  (1);
    gr_ring_EBPlus_min_eta_Product->SetMarkerStyle (24);                          gr_ring_EBMinus_min_eta_Product->SetMarkerStyle (22);
    gr_ring_EBPlus_min_eta_Product->SetMarkerColor (kBlue);                       gr_ring_EBMinus_min_eta_Product->SetMarkerColor (kBlue);
    gr_ring_EBPlus_min_eta_Product->SetLineWidth (1);                             gr_ring_EBMinus_min_eta_Product->SetLineWidth (1);
    gr_ring_EBPlus_min_eta_Product->SetLineColor (kBlue);                         gr_ring_EBMinus_min_eta_Product->SetLineColor (kBlue);


    TCanvas* ccRing_EB_ring_minmax = new TCanvas ("ccRing_EB_ring_minmax","",800,600);

    TLegend* legend3 = new TLegend(0.1,0.7,0.48,0.9);
    legend3->AddEntry(gr_ring_EBMinus_eta_Product,"average EB-","lep");
    legend3->AddEntry(gr_ring_EBPlus_eta_Product, "average EB+","lep");
    legend3->AddEntry(gr_ring_EBMinus_min_eta_Product, "min/max EB-","lep");
    legend3->AddEntry(gr_ring_EBPlus_min_eta_Product, "min/max EB+","lep");


    gr_ring_EBPlus_max_eta_Product->Draw("AP");
    gr_ring_EBPlus_min_eta_Product->Draw("P");

    gr_ring_EBPlus_eta_Product->Draw("PL");
    gr_ring_EBMinus_eta_Product->Draw("PL");

    gr_ring_EBMinus_eta_Product->Draw("P");
    gr_ring_EBMinus_max_eta_Product->Draw("P");
    gr_ring_EBMinus_min_eta_Product->Draw("P");

    gr_ring_EBPlus_max_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
    gr_ring_EBPlus_max_eta_Product->GetXaxis()->SetTitle("i#eta");

    legend3->Draw();

    ccRing_EB_ring_minmax->SaveAs("NoiseRing_EB_ring_minmax_ring_EE.png");
    ccRing_EB_ring_minmax->SaveAs("NoiseRing_EB_ring_minmax_ring_EE.root");






    //---- EE ----

    //
    //     iring = sqrt ((ix-100)^2 + (iy-100)^2)
    //


    // ---- read iring definition from file
    //        (ix, iy, iz) -> ring
    //


    std::map < std::pair<int, int> , int > iring_map_plus;
    std::map < std::pair<int, int> , int > iring_map_minus;

    std::ifstream fileEEring ("../usercode/data/eerings.dat");

    std::string buffer;
//    int num;
    if (!fileEEring.is_open()) {
        std::cerr << "** ERROR: Can't open for input" << std::endl;
        //return false;
    }

    while(!fileEEring.eof()) {
        getline(fileEEring,buffer);
        if (buffer != ""){ ///---> save from empty line at the end!
            int ix;
            int iy;
            int iz;
            int iring;

            //       std::cout << " buffer = " << buffer << std::endl;

            std::stringstream line( buffer );
            line >> ix;
            line >> iy;
            line >> iz;
            line >> iring;

            std::pair<int, int> ixiy (ix, iy);
            //       if (iz>0) iring_map_plus  [ixiy] = 38 - iring;
            //       else      iring_map_minus [ixiy] = 38 - iring;

            if (iz>0) iring_map_plus  [ixiy] = iring;
            else      iring_map_minus [ixiy] = iring;

        }
    }


    TH2F *iring_map = (TH2F*) icEE->Clone("iring_map");

    for(int ix=1; ix<=productEE->GetNbinsX(); ix++){
        for(int iy=1; iy<=productEE->GetNbinsY(); iy++){
            int iring = -99;
            std::pair<int, int> ixiy (ix%100, iy%100);
            if ((1*(ix>=100)-1*(ix<100)) > 0) {
                //         iring = 38 - iring_map_plus [ixiy];
                iring = iring_map_plus [ixiy];
            }
            else {
                //         iring = 38 - iring_map_minus [ixiy];
                iring = iring_map_minus [ixiy];
            }

           // std::cout << " ix,iy,ring = " << ix%100 << "  " << iy%100 << "  " << iring << std::endl;

            iring_map->SetBinContent(ix,iy,iring);

        }
    }


//    TCanvas* ccRingMap = new TCanvas ("ccRingMap","",800,600);
    iring_map->Draw("colz");



    std::vector<int> ix_ieta;
    std::vector<int> iy_iphi;
    std::vector<int> iz;

    std::vector<float> v_Product;
    std::vector<float> v_ProductChinese;
    std::vector<float> v_ProductRussian;



    std::vector<float> ringPlus_max_Product;      std::vector<float> ringMinus_max_Product;
    std::vector<float> ringPlus_min_Product;      std::vector<float> ringMinus_min_Product;

    std::vector<float> ringPlus_Product;          std::vector<float> ringMinus_Product;
    std::vector<float> ringPlus_ProductSpread;    std::vector<float> ringMinus_ProductSpread;
    std::vector<int> ringPlus_count;              std::vector<int> ringMinus_count;

    std::vector<float> ringPlus_max_ProductChinese;      std::vector<float> ringMinus_max_ProductChinese;
    std::vector<float> ringPlus_min_ProductChinese;      std::vector<float> ringMinus_min_ProductChinese;

    std::vector<float> ringPlus_ProductChinese;          std::vector<float> ringMinus_ProductChinese;
    std::vector<float> ringPlus_ProductChineseSpread;    std::vector<float> ringMinus_ProductChineseSpread;
    std::vector<int> ringPlus_countChinese;              std::vector<int> ringMinus_countChinese;

    std::vector<float> ringPlus_max_ProductRussian;      std::vector<float> ringMinus_max_ProductRussian;
    std::vector<float> ringPlus_min_ProductRussian;      std::vector<float> ringMinus_min_ProductRussian;

    std::vector<float> ringPlus_ProductRussian;          std::vector<float> ringMinus_ProductRussian;
    std::vector<float> ringPlus_ProductRussianSpread;    std::vector<float> ringMinus_ProductRussianSpread;
    std::vector<int> ringPlus_countRussian;              std::vector<int> ringMinus_countRussian;

    for (int iter = 0; iter < (50-11+3); iter++) {
        ringPlus_max_Product.push_back(0);        ringMinus_max_Product.push_back(0);
        ringPlus_min_Product.push_back(99);       ringMinus_min_Product.push_back(99);
        ringPlus_Product.push_back(0);            ringMinus_Product.push_back(0);
        ringPlus_ProductSpread.push_back(0);      ringMinus_ProductSpread.push_back(0);
        ringPlus_count.push_back(0);              ringMinus_count.push_back(0);

        ringPlus_max_ProductChinese.push_back(0);        ringMinus_max_ProductChinese.push_back(0);
        ringPlus_min_ProductChinese.push_back(99);       ringMinus_min_ProductChinese.push_back(99);
        ringPlus_ProductChinese.push_back(0);            ringMinus_ProductChinese.push_back(0);
        ringPlus_ProductChineseSpread.push_back(0);      ringMinus_ProductChineseSpread.push_back(0);
        ringPlus_countChinese.push_back(0);              ringMinus_countChinese.push_back(0);

        ringPlus_max_ProductRussian.push_back(0);        ringMinus_max_ProductRussian.push_back(0);
        ringPlus_min_ProductRussian.push_back(99);       ringMinus_min_ProductRussian.push_back(99);
        ringPlus_ProductRussian.push_back(0);            ringMinus_ProductRussian.push_back(0);
        ringPlus_ProductRussianSpread.push_back(0);      ringMinus_ProductRussianSpread.push_back(0);
        ringPlus_countRussian.push_back(0);              ringMinus_countRussian.push_back(0);


    }



    for(int ix=1; ix<=productEE->GetNbinsX(); ix++){
        for(int iy=1; iy<=productEE->GetNbinsY(); iy++){

            ix_ieta.push_back(ix%100);
            iy_iphi.push_back(iy%100);
            iz.push_back     (1*(ix>=100)-1*(ix<100));

            if (productEE->GetBinContent(ix,iy) > 0 ) {
                v_Product.push_back( productEE->GetBinContent(ix,iy) );
                if(productEEChinese->GetBinContent(ix,iy) > 0)
                {
                   v_ProductChinese.push_back( productEEChinese->GetBinContent(ix,iy) );
                }
                if(productEEChinese->GetBinContent(ix,iy) <= 0)
                {
                    v_ProductChinese.push_back( 0);
                }
                if(productEERussian->GetBinContent(ix,iy) > 0)
                {
                   v_ProductRussian.push_back( productEERussian->GetBinContent(ix,iy) );
                }
                if(productEERussian->GetBinContent(ix,iy) <= 0)
                {
                    v_ProductRussian.push_back( 0);
                }
            }
            else {
                v_Product.push_back( 0 );
                v_ProductChinese.push_back( 0);
                v_ProductRussian.push_back( 0);

            }
        }
    }





    for (int iter = 0; iter < (int) ix_ieta.size(); iter++) {
        if (iz.at(iter) != 0) {

            float dx = ix_ieta.at(iter) - 50;
            float dy = iy_iphi.at(iter) - 50;

            float ring = sqrt( dx*dx + dy*dy );

            int iring = round(ring) - 11;  //---- 12 [ = (62 - 50 - 1) from the 2D plot] is the first ring

            iring = -99;
            std::pair<int, int> ixiy (ix_ieta.at(iter), iy_iphi.at(iter));

            //       std::map<std::pair<int, int>,int>::const_iterator it = iring_map_plus.find(ixiy);
            //       if (it==iring_map_plus.end() ) {
            //         continue;
            //       }
            //

            if (iz.at(iter) > 0) {
                iring = iring_map_plus [ixiy];
            }
            else {
                iring = iring_map_minus [ixiy];
            }

            //       std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;

            //       if (iring < 0 ) continue;

            //       if (iring > 38 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << std::endl;
            //       if (iring > 37 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << std::endl;

            if (v_Product.at(iter) > 0 ) if (iring == 0 )  std::cout << " ix, iy, iring = " << ix_ieta.at(iter) << "  " <<  iy_iphi.at(iter)  << "  " << iring << " --> " << v_Product.at(iter) << " tot = " << ringPlus_count.at(iring) << std::endl;


            if (v_Product.at(iter) > 0 ) {
                if (iring > (50-11+2) || iring < 0) std::cout << " what ?!?   iring = " << iring << " dx = " << dx << " dy = " << dy << " :::: ix = " << ix_ieta.at(iter) << "  iy = " << iy_iphi.at(iter) << " prod = " << v_Product.at(iter) << std::endl;

                if (iz.at(iter) > 0) {

                    if (ringPlus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringPlus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringPlus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringPlus_min_Product.at(iring)   =  v_Product.at(iter);
                    }
                    ringPlus_Product.at(iring)       = ringPlus_Product.at(iring) + v_Product.at(iter);
                    ringPlus_ProductSpread.at(iring) = ringPlus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringPlus_count.at(iring)         = ringPlus_count.at(iring) + 1 ;

                    if (ringPlus_max_ProductChinese.at(iring) < v_ProductChinese.at(iter)) {
                        ringPlus_max_ProductChinese.at(iring)   =  v_ProductChinese.at(iter);
                    }
                    if (ringPlus_min_ProductChinese.at(iring) > v_ProductChinese.at(iter)) {
                        ringPlus_min_ProductChinese.at(iring)   =  v_ProductChinese.at(iter);
                    }

                    ringPlus_ProductChinese.at(iring)       = ringPlus_ProductChinese.at(iring) + v_ProductChinese.at(iter);
                    ringPlus_ProductChineseSpread.at(iring) = ringPlus_ProductChineseSpread.at(iring) + (v_ProductChinese.at(iter)) * (v_ProductChinese.at(iter));
                    ringPlus_countChinese.at(iring)         = ringPlus_countChinese.at(iring) + 1 ;

                    if (ringPlus_max_ProductRussian.at(iring) < v_ProductRussian.at(iter)) {
                        ringPlus_max_ProductRussian.at(iring)   =  v_ProductRussian.at(iter);
                    }
                    if (ringPlus_min_ProductRussian.at(iring) > v_ProductRussian.at(iter)) {
                        ringPlus_min_ProductRussian.at(iring)   =  v_ProductRussian.at(iter);
                    }

                    ringPlus_ProductRussian.at(iring)       = ringPlus_ProductRussian.at(iring) + v_ProductRussian.at(iter);
                    ringPlus_ProductRussianSpread.at(iring) = ringPlus_ProductRussianSpread.at(iring) + (v_ProductRussian.at(iter)) * (v_ProductRussian.at(iter));
                    ringPlus_countRussian.at(iring)         = ringPlus_countRussian.at(iring) + 1 ;

                }
                else {

                    if (ringMinus_max_Product.at(iring) < v_Product.at(iter)) {
                        ringMinus_max_Product.at(iring)   =  v_Product.at(iter);
                    }
                    if (ringMinus_min_Product.at(iring) > v_Product.at(iter)) {
                        ringMinus_min_Product.at(iring)   =  v_Product.at(iter);
                    }

                    ringMinus_Product.at(iring)       = ringMinus_Product.at(iring) + v_Product.at(iter);
                    ringMinus_ProductSpread.at(iring) = ringMinus_ProductSpread.at(iring) + (v_Product.at(iter)) * (v_Product.at(iter));
                    ringMinus_count.at(iring)         = ringMinus_count.at(iring) + 1 ;

                    if (ringMinus_max_ProductChinese.at(iring) < v_ProductChinese.at(iter)) {
                        ringMinus_max_ProductChinese.at(iring)   =  v_ProductChinese.at(iter);
                    }
                    if (ringMinus_min_ProductChinese.at(iring) > v_ProductChinese.at(iter)) {
                        ringMinus_min_ProductChinese.at(iring)   =  v_ProductChinese.at(iter);
                    }

                    ringMinus_ProductChinese.at(iring)       = ringMinus_ProductChinese.at(iring) + v_ProductChinese.at(iter);
                    ringMinus_ProductChineseSpread.at(iring) = ringMinus_ProductChineseSpread.at(iring) + (v_ProductChinese.at(iter)) * (v_ProductChinese.at(iter));
                    ringMinus_countChinese.at(iring)         = ringMinus_countChinese.at(iring) + 1 ;

                    if (ringMinus_max_ProductRussian.at(iring) < v_ProductRussian.at(iter)) {
                        ringMinus_max_ProductRussian.at(iring)   =  v_ProductRussian.at(iter);
                    }
                    if (ringMinus_min_ProductRussian.at(iring) > v_ProductRussian.at(iter)) {
                        ringMinus_min_ProductRussian.at(iring)   =  v_ProductRussian.at(iter);
                    }

                    ringMinus_ProductRussian.at(iring)       = ringMinus_ProductRussian.at(iring) + v_ProductRussian.at(iter);
                    ringMinus_ProductRussianSpread.at(iring) = ringMinus_ProductRussianSpread.at(iring) + (v_ProductRussian.at(iter)) * (v_ProductRussian.at(iter));
                    ringMinus_countRussian.at(iring)         = ringMinus_countRussian.at(iring) + 1 ;
                }
            }

        }
    }

    for (int iring = 0; iring < (50-11+3); iring++) {
        //     std::cout << " ringPlus_ProductSpread.at(" << iring << ") = " << ringPlus_ProductSpread.at(iring) << " ---> " << ringPlus_ProductSpread.at(iring)  << " - " <<  ringPlus_Product.at(iring)  * ringPlus_Product.at(iring) << std::endl;
        ringMinus_ProductSpread.at(iring) = sqrt(ringMinus_ProductSpread.at(iring) -  ringMinus_Product.at(iring) / ringMinus_count.at(iring) * ringMinus_Product.at(iring)) / ringMinus_count.at(iring);
        ringPlus_ProductSpread.at(iring)  = sqrt(ringPlus_ProductSpread.at(iring)  -  ringPlus_Product.at(iring)  / ringPlus_count.at(iring)  * ringPlus_Product.at(iring))  / ringPlus_count.at(iring);

        ringMinus_ProductChineseSpread.at(iring) = sqrt(ringMinus_ProductChineseSpread.at(iring) -  ringMinus_ProductChinese.at(iring) / ringMinus_countChinese.at(iring) * ringMinus_ProductChinese.at(iring)) / ringMinus_countChinese.at(iring);
        ringPlus_ProductChineseSpread.at(iring)  = sqrt(ringPlus_ProductChineseSpread.at(iring)  -  ringPlus_ProductChinese.at(iring)  / ringPlus_countChinese.at(iring)  * ringPlus_ProductChinese.at(iring))  / ringPlus_countChinese.at(iring);

        ringMinus_ProductRussianSpread.at(iring) = sqrt(ringMinus_ProductRussianSpread.at(iring) -  ringMinus_ProductRussian.at(iring) / ringMinus_countRussian.at(iring) * ringMinus_ProductRussian.at(iring)) / ringMinus_countRussian.at(iring);
        ringPlus_ProductRussianSpread.at(iring)  = sqrt(ringPlus_ProductRussianSpread.at(iring)  -  ringPlus_ProductRussian.at(iring)  / ringPlus_countRussian.at(iring)  * ringPlus_ProductRussian.at(iring))  / ringPlus_countRussian.at(iring);


    }


    TGraphErrors *gr_EEPlus_ET_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_ET_eta_Product = new TGraphErrors();

    TGraphErrors *gr_EEPlus_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_eta_Product = new TGraphErrors();
    TGraphErrors *gr_EEPlus_max_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_eta_Product = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_eta_Product = new TGraphErrors();

    TGraphErrors *gr_ring_EEPlus_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_eta_Product = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_max_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_max_eta_Product = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_min_eta_Product = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_min_eta_Product = new TGraphErrors();


    TGraphErrors *gr_EEPlus_max_ET_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_ET_eta_Product = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_ET_eta_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_ET_eta_Product = new TGraphErrors();



    TGraphErrors *gr_EEPlus_Product = new TGraphErrors();     TGraphErrors *gr_EEMinus_Product = new TGraphErrors();
    TGraphErrors *gr_EEPlus_ProductSpread = new TGraphErrors();     TGraphErrors *gr_EEMinus_ProductSpread = new TGraphErrors();


    TGraphErrors *gr_EEPlus_ET_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_ET_eta_ProductChinese = new TGraphErrors();

    TGraphErrors *gr_EEPlus_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_eta_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_EEPlus_max_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_eta_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_eta_ProductChinese = new TGraphErrors();

    TGraphErrors *gr_ring_EEPlus_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_eta_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_max_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_max_eta_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_min_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_min_eta_ProductChinese = new TGraphErrors();


    TGraphErrors *gr_EEPlus_max_ET_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_ET_eta_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_ET_eta_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_ET_eta_ProductChinese = new TGraphErrors();



    TGraphErrors *gr_EEPlus_ProductChinese = new TGraphErrors();     TGraphErrors *gr_EEMinus_ProductChinese = new TGraphErrors();
    TGraphErrors *gr_EEPlus_ProductChineseSpread = new TGraphErrors();     TGraphErrors *gr_EEMinus_ProductChineseSpread = new TGraphErrors();


    TGraphErrors *gr_EEPlus_ET_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_ET_eta_ProductRussian = new TGraphErrors();

    TGraphErrors *gr_EEPlus_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_eta_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_EEPlus_max_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_eta_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_eta_ProductRussian = new TGraphErrors();

    TGraphErrors *gr_ring_EEPlus_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_eta_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_max_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_max_eta_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_ring_EEPlus_min_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_ring_EEMinus_min_eta_ProductRussian = new TGraphErrors();


    TGraphErrors *gr_EEPlus_max_ET_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_max_ET_eta_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_EEPlus_min_ET_eta_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_min_ET_eta_ProductRussian = new TGraphErrors();



    TGraphErrors *gr_EEPlus_ProductRussian = new TGraphErrors();     TGraphErrors *gr_EEMinus_ProductRussian = new TGraphErrors();
    TGraphErrors *gr_EEPlus_ProductRussianSpread = new TGraphErrors();     TGraphErrors *gr_EEMinus_ProductRussianSpread = new TGraphErrors();


    for (int iter = 0; iter < (39); iter++) {
        //     for (int iter = 0; iter < (50-11+3); iter++) {

        //std::cout << " ringPlus_Product.at(" << iter << ") = " << ringPlus_Product.at(iter) << " N = " << ringPlus_count.at(iter) <<  " --> " << ringPlus_Product.at(iter) / ringPlus_count.at(iter) << std::endl;

        gr_EEPlus_Product-> SetPoint (iter,  iter,   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_EEMinus_Product-> SetPoint (iter,  iter,   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;

        gr_EEPlus_ProductSpread-> SetPoint (iter,   iter,   ringPlus_count.at(iter) ?  ringPlus_ProductSpread.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_EEMinus_ProductSpread-> SetPoint (iter,   iter,   ringMinus_count.at(iter) ?  ringMinus_ProductSpread.at(iter) / ringMinus_count.at(iter) : 0 ) ;

        gr_EEPlus_Product-> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter) ) ;
        gr_EEMinus_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) ) ;


        gr_EEPlus_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_EEMinus_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;

        gr_EEPlus_eta_Product -> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter) ) ;
        gr_EEMinus_eta_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) ) ;


        gr_EEPlus_max_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_Product.at(iter)) ;
        gr_EEMinus_max_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_Product.at(iter) ) ;

        gr_EEPlus_min_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_Product.at(iter)) ;
        gr_EEMinus_min_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_Product.at(iter) ) ;




        //     gr_ring_EEPlus_eta_Product -> SetPoint (iter,  42-3 -  iter,   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        //     gr_ring_EEMinus_eta_Product-> SetPoint (iter,  42-3 -  iter,   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;
        //
        //     gr_ring_EEPlus_eta_Product -> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter) ) ;
        //     gr_ring_EEMinus_eta_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) ) ;
        //
        //
        //     gr_ring_EEPlus_max_eta_Product -> SetPoint (iter,  42-3 -  iter,   ringPlus_max_Product.at(iter)) ;
        //     gr_ring_EEMinus_max_eta_Product-> SetPoint (iter,  42-3 -  iter,   ringMinus_max_Product.at(iter) ) ;
        //
        //     gr_ring_EEPlus_min_eta_Product -> SetPoint (iter,  42-3 -  iter,   ringPlus_min_Product.at(iter)) ;
        //     gr_ring_EEMinus_min_eta_Product-> SetPoint (iter,  42-3 -  iter,   ringMinus_min_Product.at(iter) ) ;

        gr_ring_EEPlus_eta_Product -> SetPoint (iter,  iter,   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 ) ;
        gr_ring_EEMinus_eta_Product-> SetPoint (iter,  iter,   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 ) ;

        gr_ring_EEPlus_eta_Product -> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter) ) ;
        gr_ring_EEMinus_eta_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) ) ;


        gr_ring_EEPlus_max_eta_Product -> SetPoint (iter,  iter,   ringPlus_max_Product.at(iter)) ;
        gr_ring_EEMinus_max_eta_Product-> SetPoint (iter,  iter,   ringMinus_max_Product.at(iter) ) ;

        gr_ring_EEPlus_min_eta_Product -> SetPoint (iter,  iter,   ringPlus_min_Product.at(iter)) ;
        gr_ring_EEMinus_min_eta_Product-> SetPoint (iter,  iter,   ringMinus_min_Product.at(iter) ) ;




        //---- p = pt * cosh(eta)
        gr_EEPlus_ET_eta_Product-> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter)    / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;
        gr_EEMinus_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;

        gr_EEPlus_ET_eta_Product-> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter)  / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;
        gr_EEMinus_ET_eta_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_ET_eta_Product-> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_count.at(iter) ? ringPlus_Product.at(iter) / ringPlus_count.at(iter) : 0 );
        //gr_EEMinus_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_count.at(iter) ?  ringMinus_Product.at(iter) / ringMinus_count.at(iter) : 0 );
        //gr_EEPlus_ET_eta_Product-> SetPointError (iter,    0,  ringPlus_ProductSpread.at(iter) );
        //gr_EEMinus_ET_eta_Product-> SetPointError (iter,   0,  ringMinus_ProductSpread.at(iter) );


        //---- p = pt * cosh(eta)
        gr_EEPlus_max_ET_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_Product.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_max_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_Product.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        gr_EEPlus_min_ET_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_Product.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_min_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_Product.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_max_ET_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_max_Product.at(iter));
        //gr_EEMinus_max_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_max_Product.at(iter));

        //gr_EEPlus_min_ET_eta_Product -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_min_Product.at(iter));
        //gr_EEMinus_min_ET_eta_Product-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_min_Product.at(iter));











        gr_EEPlus_ProductChinese-> SetPoint (iter,  iter,   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter) : 0 ) ;
        gr_EEMinus_ProductChinese-> SetPoint (iter,  iter,   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) : 0 ) ;

        gr_EEPlus_ProductChineseSpread-> SetPoint (iter,   iter,   ringPlus_countChinese.at(iter) ?  ringPlus_ProductChineseSpread.at(iter) / ringPlus_countChinese.at(iter) : 0 ) ;
        gr_EEMinus_ProductChineseSpread-> SetPoint (iter,   iter,   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChineseSpread.at(iter) / ringMinus_countChinese.at(iter) : 0 ) ;

        gr_EEPlus_ProductChinese-> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter) ) ;
        gr_EEMinus_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) ) ;


        gr_EEPlus_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter) : 0 ) ;
        gr_EEMinus_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) : 0 ) ;

        gr_EEPlus_eta_ProductChinese -> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter) ) ;
        gr_EEMinus_eta_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) ) ;


        gr_EEPlus_max_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_ProductChinese.at(iter)) ;
        gr_EEMinus_max_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_ProductChinese.at(iter) ) ;

        gr_EEPlus_min_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_ProductChinese.at(iter)) ;
        gr_EEMinus_min_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_ProductChinese.at(iter) ) ;




        //     gr_ring_EEPlus_eta_ProductChinese -> SetPoint (iter,  42-3 -  iter,   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter) : 0 ) ;
        //     gr_ring_EEMinus_eta_ProductChinese-> SetPoint (iter,  42-3 -  iter,   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) : 0 ) ;
        //
        //     gr_ring_EEPlus_eta_ProductChinese -> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter) ) ;
        //     gr_ring_EEMinus_eta_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) ) ;
        //
        //
        //     gr_ring_EEPlus_max_eta_ProductChinese -> SetPoint (iter,  42-3 -  iter,   ringPlus_max_ProductChinese.at(iter)) ;
        //     gr_ring_EEMinus_max_eta_ProductChinese-> SetPoint (iter,  42-3 -  iter,   ringMinus_max_ProductChinese.at(iter) ) ;
        //
        //     gr_ring_EEPlus_min_eta_ProductChinese -> SetPoint (iter,  42-3 -  iter,   ringPlus_min_ProductChinese.at(iter)) ;
        //     gr_ring_EEMinus_min_eta_ProductChinese-> SetPoint (iter,  42-3 -  iter,   ringMinus_min_ProductChinese.at(iter) ) ;

        gr_ring_EEPlus_eta_ProductChinese -> SetPoint (iter,  iter,   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter) : 0 ) ;
        gr_ring_EEMinus_eta_ProductChinese-> SetPoint (iter,  iter,   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) : 0 ) ;

        gr_ring_EEPlus_eta_ProductChinese -> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter) ) ;
        gr_ring_EEMinus_eta_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) ) ;


        gr_ring_EEPlus_max_eta_ProductChinese -> SetPoint (iter,  iter,   ringPlus_max_ProductChinese.at(iter)) ;
        gr_ring_EEMinus_max_eta_ProductChinese-> SetPoint (iter,  iter,   ringMinus_max_ProductChinese.at(iter) ) ;

        gr_ring_EEPlus_min_eta_ProductChinese -> SetPoint (iter,  iter,   ringPlus_min_ProductChinese.at(iter)) ;
        gr_ring_EEMinus_min_eta_ProductChinese-> SetPoint (iter,  iter,   ringMinus_min_ProductChinese.at(iter) ) ;




        //---- p = pt * cosh(eta)
        gr_EEPlus_ET_eta_ProductChinese-> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter)    / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;
        gr_EEMinus_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;

        gr_EEPlus_ET_eta_ProductChinese-> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter)  / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;
        gr_EEMinus_ET_eta_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_ET_eta_ProductChinese-> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_countChinese.at(iter) ? ringPlus_ProductChinese.at(iter) / ringPlus_countChinese.at(iter) : 0 );
        //gr_EEMinus_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_countChinese.at(iter) ?  ringMinus_ProductChinese.at(iter) / ringMinus_countChinese.at(iter) : 0 );
        //gr_EEPlus_ET_eta_ProductChinese-> SetPointError (iter,    0,  ringPlus_ProductChineseSpread.at(iter) );
        //gr_EEMinus_ET_eta_ProductChinese-> SetPointError (iter,   0,  ringMinus_ProductChineseSpread.at(iter) );


        //---- p = pt * cosh(eta)
        gr_EEPlus_max_ET_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_ProductChinese.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_max_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_ProductChinese.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        gr_EEPlus_min_ET_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_ProductChinese.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_min_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_ProductChinese.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_max_ET_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_max_ProductChinese.at(iter));
        //gr_EEMinus_max_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_max_ProductChinese.at(iter));

        //gr_EEPlus_min_ET_eta_ProductChinese -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_min_ProductChinese.at(iter));
        //gr_EEMinus_min_ET_eta_ProductChinese-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_min_ProductChinese.at(iter));













        gr_EEPlus_ProductRussian-> SetPoint (iter,  iter,   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter) : 0 ) ;
        gr_EEMinus_ProductRussian-> SetPoint (iter,  iter,   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) : 0 ) ;

        gr_EEPlus_ProductRussianSpread-> SetPoint (iter,   iter,   ringPlus_countRussian.at(iter) ?  ringPlus_ProductRussianSpread.at(iter) / ringPlus_countRussian.at(iter) : 0 ) ;
        gr_EEMinus_ProductRussianSpread-> SetPoint (iter,   iter,   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussianSpread.at(iter) / ringMinus_countRussian.at(iter) : 0 ) ;

        gr_EEPlus_ProductRussian-> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter) ) ;
        gr_EEMinus_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) ) ;


        gr_EEPlus_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter) : 0 ) ;
        gr_EEMinus_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) : 0 ) ;

        gr_EEPlus_eta_ProductRussian -> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter) ) ;
        gr_EEMinus_eta_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) ) ;


        gr_EEPlus_max_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_ProductRussian.at(iter)) ;
        gr_EEMinus_max_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_ProductRussian.at(iter) ) ;

        gr_EEPlus_min_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_ProductRussian.at(iter)) ;
        gr_EEMinus_min_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_ProductRussian.at(iter) ) ;




        //     gr_ring_EEPlus_eta_ProductRussian -> SetPoint (iter,  42-3 -  iter,   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter) : 0 ) ;
        //     gr_ring_EEMinus_eta_ProductRussian-> SetPoint (iter,  42-3 -  iter,   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) : 0 ) ;
        //
        //     gr_ring_EEPlus_eta_ProductRussian -> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter) ) ;
        //     gr_ring_EEMinus_eta_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) ) ;
        //
        //
        //     gr_ring_EEPlus_max_eta_ProductRussian -> SetPoint (iter,  42-3 -  iter,   ringPlus_max_ProductRussian.at(iter)) ;
        //     gr_ring_EEMinus_max_eta_ProductRussian-> SetPoint (iter,  42-3 -  iter,   ringMinus_max_ProductRussian.at(iter) ) ;
        //
        //     gr_ring_EEPlus_min_eta_ProductRussian -> SetPoint (iter,  42-3 -  iter,   ringPlus_min_ProductRussian.at(iter)) ;
        //     gr_ring_EEMinus_min_eta_ProductRussian-> SetPoint (iter,  42-3 -  iter,   ringMinus_min_ProductRussian.at(iter) ) ;

        gr_ring_EEPlus_eta_ProductRussian -> SetPoint (iter,  iter,   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter) : 0 ) ;
        gr_ring_EEMinus_eta_ProductRussian-> SetPoint (iter,  iter,   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) : 0 ) ;

        gr_ring_EEPlus_eta_ProductRussian -> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter) ) ;
        gr_ring_EEMinus_eta_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) ) ;


        gr_ring_EEPlus_max_eta_ProductRussian -> SetPoint (iter,  iter,   ringPlus_max_ProductRussian.at(iter)) ;
        gr_ring_EEMinus_max_eta_ProductRussian-> SetPoint (iter,  iter,   ringMinus_max_ProductRussian.at(iter) ) ;

        gr_ring_EEPlus_min_eta_ProductRussian -> SetPoint (iter,  iter,   ringPlus_min_ProductRussian.at(iter)) ;
        gr_ring_EEMinus_min_eta_ProductRussian-> SetPoint (iter,  iter,   ringMinus_min_ProductRussian.at(iter) ) ;




        //---- p = pt * cosh(eta)
        gr_EEPlus_ET_eta_ProductRussian-> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter)    / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;
        gr_EEMinus_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) / cosh(-log(tan(atan(1711./3630.*(iter+10)/50.)/2))) : 0 ) ;

        gr_EEPlus_ET_eta_ProductRussian-> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter)  / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;
        gr_EEMinus_ET_eta_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) / cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_ET_eta_ProductRussian-> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_countRussian.at(iter) ? ringPlus_ProductRussian.at(iter) / ringPlus_countRussian.at(iter) : 0 );
        //gr_EEMinus_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_countRussian.at(iter) ?  ringMinus_ProductRussian.at(iter) / ringMinus_countRussian.at(iter) : 0 );
        //gr_EEPlus_ET_eta_ProductRussian-> SetPointError (iter,    0,  ringPlus_ProductRussianSpread.at(iter) );
        //gr_EEMinus_ET_eta_ProductRussian-> SetPointError (iter,   0,  ringMinus_ProductRussianSpread.at(iter) );


        //---- p = pt * cosh(eta)
        gr_EEPlus_max_ET_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_max_ProductRussian.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_max_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_max_ProductRussian.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        gr_EEPlus_min_ET_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringPlus_min_ProductRussian.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)))) ;
        gr_EEMinus_min_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2)),   ringMinus_min_ProductRussian.at(iter)/ cosh(-log(tan(atan(1711./3630.*(39-iter+10)/50.)/2))) ) ;

        //gr_EEPlus_max_ET_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_max_ProductRussian.at(iter));
        //gr_EEMinus_max_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_max_ProductRussian.at(iter));

        //gr_EEPlus_min_ET_eta_ProductRussian -> SetPoint (iter,  -log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringPlus_min_ProductRussian.at(iter));
        //gr_EEMinus_min_ET_eta_ProductRussian-> SetPoint (iter,-log(tan(atan(1711./3630.*(iter+10)/50.)/2)),   ringMinus_min_ProductRussian.at(iter));


    }


    //   ---- style ----

    gr_EEPlus_Product->SetMarkerSize  (1);                           gr_EEMinus_Product->SetMarkerSize  (1);
    gr_EEPlus_Product->SetMarkerStyle (24);                          gr_EEMinus_Product->SetMarkerStyle (22);
    gr_EEPlus_Product->SetMarkerColor (kRed);                        gr_EEMinus_Product->SetMarkerColor (kRed);
    gr_EEPlus_Product->SetLineWidth (1);                             gr_EEMinus_Product->SetLineWidth (1);
    gr_EEPlus_Product->SetLineColor (kRed);                          gr_EEMinus_Product->SetLineColor (kRed);

    gr_EEPlus_ProductSpread->SetMarkerSize  (1);                           gr_EEMinus_ProductSpread->SetMarkerSize  (1);
    gr_EEPlus_ProductSpread->SetMarkerStyle (24);                          gr_EEMinus_ProductSpread->SetMarkerStyle (22);
    gr_EEPlus_ProductSpread->SetMarkerColor (kRed);                        gr_EEMinus_ProductSpread->SetMarkerColor (kRed);
    gr_EEPlus_ProductSpread->SetLineWidth (1);                             gr_EEMinus_ProductSpread->SetLineWidth (1);
    gr_EEPlus_ProductSpread->SetLineColor (kRed);                          gr_EEMinus_ProductSpread->SetLineColor (kRed);

    gr_EEPlus_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_eta_Product->SetMarkerColor (kRed);                        gr_EEMinus_eta_Product->SetMarkerColor (kRed);
    gr_EEPlus_eta_Product->SetLineWidth (1);                             gr_EEMinus_eta_Product->SetLineWidth (1);
    gr_EEPlus_eta_Product->SetLineColor (kRed);                          gr_EEMinus_eta_Product->SetLineColor (kRed);

    gr_ring_EEPlus_eta_Product->SetMarkerSize  (1);                           gr_ring_EEMinus_eta_Product->SetMarkerSize  (1);
    gr_ring_EEPlus_eta_Product->SetMarkerStyle (24);                          gr_ring_EEMinus_eta_Product->SetMarkerStyle (22);
    gr_ring_EEPlus_eta_Product->SetMarkerColor (kRed);                        gr_ring_EEMinus_eta_Product->SetMarkerColor (kRed);
    gr_ring_EEPlus_eta_Product->SetLineWidth (1);                             gr_ring_EEMinus_eta_Product->SetLineWidth (1);
    gr_ring_EEPlus_eta_Product->SetLineColor (kRed);                          gr_ring_EEMinus_eta_Product->SetLineColor (kRed);

    gr_EEPlus_ET_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_ET_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_ET_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_ET_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_ET_eta_Product->SetMarkerColor (kRed);                        gr_EEMinus_ET_eta_Product->SetMarkerColor (kRed);
    gr_EEPlus_ET_eta_Product->SetLineWidth (1);                             gr_EEMinus_ET_eta_Product->SetLineWidth (1);
    gr_EEPlus_ET_eta_Product->SetLineColor (kRed);                          gr_EEMinus_ET_eta_Product->SetLineColor (kRed);


    gr_EEPlus_max_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_max_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_max_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_max_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_max_eta_Product->SetMarkerColor (kBlue);                        gr_EEMinus_max_eta_Product->SetMarkerColor (kBlue);
    gr_EEPlus_max_eta_Product->SetLineWidth (1);                             gr_EEMinus_max_eta_Product->SetLineWidth (1);
    gr_EEPlus_max_eta_Product->SetLineColor (kBlue);                          gr_EEMinus_max_eta_Product->SetLineColor (kBlue);

    gr_EEPlus_min_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_min_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_min_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_min_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_min_eta_Product->SetMarkerColor (kBlue);                        gr_EEMinus_min_eta_Product->SetMarkerColor (kBlue);
    gr_EEPlus_min_eta_Product->SetLineWidth (1);                             gr_EEMinus_min_eta_Product->SetLineWidth (1);
    gr_EEPlus_min_eta_Product->SetLineColor (kBlue);                          gr_EEMinus_min_eta_Product->SetLineColor (kBlue);


    gr_ring_EEPlus_max_eta_Product->SetMarkerSize  (1);                           gr_ring_EEMinus_max_eta_Product->SetMarkerSize  (1);
    gr_ring_EEPlus_max_eta_Product->SetMarkerStyle (24);                          gr_ring_EEMinus_max_eta_Product->SetMarkerStyle (22);
    gr_ring_EEPlus_max_eta_Product->SetMarkerColor (kBlue);                        gr_ring_EEMinus_max_eta_Product->SetMarkerColor (kBlue);
    gr_ring_EEPlus_max_eta_Product->SetLineWidth (1);                             gr_ring_EEMinus_max_eta_Product->SetLineWidth (1);
    gr_ring_EEPlus_max_eta_Product->SetLineColor (kBlue);                          gr_ring_EEMinus_max_eta_Product->SetLineColor (kBlue);

    gr_ring_EEPlus_min_eta_Product->SetMarkerSize  (1);                           gr_ring_EEMinus_min_eta_Product->SetMarkerSize  (1);
    gr_ring_EEPlus_min_eta_Product->SetMarkerStyle (24);                          gr_ring_EEMinus_min_eta_Product->SetMarkerStyle (22);
    gr_ring_EEPlus_min_eta_Product->SetMarkerColor (kBlue);                        gr_ring_EEMinus_min_eta_Product->SetMarkerColor (kBlue);
    gr_ring_EEPlus_min_eta_Product->SetLineWidth (1);                             gr_ring_EEMinus_min_eta_Product->SetLineWidth (1);
    gr_ring_EEPlus_min_eta_Product->SetLineColor (kBlue);                          gr_ring_EEMinus_min_eta_Product->SetLineColor (kBlue);


    gr_EEPlus_max_ET_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_max_ET_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_max_ET_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_max_ET_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_max_ET_eta_Product->SetMarkerColor (kBlue);                        gr_EEMinus_max_ET_eta_Product->SetMarkerColor (kBlue);
    gr_EEPlus_max_ET_eta_Product->SetLineWidth (1);                             gr_EEMinus_max_ET_eta_Product->SetLineWidth (1);
    gr_EEPlus_max_ET_eta_Product->SetLineColor (kBlue);                          gr_EEMinus_max_ET_eta_Product->SetLineColor (kBlue);

    gr_EEPlus_min_ET_eta_Product->SetMarkerSize  (1);                           gr_EEMinus_min_ET_eta_Product->SetMarkerSize  (1);
    gr_EEPlus_min_ET_eta_Product->SetMarkerStyle (24);                          gr_EEMinus_min_ET_eta_Product->SetMarkerStyle (22);
    gr_EEPlus_min_ET_eta_Product->SetMarkerColor (kBlue);                        gr_EEMinus_min_ET_eta_Product->SetMarkerColor (kBlue);
    gr_EEPlus_min_ET_eta_Product->SetLineWidth (1);                             gr_EEMinus_min_ET_eta_Product->SetLineWidth (1);
    gr_EEPlus_min_ET_eta_Product->SetLineColor (kBlue);                          gr_EEMinus_min_ET_eta_Product->SetLineColor (kBlue);


    //   ---- style (end) ----

    TLegend* legend = new TLegend(0.1,0.7,0.48,0.9);
      legend->AddEntry(gr_EEMinus_Product,"EE-","lep");
      legend->AddEntry(gr_EEPlus_Product, "EE+","lep");


      TCanvas* ccRing = new TCanvas ("ccRing","",800,600);
      ccRing->Divide(2,1);

      ccRing->cd(1);

      gr_EEMinus_Product->Draw("APL");
      gr_EEPlus_Product->Draw("PL");
      legend->Draw();


      ccRing->cd(2);

      gr_EEPlus_ProductSpread->Draw("APL");

      gr_EEMinus_ProductSpread->Draw("PL");


      gr_EEMinus_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_EEMinus_Product->GetXaxis()->SetTitle("iRing");

      gr_EEPlus_ProductSpread->GetYaxis()->SetTitle("Noise RMS [GeV]");
      gr_EEPlus_ProductSpread->GetXaxis()->SetTitle("iRing");

      ccRing->SaveAs("NoiseRing_EE.png");
      ccRing->SaveAs("NoiseRing_EE.root");


      TCanvas* ccRing2 = new TCanvas ("ccRing2","",800,600);
      gr_EEMinus_Product->Draw("APL");
      gr_EEPlus_Product->Draw("PL");
      legend->Draw();
      ccRing2->SaveAs("NoiseRing2_EE.png");
      ccRing2->SaveAs("NoiseRing2_EE.root");



      TCanvas* ccRing_eta = new TCanvas ("ccRing_eta","",800,600);
      gr_EEMinus_eta_Product->Draw("APL");
      gr_EEPlus_eta_Product->Draw("PL");
      legend->Draw();

      gr_EEMinus_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_EEMinus_eta_Product->GetXaxis()->SetTitle("#eta");

      ccRing_eta->SaveAs("NoiseRing_eta_EE.png");
      ccRing_eta->SaveAs("NoiseRing_eta_EE.root");

      TF1 *f1 = new TF1("f1","[0]+[1]*exp([2]*x)",1.5,3.1);
      f1->SetParameter(0,0.1);
      f1->SetParameter(1,0.);
      f1->SetParameter(2,0.1);
      //cout << "Fit to minus, average" << endl;
      //gr_EEMinus_eta_Product->Fit("f1","","",2.4,3.0);
      //cout << endl;
      //cout << "Fit to plus, average" << endl;
      //gr_EEPlus_eta_Product->Fit("expo","","",2.3,3.0);

      TCanvas* ccRing_eta_minmax = new TCanvas ("ccRing_eta_minmax","",800,600);

      gr_EEPlus_max_eta_Product->Draw("AP");

      gr_EEMinus_eta_Product->Draw("PL");
      gr_EEPlus_eta_Product->Draw("PL");

      gr_EEPlus_min_eta_Product->Draw("P");
      gr_EEMinus_min_eta_Product->Draw("P");
      gr_EEMinus_max_eta_Product->Draw("P");

      legend->Draw();

      gr_EEPlus_max_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_EEPlus_max_eta_Product->GetXaxis()->SetTitle("#eta");

      ccRing_eta_minmax->SaveAs("NoiseRing_eta_minmax_EE.png");
      ccRing_eta_minmax->SaveAs("NoiseRing_eta_minmax_EE.root");




      TLegend* legend2 = new TLegend(0.1,0.7,0.48,0.9);
      legend2->AddEntry(gr_EEMinus_ET_eta_Product,"average EE-","lep");
      legend2->AddEntry(gr_EEPlus_ET_eta_Product, "average EE+","lep");
      legend2->AddEntry(gr_EEMinus_min_ET_eta_Product, "min/max EE-","lep");
      legend2->AddEntry(gr_EEPlus_max_ET_eta_Product, "min/max EE+","lep");



      TCanvas* ccRing_ring_minmax = new TCanvas ("ccRing_ring_minmax","",800,600);

      gr_ring_EEPlus_max_eta_Product->Draw("AP");

      gr_ring_EEMinus_eta_Product->Draw("PL");
      gr_ring_EEPlus_eta_Product->Draw("PL");

      gr_ring_EEPlus_min_eta_Product->Draw("P");
      gr_ring_EEMinus_min_eta_Product->Draw("P");
      gr_ring_EEMinus_max_eta_Product->Draw("P");

      legend2->Draw();

      gr_ring_EEPlus_max_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_ring_EEPlus_max_eta_Product->GetXaxis()->SetTitle("ring");

      ccRing_ring_minmax->SaveAs("NoiseRing_ring_minmax_ring_EE.png");
      ccRing_ring_minmax->SaveAs("NoiseRing_ring_minmax_ring_EE.root");





      TCanvas* ccRing_ET_eta = new TCanvas ("ccRing_ET_eta","",800,600);
      gr_EEMinus_ET_eta_Product->Draw("APL");
      gr_EEPlus_ET_eta_Product->Draw("PL");
      legend->Draw();

      gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
      //gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_EEMinus_ET_eta_Product->GetXaxis()->SetTitle("#eta");

      ccRing_ET_eta->SaveAs("NoiseRing_ET_eta_EE.png");
      ccRing_ET_eta->SaveAs("NoiseRing_ET_eta_EE.root");






      TCanvas* ccRing_ET_eta_minmax = new TCanvas ("ccRing_ET_eta_minmax","",800,600);

      gr_EEPlus_max_ET_eta_Product->Draw("AP");

      gr_EEMinus_ET_eta_Product->Draw("PL");
      gr_EEPlus_ET_eta_Product->Draw("PL");

      gr_EEPlus_min_ET_eta_Product->Draw("P");
      gr_EEMinus_min_ET_eta_Product->Draw("P");
      gr_EEMinus_max_ET_eta_Product->Draw("P");

      legend2->Draw();

      gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
      //gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
      gr_EEPlus_max_ET_eta_Product->GetXaxis()->SetTitle("#eta");

      ccRing_ET_eta_minmax->SaveAs("NoiseRing_ET_eta_minmax_EE.png");
      ccRing_ET_eta_minmax->SaveAs("NoiseRing_ET_eta_minmax_EE.root");


//Chinese


    //   ---- style ----

    gr_EEPlus_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_ProductChinese->SetMarkerColor (kRed);                        gr_EEMinus_ProductChinese->SetMarkerColor (kRed);
    gr_EEPlus_ProductChinese->SetLineWidth (1);                             gr_EEMinus_ProductChinese->SetLineWidth (1);
    gr_EEPlus_ProductChinese->SetLineColor (kRed);                          gr_EEMinus_ProductChinese->SetLineColor (kRed);

    gr_EEPlus_ProductChineseSpread->SetMarkerSize  (1);                           gr_EEMinus_ProductChineseSpread->SetMarkerSize  (1);
    gr_EEPlus_ProductChineseSpread->SetMarkerStyle (24);                          gr_EEMinus_ProductChineseSpread->SetMarkerStyle (22);
    gr_EEPlus_ProductChineseSpread->SetMarkerColor (kRed);                        gr_EEMinus_ProductChineseSpread->SetMarkerColor (kRed);
    gr_EEPlus_ProductChineseSpread->SetLineWidth (1);                             gr_EEMinus_ProductChineseSpread->SetLineWidth (1);
    gr_EEPlus_ProductChineseSpread->SetLineColor (kRed);                          gr_EEMinus_ProductChineseSpread->SetLineColor (kRed);

    gr_EEPlus_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_eta_ProductChinese->SetMarkerColor (kRed);                        gr_EEMinus_eta_ProductChinese->SetMarkerColor (kRed);
    gr_EEPlus_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_eta_ProductChinese->SetLineColor (kRed);                          gr_EEMinus_eta_ProductChinese->SetLineColor (kRed);

    gr_ring_EEPlus_eta_ProductChinese->SetMarkerSize  (1);                           gr_ring_EEMinus_eta_ProductChinese->SetMarkerSize  (1);
    gr_ring_EEPlus_eta_ProductChinese->SetMarkerStyle (24);                          gr_ring_EEMinus_eta_ProductChinese->SetMarkerStyle (22);
    gr_ring_EEPlus_eta_ProductChinese->SetMarkerColor (kRed);                        gr_ring_EEMinus_eta_ProductChinese->SetMarkerColor (kRed);
    gr_ring_EEPlus_eta_ProductChinese->SetLineWidth (1);                             gr_ring_EEMinus_eta_ProductChinese->SetLineWidth (1);
    gr_ring_EEPlus_eta_ProductChinese->SetLineColor (kRed);                          gr_ring_EEMinus_eta_ProductChinese->SetLineColor (kRed);

    gr_EEPlus_ET_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_ET_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_ET_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_ET_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_ET_eta_ProductChinese->SetMarkerColor (kRed);                        gr_EEMinus_ET_eta_ProductChinese->SetMarkerColor (kRed);
    gr_EEPlus_ET_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_ET_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_ET_eta_ProductChinese->SetLineColor (kRed);                          gr_EEMinus_ET_eta_ProductChinese->SetLineColor (kRed);


    gr_EEPlus_max_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_max_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_max_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_max_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_max_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_EEMinus_max_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_EEPlus_max_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_max_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_max_eta_ProductChinese->SetLineColor (kBlue);                          gr_EEMinus_max_eta_ProductChinese->SetLineColor (kBlue);

    gr_EEPlus_min_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_min_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_min_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_min_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_min_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_EEMinus_min_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_EEPlus_min_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_min_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_min_eta_ProductChinese->SetLineColor (kBlue);                          gr_EEMinus_min_eta_ProductChinese->SetLineColor (kBlue);


    gr_ring_EEPlus_max_eta_ProductChinese->SetMarkerSize  (1);                           gr_ring_EEMinus_max_eta_ProductChinese->SetMarkerSize  (1);
    gr_ring_EEPlus_max_eta_ProductChinese->SetMarkerStyle (24);                          gr_ring_EEMinus_max_eta_ProductChinese->SetMarkerStyle (22);
    gr_ring_EEPlus_max_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_ring_EEMinus_max_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_ring_EEPlus_max_eta_ProductChinese->SetLineWidth (1);                             gr_ring_EEMinus_max_eta_ProductChinese->SetLineWidth (1);
    gr_ring_EEPlus_max_eta_ProductChinese->SetLineColor (kBlue);                          gr_ring_EEMinus_max_eta_ProductChinese->SetLineColor (kBlue);

    gr_ring_EEPlus_min_eta_ProductChinese->SetMarkerSize  (1);                           gr_ring_EEMinus_min_eta_ProductChinese->SetMarkerSize  (1);
    gr_ring_EEPlus_min_eta_ProductChinese->SetMarkerStyle (24);                          gr_ring_EEMinus_min_eta_ProductChinese->SetMarkerStyle (22);
    gr_ring_EEPlus_min_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_ring_EEMinus_min_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_ring_EEPlus_min_eta_ProductChinese->SetLineWidth (1);                             gr_ring_EEMinus_min_eta_ProductChinese->SetLineWidth (1);
    gr_ring_EEPlus_min_eta_ProductChinese->SetLineColor (kBlue);                          gr_ring_EEMinus_min_eta_ProductChinese->SetLineColor (kBlue);


    gr_EEPlus_max_ET_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_max_ET_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_max_ET_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_max_ET_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_max_ET_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_EEMinus_max_ET_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_EEPlus_max_ET_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_max_ET_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_max_ET_eta_ProductChinese->SetLineColor (kBlue);                          gr_EEMinus_max_ET_eta_ProductChinese->SetLineColor (kBlue);

    gr_EEPlus_min_ET_eta_ProductChinese->SetMarkerSize  (1);                           gr_EEMinus_min_ET_eta_ProductChinese->SetMarkerSize  (1);
    gr_EEPlus_min_ET_eta_ProductChinese->SetMarkerStyle (24);                          gr_EEMinus_min_ET_eta_ProductChinese->SetMarkerStyle (22);
    gr_EEPlus_min_ET_eta_ProductChinese->SetMarkerColor (kBlue);                        gr_EEMinus_min_ET_eta_ProductChinese->SetMarkerColor (kBlue);
    gr_EEPlus_min_ET_eta_ProductChinese->SetLineWidth (1);                             gr_EEMinus_min_ET_eta_ProductChinese->SetLineWidth (1);
    gr_EEPlus_min_ET_eta_ProductChinese->SetLineColor (kBlue);                          gr_EEMinus_min_ET_eta_ProductChinese->SetLineColor (kBlue);


    //   ---- style (end) ----

    TLegend* legendChinese = new TLegend(0.1,0.7,0.48,0.9);
    legendChinese ->AddEntry(gr_EEMinus_ProductChinese,"EE-","lep");
    legendChinese ->AddEntry(gr_EEPlus_ProductChinese, "EE+","lep");


    TCanvas* ccRingChinese = new TCanvas ("ccRingChinese","",800,600);
    ccRingChinese->Divide(2,1);

    ccRingChinese->cd(1);

    gr_EEMinus_ProductChinese->Draw("APL");
    gr_EEPlus_ProductChinese->Draw("PL");
    legendChinese ->Draw();


    ccRingChinese->cd(2);

    gr_EEPlus_ProductChineseSpread->Draw("APL");

    gr_EEMinus_ProductChineseSpread->Draw("PL");


    gr_EEMinus_ProductChinese->GetYaxis()->SetTitle("Noise [GeV]");
    gr_EEMinus_ProductChinese->GetXaxis()->SetTitle("iRing");

    gr_EEPlus_ProductChineseSpread->GetYaxis()->SetTitle("Noise RMS [GeV]");
    gr_EEPlus_ProductChineseSpread->GetXaxis()->SetTitle("iRing");

    ccRingChinese->SaveAs("NoiseRingChinese_EE.png");
    ccRingChinese->SaveAs("NoiseRingChinese_EE.root");


    TCanvas* ccRingChinese2 = new TCanvas ("ccRingChinese2","",800,600);
    gr_EEMinus_ProductChinese->Draw("APL");
    gr_EEPlus_ProductChinese->Draw("PL");
    legendChinese ->Draw();
    ccRingChinese2->SaveAs("NoiseRingChinese2_EE.png");
    ccRingChinese2->SaveAs("NoiseRingChinese2_EE.root");



    TCanvas* ccRingChinese_eta = new TCanvas ("ccRingChinese_eta","",800,600);
    gr_EEMinus_eta_ProductChinese->Draw("APL");
    gr_EEPlus_eta_ProductChinese->Draw("PL");
    legendChinese ->Draw();

    gr_EEMinus_eta_ProductChinese->GetYaxis()->SetTitle("Noise [GeV]");
    gr_EEMinus_eta_ProductChinese->GetXaxis()->SetTitle("#eta");

    ccRingChinese_eta->SaveAs("NoiseRingChinese_eta_EE.png");
    ccRingChinese_eta->SaveAs("NoiseRingChinese_eta_EE.root");

    TF1 *f1Chinese = new TF1("f1Chinese","[0]+[1]*exp([2]*x)",1.5,3.1);
    f1Chinese->SetParameter(0,0.1);
    f1Chinese->SetParameter(1,0.);
    f1Chinese->SetParameter(2,0.1);
    //cout << "Fit to minus, average" << endl;
    //gr_EEMinus_eta_ProductChinese->Fit("f1Chinese","","",2.4,3.0);
    //cout << endl;
    //cout << "Fit to plus, average" << endl;
    //gr_EEPlus_eta_ProductChinese->Fit("expo","","",2.3,3.0);

    TCanvas* ccRingChinese_eta_minmax = new TCanvas ("ccRingChinese_eta_minmax","",800,600);

    gr_EEPlus_max_eta_ProductChinese->Draw("AP");

    gr_EEMinus_eta_ProductChinese->Draw("PL");
    gr_EEPlus_eta_ProductChinese->Draw("PL");

    gr_EEPlus_min_eta_ProductChinese->Draw("P");
    gr_EEMinus_min_eta_ProductChinese->Draw("P");
    gr_EEMinus_max_eta_ProductChinese->Draw("P");

    legendChinese ->Draw();

    gr_EEPlus_max_eta_ProductChinese->GetYaxis()->SetTitle("Noise [GeV]");
    gr_EEPlus_max_eta_ProductChinese->GetXaxis()->SetTitle("#eta");

    ccRingChinese_eta_minmax->SaveAs("NoiseRingChinese_eta_minmax_EE.png");
    ccRingChinese_eta_minmax->SaveAs("NoiseRingChinese_eta_minmax_EE.root");




    TLegend* legendChinese2 = new TLegend(0.1,0.7,0.48,0.9);
    legendChinese2->AddEntry(gr_EEMinus_ET_eta_ProductChinese,"average EE-","lep");
    legendChinese2->AddEntry(gr_EEPlus_ET_eta_ProductChinese, "average EE+","lep");
    legendChinese2->AddEntry(gr_EEMinus_min_ET_eta_ProductChinese, "min/max EE-","lep");
    legendChinese2->AddEntry(gr_EEPlus_max_ET_eta_ProductChinese, "min/max EE+","lep");



    TCanvas* ccRingChinese_ring_minmax = new TCanvas ("ccRingChinese_ring_minmax","",800,600);

    gr_ring_EEPlus_max_eta_ProductChinese->Draw("AP");

    gr_ring_EEMinus_eta_ProductChinese->Draw("PL");
    gr_ring_EEPlus_eta_ProductChinese->Draw("PL");

    gr_ring_EEPlus_min_eta_ProductChinese->Draw("P");
    gr_ring_EEMinus_min_eta_ProductChinese->Draw("P");
    gr_ring_EEMinus_max_eta_ProductChinese->Draw("P");

    legendChinese2->Draw();

    gr_ring_EEPlus_max_eta_ProductChinese->GetYaxis()->SetTitle("Noise [GeV]");
    gr_ring_EEPlus_max_eta_ProductChinese->GetXaxis()->SetTitle("ring");

    ccRingChinese_ring_minmax->SaveAs("NoiseRingChinese_ring_minmax_ring_EE.png");
    ccRingChinese_ring_minmax->SaveAs("NoiseRingChinese_ring_minmax_ring_EE.root");



    TCanvas* ccRingChinese_ET_eta = new TCanvas ("ccRingChinese_ET_eta","",800,600);
    gr_EEMinus_ET_eta_Product->Draw("APL");
    gr_EEPlus_ET_eta_Product->Draw("PL");
    legendChinese->Draw();

    gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
    //gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
    gr_EEMinus_ET_eta_Product->GetXaxis()->SetTitle("#eta");

    ccRingChinese_ET_eta->SaveAs("NoiseRingChinese_ET_eta_EE.png");
    ccRingChinese_ET_eta->SaveAs("NoiseRingChinese_ET_eta_EE.root");






    TCanvas* ccRingChinese_ET_eta_minmax = new TCanvas ("ccRingChinese_ET_eta_minmax","",800,600);

    gr_EEPlus_max_ET_eta_Product->Draw("AP");

    gr_EEMinus_ET_eta_Product->Draw("PL");
    gr_EEPlus_ET_eta_Product->Draw("PL");

    gr_EEPlus_min_ET_eta_Product->Draw("P");
    gr_EEMinus_min_ET_eta_Product->Draw("P");
    gr_EEMinus_max_ET_eta_Product->Draw("P");

    legendChinese2->Draw();

    gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
    //gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
    gr_EEPlus_max_ET_eta_Product->GetXaxis()->SetTitle("#eta");

    ccRingChinese_ET_eta_minmax->SaveAs("NoiseRingChinese_ET_eta_minmax_EE.png");
    ccRingChinese_ET_eta_minmax->SaveAs("NoiseRingChinese_ET_eta_minmax_EE.root");




    //Russian


        //   ---- style ----

        gr_EEPlus_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_ProductRussian->SetMarkerColor (kRed);                        gr_EEMinus_ProductRussian->SetMarkerColor (kRed);
        gr_EEPlus_ProductRussian->SetLineWidth (1);                             gr_EEMinus_ProductRussian->SetLineWidth (1);
        gr_EEPlus_ProductRussian->SetLineColor (kRed);                          gr_EEMinus_ProductRussian->SetLineColor (kRed);

        gr_EEPlus_ProductRussianSpread->SetMarkerSize  (1);                           gr_EEMinus_ProductRussianSpread->SetMarkerSize  (1);
        gr_EEPlus_ProductRussianSpread->SetMarkerStyle (24);                          gr_EEMinus_ProductRussianSpread->SetMarkerStyle (22);
        gr_EEPlus_ProductRussianSpread->SetMarkerColor (kRed);                        gr_EEMinus_ProductRussianSpread->SetMarkerColor (kRed);
        gr_EEPlus_ProductRussianSpread->SetLineWidth (1);                             gr_EEMinus_ProductRussianSpread->SetLineWidth (1);
        gr_EEPlus_ProductRussianSpread->SetLineColor (kRed);                          gr_EEMinus_ProductRussianSpread->SetLineColor (kRed);

        gr_EEPlus_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_eta_ProductRussian->SetMarkerColor (kRed);                        gr_EEMinus_eta_ProductRussian->SetMarkerColor (kRed);
        gr_EEPlus_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_eta_ProductRussian->SetLineColor (kRed);                          gr_EEMinus_eta_ProductRussian->SetLineColor (kRed);

        gr_ring_EEPlus_eta_ProductRussian->SetMarkerSize  (1);                           gr_ring_EEMinus_eta_ProductRussian->SetMarkerSize  (1);
        gr_ring_EEPlus_eta_ProductRussian->SetMarkerStyle (24);                          gr_ring_EEMinus_eta_ProductRussian->SetMarkerStyle (22);
        gr_ring_EEPlus_eta_ProductRussian->SetMarkerColor (kRed);                        gr_ring_EEMinus_eta_ProductRussian->SetMarkerColor (kRed);
        gr_ring_EEPlus_eta_ProductRussian->SetLineWidth (1);                             gr_ring_EEMinus_eta_ProductRussian->SetLineWidth (1);
        gr_ring_EEPlus_eta_ProductRussian->SetLineColor (kRed);                          gr_ring_EEMinus_eta_ProductRussian->SetLineColor (kRed);

        gr_EEPlus_ET_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_ET_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_ET_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_ET_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_ET_eta_ProductRussian->SetMarkerColor (kRed);                        gr_EEMinus_ET_eta_ProductRussian->SetMarkerColor (kRed);
        gr_EEPlus_ET_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_ET_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_ET_eta_ProductRussian->SetLineColor (kRed);                          gr_EEMinus_ET_eta_ProductRussian->SetLineColor (kRed);


        gr_EEPlus_max_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_max_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_max_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_max_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_max_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_EEMinus_max_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_EEPlus_max_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_max_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_max_eta_ProductRussian->SetLineColor (kBlue);                          gr_EEMinus_max_eta_ProductRussian->SetLineColor (kBlue);

        gr_EEPlus_min_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_min_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_min_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_min_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_min_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_EEMinus_min_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_EEPlus_min_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_min_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_min_eta_ProductRussian->SetLineColor (kBlue);                          gr_EEMinus_min_eta_ProductRussian->SetLineColor (kBlue);


        gr_ring_EEPlus_max_eta_ProductRussian->SetMarkerSize  (1);                           gr_ring_EEMinus_max_eta_ProductRussian->SetMarkerSize  (1);
        gr_ring_EEPlus_max_eta_ProductRussian->SetMarkerStyle (24);                          gr_ring_EEMinus_max_eta_ProductRussian->SetMarkerStyle (22);
        gr_ring_EEPlus_max_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_ring_EEMinus_max_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_ring_EEPlus_max_eta_ProductRussian->SetLineWidth (1);                             gr_ring_EEMinus_max_eta_ProductRussian->SetLineWidth (1);
        gr_ring_EEPlus_max_eta_ProductRussian->SetLineColor (kBlue);                          gr_ring_EEMinus_max_eta_ProductRussian->SetLineColor (kBlue);

        gr_ring_EEPlus_min_eta_ProductRussian->SetMarkerSize  (1);                           gr_ring_EEMinus_min_eta_ProductRussian->SetMarkerSize  (1);
        gr_ring_EEPlus_min_eta_ProductRussian->SetMarkerStyle (24);                          gr_ring_EEMinus_min_eta_ProductRussian->SetMarkerStyle (22);
        gr_ring_EEPlus_min_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_ring_EEMinus_min_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_ring_EEPlus_min_eta_ProductRussian->SetLineWidth (1);                             gr_ring_EEMinus_min_eta_ProductRussian->SetLineWidth (1);
        gr_ring_EEPlus_min_eta_ProductRussian->SetLineColor (kBlue);                          gr_ring_EEMinus_min_eta_ProductRussian->SetLineColor (kBlue);


        gr_EEPlus_max_ET_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_max_ET_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_max_ET_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_max_ET_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_max_ET_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_EEMinus_max_ET_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_EEPlus_max_ET_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_max_ET_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_max_ET_eta_ProductRussian->SetLineColor (kBlue);                          gr_EEMinus_max_ET_eta_ProductRussian->SetLineColor (kBlue);

        gr_EEPlus_min_ET_eta_ProductRussian->SetMarkerSize  (1);                           gr_EEMinus_min_ET_eta_ProductRussian->SetMarkerSize  (1);
        gr_EEPlus_min_ET_eta_ProductRussian->SetMarkerStyle (24);                          gr_EEMinus_min_ET_eta_ProductRussian->SetMarkerStyle (22);
        gr_EEPlus_min_ET_eta_ProductRussian->SetMarkerColor (kBlue);                        gr_EEMinus_min_ET_eta_ProductRussian->SetMarkerColor (kBlue);
        gr_EEPlus_min_ET_eta_ProductRussian->SetLineWidth (1);                             gr_EEMinus_min_ET_eta_ProductRussian->SetLineWidth (1);
        gr_EEPlus_min_ET_eta_ProductRussian->SetLineColor (kBlue);                          gr_EEMinus_min_ET_eta_ProductRussian->SetLineColor (kBlue);


        //   ---- style (end) ----

        TLegend* legendRussian = new TLegend(0.1,0.7,0.48,0.9);
        legendRussian ->AddEntry(gr_EEMinus_ProductRussian,"EE-","lep");
        legendRussian ->AddEntry(gr_EEPlus_ProductRussian, "EE+","lep");


        TCanvas* ccRingRussian = new TCanvas ("ccRingRussian","",800,600);
        ccRingRussian->Divide(2,1);

        ccRingRussian->cd(1);

        gr_EEMinus_ProductRussian->Draw("APL");
        gr_EEPlus_ProductRussian->Draw("PL");
        legendRussian ->Draw();


        ccRingRussian->cd(2);

        gr_EEPlus_ProductRussianSpread->Draw("APL");

        gr_EEMinus_ProductRussianSpread->Draw("PL");


        gr_EEMinus_ProductRussian->GetYaxis()->SetTitle("Noise [GeV]");
        gr_EEMinus_ProductRussian->GetXaxis()->SetTitle("iRing");

        gr_EEPlus_ProductRussianSpread->GetYaxis()->SetTitle("Noise RMS [GeV]");
        gr_EEPlus_ProductRussianSpread->GetXaxis()->SetTitle("iRing");

        ccRingRussian->SaveAs("NoiseRingRussian_EE.png");
        ccRingRussian->SaveAs("NoiseRingRussian_EE.root");


        TCanvas* ccRingRussian2 = new TCanvas ("ccRingRussian2","",800,600);
        gr_EEMinus_ProductRussian->Draw("APL");
        gr_EEPlus_ProductRussian->Draw("PL");
        legendRussian ->Draw();
        ccRingRussian2->SaveAs("NoiseRingRussian2_EE.png");
        ccRingRussian2->SaveAs("NoiseRingRussian2_EE.root");



        TCanvas* ccRingRussian_eta = new TCanvas ("ccRingRussian_eta","",800,600);
        gr_EEMinus_eta_ProductRussian->Draw("APL");
        gr_EEPlus_eta_ProductRussian->Draw("PL");
        legendRussian ->Draw();

        gr_EEMinus_eta_ProductRussian->GetYaxis()->SetTitle("Noise [GeV]");
        gr_EEMinus_eta_ProductRussian->GetXaxis()->SetTitle("#eta");

        ccRingRussian_eta->SaveAs("NoiseRingRussian_eta_EE.png");
        ccRingRussian_eta->SaveAs("NoiseRingRussian_eta_EE.root");

        TF1 *f1Russian = new TF1("f1Russian","[0]+[1]*exp([2]*x)",1.5,3.1);
        f1Russian->SetParameter(0,0.1);
        f1Russian->SetParameter(1,0.);
        f1Russian->SetParameter(2,0.1);
        //cout << "Fit to minus, average" << endl;
        //gr_EEMinus_eta_ProductRussian->Fit("f1Russian","","",2.4,3.0);
        //cout << endl;
        //cout << "Fit to plus, average" << endl;
        //gr_EEPlus_eta_ProductRussian->Fit("expo","","",2.3,3.0);

        TCanvas* ccRingRussian_eta_minmax = new TCanvas ("ccRingRussian_eta_minmax","",800,600);

        gr_EEPlus_max_eta_ProductRussian->Draw("AP");

        gr_EEMinus_eta_ProductRussian->Draw("PL");
        gr_EEPlus_eta_ProductRussian->Draw("PL");

        gr_EEPlus_min_eta_ProductRussian->Draw("P");
        gr_EEMinus_min_eta_ProductRussian->Draw("P");
        gr_EEMinus_max_eta_ProductRussian->Draw("P");

        legendRussian ->Draw();

        gr_EEPlus_max_eta_ProductRussian->GetYaxis()->SetTitle("Noise [GeV]");
        gr_EEPlus_max_eta_ProductRussian->GetXaxis()->SetTitle("#eta");

        ccRingRussian_eta_minmax->SaveAs("NoiseRingRussian_eta_minmax_EE.png");
        ccRingRussian_eta_minmax->SaveAs("NoiseRingRussian_eta_minmax_EE.root");




        TLegend* legendRussian2 = new TLegend(0.1,0.7,0.48,0.9);
        legendRussian2->AddEntry(gr_EEMinus_ET_eta_ProductRussian,"average EE-","lep");
        legendRussian2->AddEntry(gr_EEPlus_ET_eta_ProductRussian, "average EE+","lep");
        legendRussian2->AddEntry(gr_EEMinus_min_ET_eta_ProductRussian, "min/max EE-","lep");
        legendRussian2->AddEntry(gr_EEPlus_max_ET_eta_ProductRussian, "min/max EE+","lep");



        TCanvas* ccRingRussian_ring_minmax = new TCanvas ("ccRingRussian_ring_minmax","",800,600);

        gr_ring_EEPlus_max_eta_ProductRussian->Draw("AP");

        gr_ring_EEMinus_eta_ProductRussian->Draw("PL");
        gr_ring_EEPlus_eta_ProductRussian->Draw("PL");

        gr_ring_EEPlus_min_eta_ProductRussian->Draw("P");
        gr_ring_EEMinus_min_eta_ProductRussian->Draw("P");
        gr_ring_EEMinus_max_eta_ProductRussian->Draw("P");

        legendRussian2->Draw();

        gr_ring_EEPlus_max_eta_ProductRussian->GetYaxis()->SetTitle("Noise [GeV]");
        gr_ring_EEPlus_max_eta_ProductRussian->GetXaxis()->SetTitle("ring");

        ccRingRussian_ring_minmax->SaveAs("NoiseRingRussian_ring_minmax_ring_EE.png");
        ccRingRussian_ring_minmax->SaveAs("NoiseRingRussian_ring_minmax_ring_EE.root");



        TCanvas* ccRingRussian_ET_eta = new TCanvas ("ccRingRussian_ET_eta","",800,600);
        gr_EEMinus_ET_eta_Product->Draw("APL");
        gr_EEPlus_ET_eta_Product->Draw("PL");
        legendRussian->Draw();

        gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
        //gr_EEMinus_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
        gr_EEMinus_ET_eta_Product->GetXaxis()->SetTitle("#eta");

        ccRingRussian_ET_eta->SaveAs("NoiseRingRussian_ET_eta_EE.png");
        ccRingRussian_ET_eta->SaveAs("NoiseRingRussian_ET_eta_EE.root");






        TCanvas* ccRingRussian_ET_eta_minmax = new TCanvas ("ccRingRussian_ET_eta_minmax","",800,600);

        gr_EEPlus_max_ET_eta_Product->Draw("AP");

        gr_EEMinus_ET_eta_Product->Draw("PL");
        gr_EEPlus_ET_eta_Product->Draw("PL");

        gr_EEPlus_min_ET_eta_Product->Draw("P");
        gr_EEMinus_min_ET_eta_Product->Draw("P");
        gr_EEMinus_max_ET_eta_Product->Draw("P");

        legendRussian2->Draw();

        gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Transverse Noise [GeV]");
        //gr_EEPlus_max_ET_eta_Product->GetYaxis()->SetTitle("Noise [GeV]");
        gr_EEPlus_max_ET_eta_Product->GetXaxis()->SetTitle("#eta");

        ccRingRussian_ET_eta_minmax->SaveAs("NoiseRingRussian_ET_eta_minmax_EE.png");
        ccRingRussian_ET_eta_minmax->SaveAs("NoiseRingRussian_ET_eta_minmax_EE.root");



}


int main(int argc, char** argvv)
{
    if (argc < 2) {
        std::cout<< "Please provide the TimeStart,TimeEnd and the Tags";
    }


    /*
     *Alpha = 1. is Chinese SIC and  1.16 is Russian BTCP
     * Make plots with different crystal country
     *
     * */


    double TimeStart = 1;
    double TimeEnd = 100000;
    double laserTime = 1530432369;
        if (argc > 2)
        {
            TimeStart = std::stod(argvv[1]);
            TimeEnd = std::stod(argvv[2]);
            laserTime = std::stod(argvv[8]);


        }


    std::string tag_EcalIntercalibConstants = "EcalIntercalibConstants_2017_V1_Bon_mc";
    std::string tag_EcalLaserAlphas = "EcalLaserAlphas_2017_mc";
    std::string tag_EcalPedestals = "EcalPedestals_2017extrap_25fb_mc";
    std::string tag_EcalADCToGeVConstant = "EcalADCToGeVConstant_2017_V1_Bon_mc";
    std::string tag_EcalLaserAPDPNRatios = "EcalLaserAPDPNRatios_2017_mc";
    if(argc > 7)
    {
        tag_EcalADCToGeVConstant = argvv[3];
        tag_EcalIntercalibConstants = argvv[4];
        tag_EcalLaserAlphas = argvv[5];
        tag_EcalPedestals = argvv[6];
        tag_EcalLaserAPDPNRatios = argvv[7];


    }
    cond::NoiseDumper<EcalIntercalibConstants> IC("EcalIntercalibConstants");
    char* argv[] =  {(char *) "hi!",(char *) "-O",(char *)"EcalIntercalibConstants",(char *) "-t",(char *) "EcalIntercalibConstants_2012ABCD_offline"};
    IC.run(5,argv);
    IC.TimeLimitStart = TimeStart;
    IC.TimeLimit = TimeEnd;
    IC.MyTime = laserTime;
    IC.executeClone(tag_EcalIntercalibConstants,IC.sessionMain);

    // char* argv2[] =  {(char *) "hi!",(char *) "-O",(char *)"EcalLaserAlpha",(char *) "-t",(char *) "EcalLaserAlphas_Legacy2016_v1"};
    cond::NoiseDumper<EcalLaserAlphas> Alpha("EcalLaserAlphas");
    //Alpha.run(5,argv2);
    Alpha.connect = "frontier://FrontierProd/CMS_CONDITIONS";
    Alpha.TimeLimitStart = TimeStart;
    Alpha.TimeLimit = TimeEnd;
    Alpha.MyTime = laserTime;
    Alpha.executeClone(tag_EcalLaserAlphas,IC.sessionMain);


    cond::NoiseDumper<EcalPedestals> Ped("EcalPedestals");
    Ped.connect = "frontier://FrontierProd/CMS_CONDITIONS";
    Ped.TimeLimitStart = TimeStart;
    Ped.TimeLimit = TimeEnd;
    Ped.MyTime = laserTime;
    Ped.executeClone(tag_EcalPedestals,IC.sessionMain);

    cond::NoiseDumper<EcalADCToGeVConstant> ADCtoGEV("EcalADCToGeVConstant");
    ADCtoGEV.connect = "frontier://FrontierProd/CMS_CONDITIONS";
    ADCtoGEV.TimeLimitStart = TimeStart;
    ADCtoGEV.TimeLimit = TimeEnd;
    ADCtoGEV.MyTime = laserTime;
    ADCtoGEV.executeClone(tag_EcalADCToGeVConstant,IC.sessionMain);

    cond::NoiseLaserDumper laser;
    //    char* argv2[] =  {(char *) "lava_db2txt",(char *) "-t", (char *) "EcalLaserAPDPNRatios_prompt_v2",(char *) "-c",(char *) "frontier://FrontierProd/CMS_CONDITIONS",(char *) "-b",(char *) "6130819125003419648",(char *) "-e",(char *) "6130823351251238912",(char *) "--dump",(char *) "APDPN_Data2015fromData.txt"};
    //    laser.session = IC.sessionMain;
    //    laser.run(11,argv2);
    laser.TimeLimitStart = TimeStart;
    laser.TimeLimit = TimeEnd;
    laser.MyTime = laserTime;

    laser.executeClone(tag_EcalLaserAPDPNRatios,IC.sessionMain);


    std::cout << (time_t) IC.a.IOVBegin << (time_t) IC.a.IOVEnd << std::endl;
    plotIC(IC);
    plotAlpha(Alpha);
    plotPedestals(Ped);
    plotLaser(laser);
    productOfTag(ADCtoGEV);





    return 0;
}

