{
    TFile * f = new TFile("plots_testGEMGeo_clean.root");
     for(unsigned int iE =  0 ; iE < 8; ++iE){
         Plotter * p = new Plotter;
       for(unsigned int iX =  1 ; iX < 7; ++iX){
           TH1 * h;
           f->GetObject(TString::Format("deltaXoSig_%u_%u",iE,iX),h);
           p->addHistLine(h,TString::Format("%u",iX),-1,1,1);
       }
       p->normalize();
       p->setYTitle(TString::Format("%u",iE));
       p->draw();
   }
}


{
    TFile * f = new TFile("plots_testGEMGeo_clean.root");
     for(unsigned int iX2 =  0 ; iX2 < 6; ++iX2){
         Plotter * p = new Plotter;
       for(unsigned int iX =  1 ; iX < 7; ++iX){
           TH1 * h;
           f->GetObject(TString::Format("deltaXoSig_a_%u_x_%u",iX2,iX),h);
           p->addHistLine(h,TString::Format("%u",iX),-1,1,1);
       }
       p->normalize();
       p->setYTitle(TString::Format("%u",iX2));
       p->draw();
   }
}


{
    TFile * f = new TFile("plots_testGEMGeo_clean.root");
    TString vars[] = {"deltaY","deltaYoSig",""};
    
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iE =  0 ; iE < 8; ++iE){
         TH1 * h;
         f->GetObject(TString::Format("%s_eta_%u",vars[iV].Data(),iE),h);
         p->addHistLine(h,TString::Format("#eta part. %u",iE),-1,1,2);         
         }  
         p->normalize();
         p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));               
    }    
}


{
    TFile * f = new TFile("plots_testGEMGeo_clean.root");
    TString vars[] = {"deltaY","deltaYoSig",""};
    
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iE =  0 ; iE < 8; ++iE){
         TH1 * h;
         f->GetObject(TString::Format("%s_eta_%u",vars[iV].Data(),iE),h);
         p->addHistLine(h,TString::Format("#eta part. %u",iE),-1,1,2);         
         }  
         p->normalize();
         p->draw(true,TString::Format("plots/%s.pdf",vars[iV].Data()));               
    }    
}


{
    // TString type = "";
    TString type = "_increaseCSC";
    
    TFile * f = new TFile(TString::Format("plots_testGEMGeo%s_clean.root",type.Data()));
    TString vars[] = {"deltaX","deltaXoSig",""};
    TString extra = "";
        // TString extra = "_const";
// TString extra =     "_const_oneStrip";
    double xBins[] = {-22,-16.5,-11,-5.5,0,5.5,11,16.5,22};
    
    for(unsigned int iV =  0 ; vars[iV][0]; ++iV){
         Plotter * p = new Plotter;
     for(unsigned int iX =  0 ; iX < 8; ++iX){
         TH1 * h;
         TString name = TString::Format("%s%s_x_%.0f_%.0f", vars[iV].Data(),extra.Data(),xBins[iX]*10,xBins[iX+1]*10);
         f->GetObject(name,h);
         TString title = TString::Format("%.1f < x_{g} < %.1f",xBins[iX],xBins[iX+1]);
         p->addHistLine(h,title,-1,1,2);         
         }  
         p->normalize();
         p->draw(true,TString::Format("plots/%s%s%s.pdf",vars[iV].Data(),extra.Data(),type.Data()));               
    }    
}

{
    gStyle->SetPaintTextFormat(".2f");
    noisemap_num->Divide(noisemap_denom);
    noisemap_num->Draw("COLZ");
    
    TH1 * h = new TH1F("noiseComp",";noise rate per strip;number of strips",50,0,.00005);    
    for(unsigned int iX= 1; iX <= noisemap_num->GetNbinsX(); ++iX ){
        for(unsigned int iY= 1; iY <= noisemap_num->GetNbinsY(); ++iY ){
            if(occmap->GetBinContent(iX,iY) > 0)
            h->Fill(noisemap_num->GetBinContent(iX,iY));
        }
    }
    new TCanvas();
    h->Draw();
    
    TFile * f = new TFile("plots_testGEMNoiseAndEff_increaseCSC_clean.root");
    
    Plotter * p = new Plotter;
for(unsigned int iE =  0 ; iE < 8; ++iE){
    TH1 * h;
    TString name = TString::Format("noiseClusterCount_eta_%u",iE);
    f->GetObject(name,h);
    p->addHistLine(h,TString::Format("#eta part. %u",iE),-1,1,2);         
    }  
    p->normalize();
    p->draw(true,TString::Format("plots/noiseClusterCount.pdf"));               
    
    
}

{
    gStyle->SetPaintTextFormat(".2f");
    effmap_num->Divide(effmap_denom);
    effmap_num->Draw("COLZTEXT");
    
    new TCanvas();
    effmap_num_int->Divide(effmap_denom_int);
    effmap_num_int->Draw();
    
}



    
