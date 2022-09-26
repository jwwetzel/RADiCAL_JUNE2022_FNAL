
/*Function to get integral*/
/*These signals start around timeslice 200*/
/*Average the first 50 to estimate 'zero' signal*/
double getIntegral(Float_t f_amp[1024])
{
    double d_min = 9999999;
    double d_ped=0;
    double d_integral = 0.0;
    
    /*Pedestal value calculation*/
    /*Average the first 50 timeslices*/
    for(int i=0; i!=50; ++i)
    {
        d_ped+= f_amp[i]/50.0;
    }
    
    /*Plot the subtracted waveform*/
    /*Calculate integral*/
    TH1F *baseLine = new TH1F("bL","baseLine",50,0,5);
    double d_t[1024], d_amp[1024];
    for(int i = 0; i!=1024; ++i)
    {
        d_t[i] = i*0.2;
        d_amp[i] = f_amp[i]-d_ped;
        if (i<50) baseLine->Fill(d_amp[i]);
        d_integral += f_amp[i]-d_ped;
    }
    TGraph *waveForm = new TGraph(1024,d_t,d_amp);
    TFile *temp = new TFile("test.root","RECREATE");
    baseLine->Write();
    waveForm->Write();
    return d_integral;
}

double GetPedestal(Float_t f_amp[1024])
{
    double d_ped=0;
    /*Pedestal value calculation*/
    /*Average the first 50 timeslices*/
    for(int i=0; i!=50; ++i)
    {
        d_ped+= f_amp[i]/50.0;
    }
    return d_ped;
}

double GetMaxAmplitude(Float_t f_amp[1024])
{
  double min = 9999999;
  double sx=0;
  for(int i=0; i!=1024; ++i)
  {
    if( f_amp[i]<min ) min = f_amp[i];
    if(i<50)
    {
     sx+= f_amp[i];
    }
  }
  return sx/50. - min;
}

void AnalyzeRAD()
{
    ROOT::EnableImplicitMT();
    
    struct {
      Float_t RAD_time[2][1024];
      Float_t RAD_amp[2][9][1024];
      Float_t BTL_amp[2][9][1024];
      Int_t   RAD_triggerno;
        
      Float_t PIX_xSlope;
      Float_t PIX_ySlope;
      Float_t PIX_xIntercept;
      Float_t PIX_yIntercept;
      Float_t PIX_chi2;
    } event;
    
    TChain *dataTree = new TChain("tree");
    dataTree->Add("../Data/RUN61733.root");
    dataTree->Add("../Data/RUN61734.root");
    dataTree->Add("../Data/RUN61735.root");
    dataTree->Add("../Data/RUN61736.root");
    dataTree->Add("../Data/RUN61737.root");
    dataTree->Add("../Data/RUN61738.root");
    dataTree->Add("../Data/RUN61739.root");
    
    dataTree->SetBranchAddress("RAD_x",event.RAD_time); //
    dataTree->SetBranchAddress("RAD_y",event.RAD_amp); //
    dataTree->SetBranchAddress("BTL_y",event.BTL_amp); //
    dataTree->SetBranchAddress("RAD_triggerno",&event.RAD_triggerno); //
    dataTree->SetBranchAddress("PIX_xIntercept",&event.PIX_xIntercept); //
    dataTree->SetBranchAddress("PIX_yIntercept",&event.PIX_yIntercept); //
    dataTree->SetBranchAddress("PIX_xSlope",&event.PIX_xSlope); //
    dataTree->SetBranchAddress("PIX_ySlope",&event.PIX_ySlope); //
    dataTree->SetBranchAddress("PIX_chi2",&event.PIX_chi2); //
    
    int nev = dataTree->GetEntries();
    
    TH1F *hist = new TH1F("name","title",10000,0,10000);
    
    TH1F *h_rmsHists[8];
    TProfile *tp_Profiles[2][8];
    TH1F *pedHists[2][8];

    char c_name[50];
    for (int i = 0; i!=8; ++i)
    {
        sprintf(c_name,"rms_hist_Low_Gain_ch_%d",i+1);
        h_rmsHists[i] = new TH1F(c_name,c_name,100,0,5);
        
        sprintf(c_name,"tprofile_Low_Gain_ch_%d",i+1);
        tp_Profiles[0][i] = new TProfile(c_name,c_name,1024,0,204.8);
        
        sprintf(c_name,"tprofile_High_Gain_ch_%d",i+1);
        tp_Profiles[1][i] = new TProfile(c_name,c_name,1024,0,204.8);
        
        sprintf(c_name,"Peds_Low_Gain_ch_%d",i+1);
        pedHists[0][i] = new TH1F(c_name,c_name,500,0,500);
        
        sprintf(c_name,"Peds_High_Gain_ch_%d",i+1);
        pedHists[1][i] = new TH1F(c_name,c_name,500,0,500);
    }
    
    for(int iev=0; iev<nev; ++iev)
    {
        if (iev%1000 == 1) cout << "Processing even: " << iev << endl;
        dataTree->GetEntry(iev);  // Get the iev'th event
        
//        if (abs(GetMaxAmplitude(event.RAD_amp[0][0])) > 50) continue;
        
        double d_eventPed[2][8];
        d_eventPed[0][0] = GetPedestal(event.BTL_amp[1][4]); //T1Low
        d_eventPed[1][0] = GetPedestal(event.BTL_amp[1][6]);
        for (int i = 1; i!=8; ++i)
        {
            d_eventPed[0][i] = GetPedestal(event.RAD_amp[0][i]);
        }

        for (int i = 1; i!=8; ++i)
        {
            d_eventPed[1][i] = GetPedestal(event.RAD_amp[1][i]);
        }
        
        for(int i = 0; i!= 8; ++i)
        {
            for(int j = 0; j != 50; ++j)
            {
                if(i == 0)
                {
                    h_rmsHists[0]->Fill(sqrt((event.BTL_amp[1][4][j]-d_eventPed[0][0])*(event.BTL_amp[1][4][j]-d_eventPed[0][0])/50.));
                }
                else
                {
                    h_rmsHists[i]->Fill(sqrt((event.RAD_amp[0][i][j]-d_eventPed[0][i])*(event.RAD_amp[0][i][j]-d_eventPed[0][i])/50.));
                }
            }
        }

        double t[1024], amp[2][8][1024];
        
        for(int i = 0; i != 8; ++i)
        {
            pedHists[0][i]->Fill(d_eventPed[0][i]);
            pedHists[1][i]->Fill(d_eventPed[1][i]);
        }
        
    
        for (int i = 0; i != 8; ++i)
        {
            for (int j = 0; j != 1024; ++j)
            {
                t[j] = j*0.2;
                if(i == 0)
                {
                    amp[0][0][j] = -1*(event.BTL_amp[1][4][j]-d_eventPed[0][0]);
                    amp[1][0][j] = -1*(event.BTL_amp[1][6][j]-d_eventPed[1][0]);
                }
                else
                {
                    amp[0][i][j] = -1*(event.RAD_amp[0][i][j]-d_eventPed[0][i]);
                }
                amp[1][i][j] = -1*(event.RAD_amp[1][i][j]-d_eventPed[1][i]);
                tp_Profiles[0][i]->Fill(t[j],amp[0][i][j]);
                tp_Profiles[1][i]->Fill(t[j],amp[1][i][j]);
            }
        }
//        TGraph *waveGraph = new TGraph(1024,t[i],amp[i]);
//        waveGraph->Draw("ACSame");
    }
        
    tp_Profiles[0][0]->SetTitle("Average Waveform Low Gain U1");
    tp_Profiles[0][1]->SetTitle("Average Waveform Low Gain U2");
    tp_Profiles[0][2]->SetTitle("Average Waveform Low Gain U3");
    tp_Profiles[0][3]->SetTitle("Average Waveform Low Gain U4");
    tp_Profiles[0][4]->SetTitle("Average Waveform Low Gain D1");
    tp_Profiles[0][5]->SetTitle("Average Waveform Low Gain D2");
    tp_Profiles[0][6]->SetTitle("Average Waveform Low Gain D3");
    tp_Profiles[0][7]->SetTitle("Average Waveform Low Gain D4");
    
    tp_Profiles[1][0]->SetTitle("Average Waveform High Gain U1");
    tp_Profiles[1][1]->SetTitle("Average Waveform High Gain U2");
    tp_Profiles[1][2]->SetTitle("Average Waveform High Gain U3");
    tp_Profiles[1][3]->SetTitle("Average Waveform High Gain U4");
    tp_Profiles[1][4]->SetTitle("Average Waveform High Gain D1");
    tp_Profiles[1][5]->SetTitle("Average Waveform High Gain D2");
    tp_Profiles[1][6]->SetTitle("Average Waveform High Gain D3");
    tp_Profiles[1][7]->SetTitle("Average Waveform High Gain D4");
    
    TFile *outputFile = new TFile("output.root","RECREATE");
    
    TCanvas *canvas_LowGainProfile = new TCanvas("lowGainProfileCanvas","lowGainProfileCanvas",1600,800);
    canvas_LowGainProfile->Divide(4,2);
    
    TCanvas *canvas_HighGainProfile = new TCanvas("highGainProfileCanvas","highGainProfileCanvas",1600,800);
    canvas_HighGainProfile->Divide(4,2);
    
    for(int i = 0; i != 8; ++i)
    {
        canvas_LowGainProfile->cd(i+1);
        tp_Profiles[0][i]->Draw();
        tp_Profiles[0][i]->GetYaxis()->SetRangeUser(-20,100);
        tp_Profiles[0][i]->GetYaxis()->SetTitle("Average Amplitude (mV)");
        tp_Profiles[0][i]->GetXaxis()->SetTitle("Time (ns)");
        h_rmsHists[i]->Write();
        tp_Profiles[0][i]->Write();
    }
    
    for(int i = 0; i != 8; ++i)
    {
        canvas_HighGainProfile->cd(i+1);
        tp_Profiles[1][i]->Draw();
        tp_Profiles[1][i]->GetYaxis()->SetRangeUser(-60,120);
        tp_Profiles[1][i]->GetYaxis()->SetTitle("Average Amplitude (mV)");
        tp_Profiles[1][i]->GetXaxis()->SetTitle("Time (ns)");
        tp_Profiles[1][i]->Write();
        
        pedHists[0][i]->Write();
        pedHists[1][i]->Write();
    }
    
    canvas_LowGainProfile->Write();
    canvas_HighGainProfile->Write();
    
    
    outputFile->Close();
}
