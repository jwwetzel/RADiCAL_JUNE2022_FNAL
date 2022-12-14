/*Get the Pedestal*/
/*First loop estimates the pedestal*/
/*Second loop excludes anomalous spikes in signal*/
/*The pedestal is needed to estimate anomalos spikes*/
/*The spikes are a 'feature' of CAEN calibration*/
double GetPedestal(Float_t f_amp[1024], int i_numSamples, double d_rmsValue, TH1F* h_rmsHistNoCut, TH1F* h_rmsHistWithCut, TH1F* h_pedHist, double d_integral, double d_pedSubbedAmps)
{
    double d_ped        =   0.0;
    double d_newPed     =   0.0;
    float f_rmsCuts     =   0.0;
    float f_rmsNoCuts   =   0.0;
    
    /*Pedestal value calculation*/
    /*Average the first N timeslices*/
    for(int i=25; i!=i_numSamples+25; ++i)
    {
        d_ped += f_amp[i]/float(i_numSamples);
    }
    
    double tslice[1024];
    double fslice[1024];
    
    for(int i = 0; i != 1024; ++i)
    {
        tslice[i] = 0.2*i;
        fslice[i] = f_amp[i];
    }
    
    /*Estimate the pedestal without the janky spikes*/
    /*If 100 samples, but 10 are bad, need to average 90*/
    /*Then sampleCounter keeps track of how many make it*/
    int i_sampleCounter = 0;
    for(int i=25; i!=i_numSamples+25; ++i)
    {
        if (abs(f_amp[i] - d_ped) < 5)
        {
            i_sampleCounter += 1;
            d_newPed += f_amp[i];
        }
    }
    if (i_sampleCounter == 0) i_sampleCounter = 1;
    d_newPed = d_newPed/i_sampleCounter;
    i_sampleCounter = 0;
    for (int i=25; i != i_numSamples+25; ++i)
    {
        h_pedHist->Fill(f_amp[i] - d_ped);
        f_rmsNoCuts  += (f_amp[i]-d_ped)*(f_amp[i]-d_ped);
        if (abs(f_amp[i] - d_ped) < 5)
        {
            f_rmsCuts += (f_amp[i]-d_newPed)*(f_amp[i]-d_newPed);
            i_sampleCounter+=1;
        }
    }
    
    if (i_sampleCounter == 0) i_sampleCounter = 1;
    h_rmsHistNoCut->Fill(sqrt(f_rmsNoCuts/i_numSamples));
    h_rmsHistWithCut->Fill(sqrt(f_rmsCuts/i_sampleCounter));
    
    d_rmsValue = sqrt(f_rmsCuts/i_sampleCounter);
    
    return d_newPed;
}

/*This estimates the maximum amplitude less the pedestal*/
/*This pedestal is adjusted for the spikes in random time slices*/
double GetMaxAmplitude(Float_t f_amp[1024])
{
    double min  =   9999999;
    double sx   =   0.0;
    double sxx  =   0.0;
 
    /*This loop estimates pedestal (sx)*/
    for(int i=0; i!=1024; ++i)
    {
        if( f_amp[i]<min ) min = f_amp[i];
        if(i>25 && i<125)
        {
         sx+= f_amp[i];
        }
    }
    
    sx = sx/100.;
    
    /*This calculates pedestal excluding spikes */
    int counter = 0;
    for(int i=5; i!=105; ++i)
    {
        if( abs(f_amp[i]-sx) < 6 )
        {
            counter +=1;
            sxx+=f_amp[i];
        }
    }
    
    return sxx/counter - min;
}


/* This is the meat and potatoes of the analysis.*/
void analyzeRADData()
{
    /* Set ROOT styles */
//    gROOT->Reset();
//    gROOT->SetStyle("Plain");
    gStyle->SetPalette(57);
//    gStyle->SetOptTitle(0);
    gStyle->SetOptLogz();

//    int style_label_font=42;
//    gStyle->SetLabelFont(style_label_font,"xyz");
//    gStyle->SetLabelSize(0.045,"xyz");
//    gStyle->SetLabelOffset(0.015,"xyz");
//    gStyle->SetStatFont(style_label_font);
//    gStyle->SetTitleFont(style_label_font,"xyz"); // axis titles
//    gStyle->SetTitleFont(style_label_font,"h"); // histogram title
//    gStyle->SetTitleSize(0.05,"xyz"); // axis titles
//    gStyle->SetTitleSize(0.05,"h"); // histogram title
//    gStyle->SetTitleOffset(1.1,"x");
//    gStyle->SetTitleOffset(1.2,"y");
//    gStyle->SetStripDecimals(kFALSE); // if we have 1.5 do not set 1.0 -> 1
//    gStyle->SetTitleX(0.12); // spot where histogram title goes
//    gStyle->SetTitleW(0.78); // width computed so that title is centered
//    TGaxis::SetMaxDigits(2); // restrict the number of digits in labels

    gStyle->SetPadTopMargin(0.08);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.14);
    
    /* Define a struct to hold related data, call it event*/
    /* To access a member, use event.RAD_x for example */
    /* Or RAD_x[1][2]*/
    struct {
        Float_t RAD_x[2][1024];         // RAD_x[down/up][timeslice]
        Float_t RAD_y[2][2][4][1024];   // RAD_y[low/high][down/up][NW/NE/SW/SE][amplitude]
        Float_t PbG_y[1024];            // PbG_y[amplitude] lead glass calorimeter
        Float_t MCP_y[1024];            // MCP_y[amplitude]

        Int_t   RAD_triggerno;

        Float_t PIX_xSlope;
        Float_t PIX_ySlope;
        Float_t PIX_xIntercept;
        Float_t PIX_yIntercept;
        Float_t PIX_chi2;
    } event;
    
    TFile *inFile   = new TFile("RAD_Data_61733-61739.root","open");
    TTree *dataTree = (TTree*)inFile->Get("radTree");
    
    dataTree->SetBranchAddress("RAD_x",&event.RAD_x);                   //
    dataTree->SetBranchAddress("RAD_y",&event.RAD_y);                   //
    dataTree->SetBranchAddress("PbG_y",&event.PbG_y);                   //
    dataTree->SetBranchAddress("MCP_y",&event.MCP_y);                   //
    dataTree->SetBranchAddress("RAD_triggerno",&event.RAD_triggerno);   //
    dataTree->SetBranchAddress("PIX_xIntercept",&event.PIX_xIntercept); //
    dataTree->SetBranchAddress("PIX_yIntercept",&event.PIX_yIntercept); //
    dataTree->SetBranchAddress("PIX_xSlope",&event.PIX_xSlope);         //
    dataTree->SetBranchAddress("PIX_ySlope",&event.PIX_ySlope);         //
    dataTree->SetBranchAddress("PIX_chi2",&event.PIX_chi2);             //
    
    /* Get the number of Entries */
    int nev = dataTree->GetEntries();
    
    /* Book the Histograms */
    TH1F *h_rmsHists[2][2][2][4];        //[low/high][nocuts/withcuts][ch1234]
    
    TProfile *tp_Profiles[2][4];
    TProfile *rms_Profiles[2][2][2][4]; //[low/high][nocuts/withcuts][down/up][channel]
    
    TH1F *pedHists[2][2][8];

    char c_name[50];
    for (int i = 0; i!=4; ++i)
    {
        //Histograms of the root mean squares (RMS) of the pedestal values
        //RMS = sqrt(sum(sig-ped)^2/numsamples)
        sprintf(c_name,"rms_hist_LowGain_noCuts_Downstream_%d",i+1);        h_rmsHists[0][0][0][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_LowGain_noCuts_Upstream_%d",i+1);          h_rmsHists[0][0][1][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_LowGain_cuts_Downstream_%d",i+1);          h_rmsHists[0][1][0][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_LowGain_cuts_Upstream_%d",i+1);            h_rmsHists[0][1][1][i]  = new TH1F(c_name,c_name,100,0,5);
        
        sprintf(c_name,"rms_hist_highGain_noCuts_Downstream_%d",i+1);       h_rmsHists[1][0][0][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_highGain_noCuts_Upstream_%d",i+1);         h_rmsHists[1][0][1][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_highGain_cuts_Downstream_%d",i+1);         h_rmsHists[1][1][0][i]  = new TH1F(c_name,c_name,100,0,5);
        sprintf(c_name,"rms_hist_highGain_cuts_Upstream_%d",i+1);           h_rmsHists[1][1][1][i]  = new TH1F(c_name,c_name,100,0,5);
        
        //TProfiles aka average waveforms for each channel
        sprintf(c_name,"tprofile_Low_Gain_Downstream_%d",i+1);              tp_Profiles[0][i]       = new TProfile(c_name,c_name,1024,0,204.8);
        sprintf(c_name,"tprofile_Low_Gain_Upstream_%d",i+1);                tp_Profiles[1][i]       = new TProfile(c_name,c_name,1024,0,204.8);
        
        /*TProfiles for waveforms RMS calculation*/
        sprintf(c_name,"RMS_tprofile_LowCut_LowGain_Downstream_%d",i+1);    rms_Profiles[0][0][0][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_LowCut_LowGain_Upstream_%d",i+1);      rms_Profiles[0][0][1][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_HighCut_LowGain_Downstream_%d",i+1);   rms_Profiles[1][0][0][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_HighCut_LowGain_Upstream_%d",i+1);     rms_Profiles[1][0][1][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_LowCut_HighGain_Downstream_%d",i+1);   rms_Profiles[0][1][0][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_LowCut_HighGain_Upstream_%d",i+1);     rms_Profiles[0][1][1][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_HighCut_HighGain_Downstream_%d",i+1);  rms_Profiles[1][1][0][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        sprintf(c_name,"RMS_tprofile_HighCut_HighGain_Upstream_%d",i+1);    rms_Profiles[1][1][1][i] = new TProfile(c_name,c_name,1024,0,204.8);    //[nocuts/withcuts][low/high][down/up][channel]
        
//
//        sprintf(c_name,"tprofile_High_Gain_Downstream_%d",i+1);
//        tp_Profiles[1][i] = new TProfile(c_name,c_name,1024,0,204.8);
//
        //Hists to look at typical empty signals.
        sprintf(c_name,"Peds_Low_Gain_Downstream_%d",i+1);                  pedHists[0][0][i] = new TH1F(c_name,c_name,500,-50,50);
        sprintf(c_name,"Peds_Low_Gain_Upstream_%d",i+1);                    pedHists[0][1][i] = new TH1F(c_name,c_name,500,-50,50);
        sprintf(c_name,"Peds_High_Gain_Downstream_%d",i+1);                 pedHists[1][0][i] = new TH1F(c_name,c_name,500,-50,50);
        sprintf(c_name,"Peds_High_Gain_Upstream_%d",i+1);                   pedHists[1][1][i] = new TH1F(c_name,c_name,500,-50,50);
    }
    
    TH2D *corrPlots[2][4];
    
    corrPlots[0][0] = new TH2D("corr00","NW Correlation - Low Gain  - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[0][1] = new TH2D("corr01","NE Correlation - Low Gain  - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[0][2] = new TH2D("corr02","SW Correlation - Low Gain  - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[0][3] = new TH2D("corr03","SE Correlation - Low Gain  - Up vs Down", 100, 0, 1000, 100, 0, 1000);
          
    corrPlots[1][0] = new TH2D("corr10","NW Correlation - High Gain - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[1][1] = new TH2D("corr11","NE Correlation - High Gain - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[1][2] = new TH2D("corr12","SW Correlation - High Gain - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    corrPlots[1][3] = new TH2D("corr13","SE Correlation - High Gain - Up vs Down", 100, 0, 1000, 100, 0, 1000);
    
    TH1F *sampleCounter[4];
    sampleCounter[0] = new TH1F("sampleCounter1", "number of Samples 1", 200,0,200);
    sampleCounter[1] = new TH1F("sampleCounter2", "number of Samples 2", 200,0,200);
    sampleCounter[2] = new TH1F("sampleCounter3", "number of Samples 3", 200,0,200);
    sampleCounter[3] = new TH1F("sampleCounter4", "number of Samples 4", 200,0,200);
    
    double amps[2][2][4];
    
    TCanvas *cutCanv[4][4];
    cutCanv[0][0] = new TCanvas("cutcanv00", "cutcanvLowDown0", 1024, 1024);
    cutCanv[0][1] = new TCanvas("cutcanv01", "cutcanvLowDown1", 1024, 1024);
    cutCanv[0][2] = new TCanvas("cutcanv02", "cutcanvLowDown2", 1024, 1024);
    cutCanv[0][3] = new TCanvas("cutcanv03", "cutcanvLowDown3", 1024, 1024);
    
    cutCanv[1][0] = new TCanvas("cutcanv10", "cutcanvLowUp0", 1024, 1024);
    cutCanv[1][1] = new TCanvas("cutcanv11", "cutcanvLowUp1", 1024, 1024);
    cutCanv[1][2] = new TCanvas("cutcanv12", "cutcanvLowUp2", 1024, 1024);
    cutCanv[1][3] = new TCanvas("cutcanv13", "cutcanvLowUp3", 1024, 1024);
    
    cutCanv[2][0] = new TCanvas("cutcanv20", "cutcanvHighDown0", 1024, 1024);
    cutCanv[2][1] = new TCanvas("cutcanv21", "cutcanvHighDown1", 1024, 1024);
    cutCanv[2][2] = new TCanvas("cutcanv22", "cutcanvHighDown2", 1024, 1024);
    cutCanv[2][3] = new TCanvas("cutcanv23", "cutcanvHighDown3", 1024, 1024);
    
    cutCanv[3][0] = new TCanvas("cutcanv30", "cutcanvHighUp0", 1024, 1024);
    cutCanv[3][1] = new TCanvas("cutcanv31", "cutcanvHighUp1", 1024, 1024);
    cutCanv[3][2] = new TCanvas("cutcanv32", "cutcanvHighUp2", 1024, 1024);
    cutCanv[3][3] = new TCanvas("cutcanv33", "cutcanvHighUp3", 1024, 1024);
    
    for (int i = 0; i != 4; ++i)
    {
        for (int j = 0; j != 4; ++j)
        {
            cutCanv[i][j]->Divide(4,4);
        }
    }
    
    TGraph *g[100];
    
    for(int iev=0; iev < nev; ++iev)
    {
        if (iev%5000 == 1) cout << "Processing event: " << iev << endl;
        dataTree->GetEntry(iev);    // Get the iev'th event
        
        double d_rmsValue[2][2][4]; //[Low/High][Down/Up][NW/NE/SW/SE]
        double d_eventPed[2][2][4]; //[Low/High][down/up][NW/NE/SW/SE]
        double d_integral[2][2][4]; //[Low/High][down/up][NW/NE/SW/SE]
        
        int i_numPedSamples = 100;
        int histCounter = 16;
        int counter = 0;
        double d_pedSubbedAmps[1024] = {0.};
        
        for(int i = 0; i!= 2; ++i)              //  Low gain / High Gain
        {
            for (int j = 0; j != 2; ++j)        //  Downstream / Upstream
            {
                for (int k = 0; k!= 4; ++k)     //  Channel 1 - 4
                {
                    d_rmsValue[i][j][k] = 0.0;
                    d_eventPed[i][j][k] = GetPedestal(event.RAD_y[i][j][k], i_numPedSamples, d_rmsValue[i][j][k], h_rmsHists[i][0][j][k],   h_rmsHists[i][1][j][k], pedHists[i][j][k], d_integral[i][j][k], d_pedSubbedAmps);

                    for (int m = 0; m != 1024; ++m)
                    {
                        if (d_rmsValue[i][j][k] < 3) rms_Profiles[0][i][j][k]->Fill(event.RAD_x[0][j],(-1*(event.RAD_y[i][j][k][m]-d_eventPed[i][j][k])));     //[nocuts/withcuts][low/high][down/up][channel]
                        if (d_rmsValue[i][j][k] > 3) rms_Profiles[1][i][j][k]->Fill(event.RAD_x[0][j],(-1*(event.RAD_y[i][j][k][m]-d_eventPed[i][j][k])));     //[nocuts/withcuts][low/high][down/up][channel]
                        d_pedSubbedAmps[m] = -1*(event.RAD_y[i][j][k][m]-d_eventPed[i][j][k]);
                    }
                    
                    if (iev < histCounter)
                    {
                        cutCanv[counter][k]->cd(iev+1);
                        g[iev] = new TGraph(1024,event.RAD_x[0],d_pedSubbedAmps);
                        g[iev]->SetMaximum(600);
                        g[iev]->SetMinimum(-600);
                        g[iev]->GetXaxis()->SetTitle("Time (ns)");
                        g[iev]->GetYaxis()->SetTitle("Amplitude (mV)");
                        g[iev]->Draw("AC Same");
                    }
                    
                    amps[i][j][k] = GetMaxAmplitude(event.RAD_y[i][j][k]);
                }
                counter +=1;
            }
        }
        for (int i = 0; i != 2; ++i)
        {
            for (int j = 0; j != 4; ++j)
            {
                corrPlots[i][j]->Fill(amps[i][0][j],amps[i][1][j]);
            }
        }
    }
    
    /* Defint the Cuts */
    bool b_PbG  = 1;//GetMaxAmplitude(event.PbG_y) < 50 ? 1 : 0;
    bool b_DNW  = 1;//GetMaxAmplitude(event.RAD_y[0][0][0]) > 200 ? 1 : 0;
    bool b_DNE  = 1;//GetMaxAmplitude(event.RAD_y[0][0][1]) > 200 ? 1 : 0;
    bool b_DSW  = 1;//GetMaxAmplitude(event.RAD_y[0][0][2]) > 200 ? 1 : 0;
    bool b_DSE  = 1;//GetMaxAmplitude(event.RAD_y[0][0][3]) > 200 ? 1 : 0;
    
//        double t[1024], amp[2][8][1024];
//
//        for(int i = 0; i != 8; ++i)
//        {
//            pedHists[0][i]->Fill(d_eventPed[0][i]);
//            pedHists[1][i]->Fill(d_eventPed[1][i]);
//        }
//
//
//        for (int i = 0; i != 8; ++i)
//        {
//            for (int j = 0; j != 1024; ++j)
//            {
//                t[j] = j*0.2;
//                if(i == 0)
//                {
//                    amp[0][0][j] = -1*(event.BTL_amp[1][4][j]-d_eventPed[0][0]);
//                    amp[1][0][j] = -1*(event.BTL_amp[1][6][j]-d_eventPed[1][0]);
//                }
//                else
//                {
//                    amp[0][i][j] = -1*(event.RAD_y[0][i][j]-d_eventPed[0][i]);
//                }
//                amp[1][i][j] = -1*(event.RAD_y[1][i][j]-d_eventPed[1][i]);
//                tp_Profiles[0][i]->Fill(t[j],amp[0][i][j]);
//                tp_Profiles[1][i]->Fill(t[j],amp[1][i][j]);
//            }
//        }
////        TGraph *waveGraph = new TGraph(1024,t[i],amp[i]);
////        waveGraph->Draw("ACSame");
//    }
//
//    tp_Profiles[0][0]->SetTitle("Average Waveform Low Gain U1");
//    tp_Profiles[0][1]->SetTitle("Average Waveform Low Gain U2");
//    tp_Profiles[0][2]->SetTitle("Average Waveform Low Gain U3");
//    tp_Profiles[0][3]->SetTitle("Average Waveform Low Gain U4");
//    tp_Profiles[0][4]->SetTitle("Average Waveform Low Gain D1");
//    tp_Profiles[0][5]->SetTitle("Average Waveform Low Gain D2");
//    tp_Profiles[0][6]->SetTitle("Average Waveform Low Gain D3");
//    tp_Profiles[0][7]->SetTitle("Average Waveform Low Gain D4");
//
//    tp_Profiles[1][0]->SetTitle("Average Waveform High Gain U1");
//    tp_Profiles[1][1]->SetTitle("Average Waveform High Gain U2");
//    tp_Profiles[1][2]->SetTitle("Average Waveform High Gain U3");
//    tp_Profiles[1][3]->SetTitle("Average Waveform High Gain U4");
//    tp_Profiles[1][4]->SetTitle("Average Waveform High Gain D1");
//    tp_Profiles[1][5]->SetTitle("Average Waveform High Gain D2");
//    tp_Profiles[1][6]->SetTitle("Average Waveform High Gain D3");
//    tp_Profiles[1][7]->SetTitle("Average Waveform High Gain D4");
//
    TFile *outputFile = new TFile("radAnalyzed.root","RECREATE");
    outputFile->cd();
    outputFile->mkdir("RMS_Calcs/");
    outputFile->mkdir("Pedestals/");
    outputFile->mkdir("TProfiles/");
    outputFile->mkdir("CorrPlots/");
    outputFile->mkdir("Canvases/");
    outputFile->mkdir("Counters/");
//
//    TCanvas *canvas_LowGainProfile = new TCanvas("lowGainProfileCanvas","lowGainProfileCanvas",1600,800);
//    canvas_LowGainProfile->Divide(4,2);
//
//    TCanvas *canvas_HighGainProfile = new TCanvas("highGainProfileCanvas","highGainProfileCanvas",1600,800);
//    canvas_HighGainProfile->Divide(4,2);
//
    outputFile->cd("Canvases/");
    TCanvas *corrCanv = new TCanvas("corrCanv", "Corr Canvas", 800,800);
    corrCanv->Divide(2,2);
//    corrCanv->SetGrayscale();
    
    for (int i = 0; i!= 4; ++i)
    {
        for (int j=0; j!= 4; ++j)
        {
            cutCanv[i][j]->Draw();
            cutCanv[i][j]->Write();
        }

        corrCanv->cd(i+1);
        gPad->SetGrid();
        switch(i)
        {
            case 0:
                corrPlots[0][0]->GetYaxis()->SetTitle("Max Amplitude NW Upstream (mV)");
                corrPlots[0][0]->GetXaxis()->SetTitle("Max Amplitude NW Downstream (mV)");
                corrPlots[1][0]->GetYaxis()->SetTitle("Max Amplitude NW Upstream (mV)");
                corrPlots[1][0]->GetXaxis()->SetTitle("Max Amplitude NW Downstream (mV)");
                break;
            case 1:
                corrPlots[0][1]->GetYaxis()->SetTitle("Max Amplitude NE Downstream (mV)");
                corrPlots[0][1]->GetXaxis()->SetTitle("Max Amplitude NE Downstream (mV)");
                corrPlots[1][1]->GetYaxis()->SetTitle("Max Amplitude NE Downstream (mV)");
                corrPlots[1][1]->GetXaxis()->SetTitle("Max Amplitude NE Downstream (mV)");
                break;
            case 2:
                corrPlots[0][2]->GetYaxis()->SetTitle("Max Amplitude SW Downstream (mV)");
                corrPlots[0][2]->GetXaxis()->SetTitle("Max Amplitude SW Downstream (mV)");
                corrPlots[1][2]->GetYaxis()->SetTitle("Max Amplitude SW Downstream (mV)");
                corrPlots[1][2]->GetXaxis()->SetTitle("Max Amplitude SW Downstream (mV)");
                break;
            case 3:
                corrPlots[0][3]->GetYaxis()->SetTitle("Max Amplitude SE Downstream (mV)");
                corrPlots[0][3]->GetXaxis()->SetTitle("Max Amplitude SE Downstream (mV)");
                corrPlots[1][3]->GetYaxis()->SetTitle("Max Amplitude SE Downstream (mV)");
                corrPlots[1][3]->GetXaxis()->SetTitle("Max Amplitude SE Downstream (mV)");
                break;
        }
//        corrPlots[i]->GetYaxis()->SetMaxDigits(3);
//        corrPlots[i]->GetYaxis()->SetMaxDigits(3);
//        corrPlots[i]->GetZaxis()->SetMaxDigits(3);
        corrPlots[0][i]->Draw("COLZ");
    }
    
    corrCanv->Draw();
    corrCanv->Write();
    
//    TCanvas *rmsCanvas = new TCanvas("rmsCanv","rms Canvas",800,800);
//    rmsCanvas->Divide(4,2);
    
    for(int i = 0; i != 4; ++i)
    {
//        canvas_LowGainProfile->cd(i+1);
//        tp_Profiles[0][i]->Draw();
//        tp_Profiles[0][i]->GetYaxis()->SetRangeUser(-20,100);
//        tp_Profiles[0][i]->GetYaxis()->SetTitle("Average Amplitude (mV)");
//        tp_Profiles[0][i]->GetXaxis()->SetTitle("Time (ns)");
        
        outputFile->cd("RMS_Calcs/");
        for (int j = 0; j!= 2; ++j)
        {
            for (int k = 0; k!= 2; ++k)
            {
                for (int l = 0; l!= 2; ++l)
                {
                    h_rmsHists[j][k][l][i]->Write();
                    rms_Profiles[j][k][l][i]->Write();
                }
            }
        }
        
//        rms_Profiles[0][0][0][i]->Write();
//        rms_Profiles[0][0][1][i]->Write();
//        rms_Profiles[1][0][0][i]->Write();
//        rms_Profiles[1][0][1][i]->Write();
//        rms_Profiles[0][1][0][i]->Write();
//        rms_Profiles[0][1][1][i]->Write();
//        rms_Profiles[1][1][0][i]->Write();
//        rms_Profiles[1][1][1][i]->Write();
        
//        rmsCanvas->cd(i);
//        rms_Profiles[0][0][0][i]->Draw();
//
//        rmsCanvas->cd(i+4);
//        rms_Profiles[0][0][1][i]->
        
        outputFile->cd("CorrPlots/");
        corrPlots[0][i]->Write();
        
        outputFile->cd("Counters/");
        sampleCounter[i]->Write();
//        tp_Profiles[0][i]->Write();
        
        outputFile->cd("Pedestals/");
        pedHists[0][0][i]->Write();
        pedHists[0][1][i]->Write();
        pedHists[1][0][i]->Write();
        pedHists[1][1][i]->Write();
    }
    

    
//    int LineColors[8] = {40,30,46,8,9,12,32,49};
    gStyle->SetPalette(kRust);
    THStack *rmsCutStack = new THStack("rmsCutStack","RMS Hists Cuts");
    rmsCutStack->SetMaximum(6000);
    
    THStack *rmsNoCutStack = new THStack("rmsNoCutStack","RMS Hists No Cuts");
    rmsNoCutStack->SetMaximum(6000);
    
    for (int i = 0; i != 4; ++i)
    {
//        h_rmsUpHists[i]->SetLineColor(LineColors[i]);
        h_rmsHists[0][0][0][i]->SetLineWidth(4);
        h_rmsHists[0][1][0][i]->SetLineWidth(4);
//        h_rmsDownHists[i]->SetLineColor(LineColors[i+4]);
        h_rmsHists[0][0][1][i]->SetLineWidth(4);
        h_rmsHists[0][1][1][i]->SetLineWidth(4);
        rmsNoCutStack->Add(h_rmsHists[0][0][0][i]);
        rmsNoCutStack->Add(h_rmsHists[0][0][1][i]);
        rmsCutStack->Add(h_rmsHists[0][1][0][i]);
        rmsCutStack->Add(h_rmsHists[0][1][1][i]);
    }
    TCanvas *rmsCutCanv = new TCanvas("rmsCutCanv","RMS Cut Canvas", 800,1600);
    rmsCutCanv->Divide(1,2);
    gPad->SetGrid();
    gPad->SetFrameLineWidth(3);
    
    rmsCutCanv->cd(1);
    rmsNoCutStack->SetTitle("Pedestal RMS No Cuts for Low Gain Channels;Pedestal RMS (mV);# of Events");
    rmsNoCutStack->Draw("nostack PLC");
    
    rmsCutCanv->cd(2);
    rmsCutStack->SetTitle("Pedestal RMS With Cuts for Low Gain Channels;Pedestal RMS (mV);# of Events");
    rmsCutStack->Draw("nostack PLC");

    rmsCutCanv->Write();
    
    
//
//    for(int i = 0; i != 8; ++i)
//    {
//        canvas_HighGainProfile->cd(i+1);
//        tp_Profiles[1][i]->Draw();
//        tp_Profiles[1][i]->GetYaxis()->SetRangeUser(-60,120);
//        tp_Profiles[1][i]->GetYaxis()->SetTitle("Average Amplitude (mV)");
//        tp_Profiles[1][i]->GetXaxis()->SetTitle("Time (ns)");
//        tp_Profiles[1][i]->Write();
//
//        pedHists[0][i]->Write();
//        pedHists[1][i]->Write();
//    }
//
//    canvas_LowGainProfile->Write();
//    canvas_HighGainProfile->Write();

    outputFile->Close();
    
}
