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

void cleanRADData()
{
    
    struct {
        Float_t RAD_x[2][1024];
        Float_t RAD_y[2][9][1024];
        Float_t BTL_y[2][9][1024];
        Int_t   RAD_triggerno;

        Float_t PIX_xSlope;
        Float_t PIX_ySlope;
        Float_t PIX_xIntercept;
        Float_t PIX_yIntercept;
        Float_t PIX_chi2;
    } event;
    
    struct {
        Float_t RAD_x[2][1024];
        Float_t RAD_y[2][2][4][1024];
        Float_t PbG_y[1024];

        Int_t   RAD_triggerno;

        Float_t PIX_xSlope;
        Float_t PIX_ySlope;
        Float_t PIX_xIntercept;
        Float_t PIX_yIntercept;
        Float_t PIX_chi2;
    } newEvent;
    
    TChain *dataTree = new TChain("tree");
    dataTree->Add("../Data/RUN61733.root");
    dataTree->Add("../Data/RUN61734.root");
    dataTree->Add("../Data/RUN61735.root");
    dataTree->Add("../Data/RUN61736.root");
    dataTree->Add("../Data/RUN61737.root");
    dataTree->Add("../Data/RUN61738.root");
    dataTree->Add("../Data/RUN61739.root");
    
    dataTree->SetBranchAddress("RAD_x",&event.RAD_x); //
    dataTree->SetBranchAddress("RAD_y",&event.RAD_y); //
    dataTree->SetBranchAddress("BTL_y",&event.BTL_y); //
    dataTree->SetBranchAddress("RAD_triggerno",&event.RAD_triggerno); //
    dataTree->SetBranchAddress("PIX_xIntercept",&event.PIX_xIntercept); //
    dataTree->SetBranchAddress("PIX_yIntercept",&event.PIX_yIntercept); //
    dataTree->SetBranchAddress("PIX_xSlope",&event.PIX_xSlope); //
    dataTree->SetBranchAddress("PIX_ySlope",&event.PIX_ySlope); //
    dataTree->SetBranchAddress("PIX_chi2",&event.PIX_chi2); //
    
    TFile *outFile = new TFile("RAD_Data.root","RECREATE");
    
    TTree *radTree = new TTree("radTree","radTree");
    radTree->Branch("RAD_x",&newEvent.RAD_x,"RAD_x[2][1024]/F"); //
    radTree->Branch("RAD_y",&newEvent.RAD_y,"RAD_y[2][2][4][1024]/F"); //
    radTree->Branch("PbG",&newEvent.PbG_y,"PbG[1024]/F");
    radTree->Branch("RAD_triggerno",&newEvent.RAD_triggerno); //
    radTree->Branch("PIX_xIntercept",&newEvent.PIX_xIntercept); //
    radTree->Branch("PIX_yIntercept",&newEvent.PIX_yIntercept); //
    radTree->Branch("PIX_xSlope",&newEvent.PIX_xSlope); //
    radTree->Branch("PIX_ySlope",&newEvent.PIX_ySlope); //
    radTree->Branch("PIX_chi2",&newEvent.PIX_chi2); //

    int nev = dataTree->GetEntries();
    
    for(int iev=0; iev<nev; ++iev)
    {
        if (iev%1000 == 1) cout << iev << endl;
        dataTree->GetEntry(iev);  // Get the iev'th event
        
        for (int i = 0; i != 1024; ++i)
        {
            newEvent.RAD_x[0][i] = event.RAD_x[0][i];
            
            //Low gain channels
            //RAD_y[low/high][down/up][NW/NE/SW/SE]
            newEvent.RAD_y[0][0][0][i] = event.RAD_y[0][6][i]; // NW_D_1
            newEvent.RAD_y[0][0][1][i] = event.RAD_y[0][5][i]; // NE_D_2
            newEvent.RAD_y[0][0][2][i] = event.RAD_y[0][7][i]; // SW_D_3
            newEvent.RAD_y[0][0][3][i] = event.RAD_y[0][4][i]; // SE_D_4
            
            newEvent.RAD_y[0][1][0][i] = event.BTL_y[1][4][i]; // NW_U_1
            newEvent.RAD_y[0][1][1][i] = event.RAD_y[0][1][i]; // NE_U_2
            newEvent.RAD_y[0][1][2][i] = event.RAD_y[0][2][i]; // SW_U_3
            newEvent.RAD_y[0][1][3][i] = event.RAD_y[0][3][i]; // SE_U_4
            
            //Hight gain channels
            newEvent.RAD_y[1][0][0][i] = event.RAD_y[1][5][i]; // NW_D_1
            newEvent.RAD_y[1][0][1][i] = event.RAD_y[1][4][i]; // NE_D_2
            newEvent.RAD_y[1][0][2][i] = event.RAD_y[1][6][i]; // SW_D_3
            newEvent.RAD_y[1][0][3][i] = event.RAD_y[1][3][i]; // SE_D_4
            
            newEvent.RAD_y[1][1][0][i] = event.BTL_y[1][6][i]; // NW_U_1
            newEvent.RAD_y[1][1][1][i] = event.RAD_y[1][0][i]; // NE_U_2
            newEvent.RAD_y[1][1][2][i] = event.RAD_y[1][1][i]; // SW_U_3
            newEvent.RAD_y[1][1][3][i] = event.RAD_y[1][2][i]; // SE_U_4
            
            newEvent.PbG_y[i] = event.RAD_y[0][0][i]; // PbG
        }
        
        newEvent.RAD_triggerno = event.RAD_triggerno;
        
        newEvent.PIX_xSlope     = event.PIX_xSlope;
        newEvent.PIX_ySlope     = event.PIX_ySlope;
        newEvent.PIX_xIntercept = event.PIX_xIntercept;
        newEvent.PIX_yIntercept = event.PIX_yIntercept;
        newEvent.PIX_chi2       = event.PIX_chi2 ;
        
        radTree->Fill();
    }
    radTree->Write();
    outFile->Close();
}
