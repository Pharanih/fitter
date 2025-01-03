#include <TFile.h>
#include <TH2D.h>
#include <iostream>

void sum_histograms(const char* root_file_path = "rebinned_observed_event_hists.root") {
    // Open the ROOT file
    TFile* root_file = TFile::Open(root_file_path);
    if (!root_file || root_file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return;
    }

    // Load the 4 histograms
    // TH2D* hist_numu = (TH2D*)root_file->Get("reduced_numu_flux");
    // TH2D* hist_nue = (TH2D*)root_file->Get("reduced_nue_flux");
    // TH2D* hist_numubar = (TH2D*)root_file->Get("reduced_numubar_flux");
    // TH2D* hist_nuebar = (TH2D*)root_file->Get("reduced_nuebar_flux");

    TH2D* hist_numu = (TH2D*)root_file->Get("observed_numu_rebinned");
    TH2D* hist_nue = (TH2D*)root_file->Get("observed_nue_rebinned");
    TH2D* hist_numubar = (TH2D*)root_file->Get("observed_numubar_rebinned");
    TH2D* hist_nuebar = (TH2D*)root_file->Get("observed_nuebar_rebinned");

    if (!hist_numu || !hist_nue || !hist_numubar || !hist_nuebar) {
        std::cerr << "Error: Could not load all histograms!" << std::endl;
        root_file->Close();
        return;
    }

    // Calculate and output the sum for each histogram
    double sum_numu = hist_numu->Integral();
    double sum_nue = hist_nue->Integral();
    double sum_numubar = hist_numubar->Integral();
    double sum_nuebar = hist_nuebar->Integral();

    std::cout << "Total sum for reduced_numu_flux: " << sum_numu << std::endl;
    std::cout << "Total sum for reduced_nue_flux: " << sum_nue << std::endl;
    std::cout << "Total sum for reduced_numubar_flux: " << sum_numubar << std::endl;
    std::cout << "Total sum for reduced_nuebar_flux: " << sum_nuebar << std::endl;

    // Close the root file
    root_file->Close();
}
