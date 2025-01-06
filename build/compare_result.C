#include <TFile.h>
#include <TH2D.h>
#include <iostream>

void sum_histograms(const char* root_file_path = "check_event_hists.root") {
    // Open the ROOT file
    TFile* root_file = TFile::Open(root_file_path);
    if (!root_file || root_file->IsZombie()) {
        std::cerr << "Error: Cannot open ROOT file!" << std::endl;
        return;
    }

    // Load the histograms
    // TH2D* hist_final_numu = (TH2D*)root_file->Get("final_numu_flux_rebinned");
    // TH2D* hist_final_nue = (TH2D*)root_file->Get("final_nue_flux_rebinned");
    // TH2D* hist_final_numubar = (TH2D*)root_file->Get("final_numubar_flux_rebinned");
    // TH2D* hist_final_nuebar = (TH2D*)root_file->Get("final_nuebar_flux_rebinned");

    TH2D* hist_final_numu = (TH2D*)root_file->Get("Prediction_hist_numu");
    TH2D* hist_final_nue = (TH2D*)root_file->Get("Prediction_hist_nue");
    TH2D* hist_final_numubar = (TH2D*)root_file->Get("Prediction_hist_numubar");
    TH2D* hist_final_nuebar = (TH2D*)root_file->Get("Prediction_hist_nuebar");

    TH2D* hist_observed_numu = (TH2D*)root_file->Get("observed_numu_rebinned");
    TH2D* hist_observed_nue = (TH2D*)root_file->Get("observed_nue_rebinned");
    TH2D* hist_observed_numubar = (TH2D*)root_file->Get("observed_numubar_rebinned");
    TH2D* hist_observed_nuebar = (TH2D*)root_file->Get("observed_nuebar_rebinned");

    // Check if histograms are loaded correctly
    if (!hist_final_numu || !hist_final_nue || !hist_final_numubar || !hist_final_nuebar ||
        !hist_observed_numu || !hist_observed_nue || !hist_observed_numubar || !hist_observed_nuebar) {
        std::cerr << "Error: Could not load all histograms!" << std::endl;
        root_file->Close();
        return;
    }

    // Calculate and output the sum for each histogram
    double sum_final_numu = hist_final_numu->Integral();
    double sum_final_nue = hist_final_nue->Integral();
    double sum_final_numubar = hist_final_numubar->Integral();
    double sum_final_nuebar = hist_final_nuebar->Integral();

    double sum_observed_numu = hist_observed_numu->Integral();
    double sum_observed_nue = hist_observed_nue->Integral();
    double sum_observed_numubar = hist_observed_numubar->Integral();
    double sum_observed_nuebar = hist_observed_nuebar->Integral();

    std::cout << "Total sum for final_numu_flux_rebinned: " << sum_final_numu << std::endl;
    std::cout << "Total sum for final_nue_flux_rebinned: " << sum_final_nue << std::endl;
    std::cout << "Total sum for final_numubar_flux_rebinned: " << sum_final_numubar << std::endl;
    std::cout << "Total sum for final_nuebar_flux_rebinned: " << sum_final_nuebar << std::endl;

    std::cout << "Total sum for observed_numu_rebinned: " << sum_observed_numu << std::endl;
    std::cout << "Total sum for observed_nue_rebinned: " << sum_observed_nue << std::endl;
    std::cout << "Total sum for observed_numubar_rebinned: " << sum_observed_numubar << std::endl;
    std::cout << "Total sum for observed_nuebar_rebinned: " << sum_observed_nuebar << std::endl;

    // Close the ROOT file
    root_file->Close();
}
