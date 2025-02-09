#include "oscillation_fitter.h"

// Define absolute paths to ROOT files
const std::string ROOT_FILES_DIR = "/home/jake/Projects/Fitter/StatOnly/sdu_xcheck/"; // Update this directory on one's own machine
const std::string EVENT_NORMAL_ROOT = ROOT_FILES_DIR + "xcheck_event_invert.root";
const std::string NUMU_FLUX_ROOT = ROOT_FILES_DIR + "numu_flux.root";
const std::string NUE_FLUX_ROOT = ROOT_FILES_DIR + "nue_flux.root";
const std::string NUMUBAR_FLUX_ROOT = ROOT_FILES_DIR + "numubar_flux.root";
const std::string NUEBAR_FLUX_ROOT = ROOT_FILES_DIR + "nuebar_flux.root";
const std::string XSEC_C12_ROOT = ROOT_FILES_DIR + "xcheck_xsec_tot_cc_C12.root";
const std::string XSEC_H1_ROOT = ROOT_FILES_DIR + "xcheck_xsec_tot_cc_H1.root";


// Global variable definitions
std::vector<TH2D*> expected_event_hists;
std::vector<TH2D*> observed_event_hists;
std::vector<TH2D*> flux_histograms;
std::vector<TH1D*> xsec_C12_histograms;
std::vector<TH1D*> xsec_H1_histograms;
std::vector<std::string> flavor_names;
std::vector<TH2D*> normal_event_hists;

std::vector<TH2D*> rebinned_expected_event_hists;
std::vector<TH2D*> rebinned_observed_event_hists;

// Function to load a TH2D histogram
TH2D* loadTH2D(const std::string& file, const std::string& hist) {
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file " << file << std::endl;
        exit(1);
    }
    TH2D* h = dynamic_cast<TH2D*>(f->Get(hist.c_str()));
    if (!h) {
        std::cerr << "Error: Could not find histogram " << hist << " in file " << file << std::endl;
        exit(1);
    }
    return h;
}

// Function to load a TH1D histogram
TH1D* loadTH1D(const std::string& file, const std::string& hist) {
    TFile* f = TFile::Open(file.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Could not open file " << file << std::endl;
        exit(1);
    }
    TH1D* h = dynamic_cast<TH1D*>(f->Get(hist.c_str()));
    if (!h) {
        std::cerr << "Error: Could not find histogram " << hist << " in file " << file << std::endl;
        exit(1);
    }
    return h;
}

// Convert a string to lowercase
// std::string to_lowercase(const std::string& str) {
//     std::string lower = str;
//     std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
//     return lower;
// }

// Initialize all necessary data
void initialize() {
    flavor_names = {"numu", "nue", "numubar", "nuebar"};

    TFile* event_normal_file = TFile::Open(EVENT_NORMAL_ROOT.c_str());
    if (!event_normal_file || event_normal_file->IsZombie()) {
        std::cerr << "Error: Could not open " << EVENT_NORMAL_ROOT << std::endl;
        exit(1);
    }

    for (const auto& flavor : flavor_names) {
        std::string hist_name = "final_" + flavor + "_flux";
        // std::string hist_name = "reduced_" + flavor + "_flux";
        TH2D* hist = dynamic_cast<TH2D*>(event_normal_file->Get(hist_name.c_str()));
        if (!hist) {
            std::cerr << "Error: Could not load histogram " << hist_name << std::endl;
            exit(1);
        }
        expected_event_hists.push_back(hist);

        // Initialize observed_event_hists with the same structure as expected_event_hists
        TH2D* observed_hist = dynamic_cast<TH2D*>(hist->Clone(("observed_" + flavor).c_str()));
        observed_hist->Reset();
        observed_event_hists.push_back(observed_hist);
    }

    // Load flux histograms
    flux_histograms = {
        loadTH2D(NUMU_FLUX_ROOT, "custom_flux_hist"),
        loadTH2D(NUE_FLUX_ROOT, "custom_flux_hist"),
        loadTH2D(NUMUBAR_FLUX_ROOT, "custom_flux_hist"),
        loadTH2D(NUEBAR_FLUX_ROOT, "custom_flux_hist")
    };

    // Load cross-section histograms
    xsec_C12_histograms = {
        loadTH1D(XSEC_C12_ROOT, "nu_mu_tot_cc"),
        loadTH1D(XSEC_C12_ROOT, "nu_e_tot_cc"),
        loadTH1D(XSEC_C12_ROOT, "nu_mu_bar_tot_cc"),
        loadTH1D(XSEC_C12_ROOT, "nu_e_bar_tot_cc")
    };

    xsec_H1_histograms = {
        loadTH1D(XSEC_H1_ROOT, "nu_mu_tot_cc"),
        loadTH1D(XSEC_H1_ROOT, "nu_e_tot_cc"),
        loadTH1D(XSEC_H1_ROOT, "nu_mu_bar_tot_cc"),
        loadTH1D(XSEC_H1_ROOT, "nu_e_bar_tot_cc")
    };

    // observed_event_rate = (TH1D*)expected_event_rate->Clone("observed_event_rate");
    // observed_event_rate->Reset();
}

// Calculate oscillation probabilities using Prob3++
void calculate_oscillation_probabilities(
    double delta_m2_21, double delta_m2_32, double sin2_theta_12,
    double sin2_theta_13, double sin2_theta_23, double delta_CP,
    std::vector<TH2D*>& osc_probabilities) {

    NeutrinoPropagator *myNu;
    BargerPropagator *bNu; 

    bNu = new BargerPropagator();
    bNu->UseMassEigenstates(false);
    bNu->SetOneMassScaleMode(false);
    bNu->SetDefaultOctant(23, 2);
    bNu->SetWarningSuppression(true);
    myNu = bNu;

    if (sin2_theta_23 < 0.5) {
        bNu->SetDefaultOctant(23, 1);
    }

    // myNu->SetPathLength(15.0); // Approximate Earth path length

    double delta_m2_31 = delta_m2_32 + delta_m2_21;
    // double theta_12 = TMath::ASin(TMath::Sqrt(sin2_theta_12));
    // double theta_13 = TMath::ASin(TMath::Sqrt(sin2_theta_13));
    // double theta_23 = TMath::ASin(TMath::Sqrt(sin2_theta_23));
    double theta_12 = sin2_theta_12;
    double theta_13 = sin2_theta_13;
    double theta_23 = sin2_theta_23;

    const int NBinsEnergy = 400;
    const int Dm12ZenithNbin = 400;
    double EnergyBins[NBinsEnergy + 1];
    double Dm12ZenithEdge[Dm12ZenithNbin + 1];
    for (int i = 0; i <= NBinsEnergy; ++i) {
       EnergyBins[i] = 0.1 * pow(10.0, double(i) * log10(200.0) / NBinsEnergy);
    }
    for (int i = 0; i <= Dm12ZenithNbin; ++i) {
       Dm12ZenithEdge[i] = -1.0 + double(i) * (2.0) / Dm12ZenithNbin;
    }

    for (size_t i = 0; i < flux_histograms.size(); ++i) {
        // Determine kNuBar (1 for neutrino, -1 for antineutrino)
        int kNuBar = (flavor_names[i] == "numu" || flavor_names[i] == "nue") ? 1 : -1;

        // dis for disappearance, app for appearance of that flavor, for example, (1,1) and (2,1) 
        TH2D* osc_prob_dis_hist = (TH2D*)flux_histograms[i]->Clone(("osc_prob_" + flavor_names[i]).c_str());
        TH2D* osc_prob_app_hist = (TH2D*)flux_histograms[i]->Clone(("osc_prob_" + flavor_names[i]).c_str());
        osc_prob_dis_hist->Reset();
        osc_prob_app_hist->Reset();

        for (int e_bin = 1; e_bin <= osc_prob_dis_hist->GetNbinsX(); ++e_bin) {
            for (int cos_bin = 1; cos_bin <= osc_prob_dis_hist->GetNbinsY(); ++cos_bin) {
                double E_low = EnergyBins[e_bin - 1];
                double E_high = EnergyBins[e_bin];
                double E = (E_low + E_high) / 2.0;
                double cos_low = Dm12ZenithEdge[cos_bin - 1];
                double cos_high = Dm12ZenithEdge[cos_bin];
                double cos = (cos_low + cos_high) / 2.0;

                // myNu->SetEnergy(E);
                myNu->SetMNS(theta_12, theta_13, theta_23, delta_m2_21, delta_m2_32, delta_CP, E, true, kNuBar);
                myNu->DefinePath(cos, 15.0);
                myNu->propagate(kNuBar);

                double prob_dis = (flavor_names[i] == "nue" || flavor_names[i] == "nuebar") ? 
                              myNu->GetProb(1, 1) : myNu->GetProb(2, 2);
                osc_prob_dis_hist->SetBinContent(e_bin, cos_bin, prob_dis);

                double prob_app = (flavor_names[i] == "nue" || flavor_names[i] == "nuebar") ? 
                              myNu->GetProb(2, 1) : myNu->GetProb(1, 2);
                osc_prob_app_hist->SetBinContent(e_bin, cos_bin, prob_app);
            }
        }
        osc_probabilities.push_back(osc_prob_dis_hist);
        osc_probabilities.push_back(osc_prob_app_hist);
    }
}

// // Chi-squared function for minimization
// void chi_squared_function(int& npar, double* grad, double& fval, double* par, int flag) {
//     fval = 0.0;

//     const double delta_m2_21   = par[0];
//     const double delta_m2_32   = -TMath::Abs(par[1]);
//     const double sin2_theta_12 = par[2];
//     const double sin2_theta_13 = par[3];
//     const double sin2_theta_23 = par[4];
//     const double delta_CP      = par[5];

//     // observed_event_rate->Reset();

//     // for (int kNuBar : {1, -1}) { // Loop over neutrino (1) and antineutrino (-1)
//         std::vector<TH2D*> osc_probabilities;
//         calculate_oscillation_probabilities(delta_m2_21, delta_m2_32, sin2_theta_12, sin2_theta_13, sin2_theta_23, delta_CP, osc_probabilities);

//         double E = 0.0;
//         double O = 0.0;

//         for (size_t i = 0; i < flux_histograms.size(); i += 2) { // Step through numu and nue together
//             TH2D* numu_flux = flux_histograms[i];
//             TH2D* nue_flux = flux_histograms[i + 1];
//             TH2D* numu_prob_dis = osc_probabilities[i*2];
//             TH2D* numu_prob_app = osc_probabilities[i*2 + 1];
//             TH2D* nue_prob_dis = osc_probabilities[i*2 + 2];
//             TH2D* nue_prob_app = osc_probabilities[i*2 + 3];
//             TH1D* numu_xsec_C12 = xsec_C12_histograms[i];
//             TH1D* numu_xsec_H1 = xsec_H1_histograms[i];
//             TH1D* nue_xsec_C12 = xsec_C12_histograms[i + 1];
//             TH1D* nue_xsec_H1 = xsec_H1_histograms[i + 1];

//             // TH1D* temp_observed_event_rate = (TH1D*)observed_event_rate->Clone("temp_observed_event_rate");
//             // temp_observed_event_rate->Reset();

//             for (int e_bin = 1; e_bin <= numu_flux->GetNbinsX(); ++e_bin) {
//                 for (int cos_bin = 1; cos_bin <= numu_flux->GetNbinsY(); ++cos_bin) {
//                     double numu_flux_val = numu_flux->GetBinContent(e_bin, cos_bin);
//                     double nue_flux_val = nue_flux->GetBinContent(e_bin, cos_bin);

//                     double numu_prob_dis_val = numu_prob_dis->GetBinContent(e_bin, cos_bin);
//                     double numu_prob_app_val = numu_prob_app->GetBinContent(e_bin, cos_bin);
//                     double nue_prob_dis_val = nue_prob_dis->GetBinContent(e_bin, cos_bin);
//                     double nue_prob_app_val = nue_prob_app->GetBinContent(e_bin, cos_bin);

//                     double numu_event_C12 = (numu_flux_val * numu_prob_dis_val + nue_flux_val * numu_prob_app_val )
//                                             * numu_xsec_C12->GetBinContent(e_bin) * C12_factor;
//                     double numu_event_H1 = (numu_flux_val * numu_prob_dis_val + nue_flux_val * numu_prob_app_val ) 
//                                             * numu_xsec_H1->GetBinContent(e_bin) * H1_factor;

//                     double nue_event_C12 = (nue_flux_val * nue_prob_dis_val + numu_flux_val * nue_prob_app_val) 
//                                             * nue_xsec_C12->GetBinContent(e_bin) * C12_factor;
//                     double nue_event_H1 = (nue_flux_val * nue_prob_dis_val + numu_flux_val * nue_prob_app_val)
//                                              * nue_xsec_H1->GetBinContent(e_bin) * H1_factor;

//                     double total_event_numu = numu_event_C12 + numu_event_H1;
//                     double total_event_nue = nue_event_C12 + nue_event_H1;
//                     // temp_observed_event_rate->Fill(numu_flux->GetXaxis()->GetBinCenter(e_bin), total_event);

//                     observed_event_hists[i]->SetBinContent(e_bin, cos_bin, total_event_numu);
//                     observed_event_hists[i+1]->SetBinContent(e_bin, cos_bin, total_event_nue);
//                 }
//             }

//             // observed_event_rate->Add(temp_observed_event_rate);
//             // delete temp_observed_event_rate;
//         }
//         for (auto* hist : osc_probabilities) delete hist; // Clean up memory
//     // }

//     rebin_histograms();

//     for (size_t i = 0; i < flux_histograms.size(); i++) {
//         for (int e_bin = 1; e_bin <= rebinned_expected_event_hists[i]->GetNbinsX(); ++e_bin) {
//             for (int cos_bin = 1; cos_bin <= rebinned_expected_event_hists[i]->GetNbinsY(); ++cos_bin) {
//                 double O = rebinned_expected_event_hists[i]->GetBinContent(e_bin, cos_bin);
//                 double E = rebinned_observed_event_hists[i]->GetBinContent(e_bin, cos_bin);
//                 if (E > 0 && O > 0) fval += 2 * (E - O + O * TMath::Log(O / E));
//             }
//         }
//     }

//     // add pull terms to chi2
//     // double pull[6];
//     // pull[0] = pow((delta_m2_21 - 7.53e-5) / (0.18e-5), 2);
//     // pull[1] = pow((delta_m2_32 + 2.529e-3) / (0.029e-3), 2);
//     // pull[2] = pow((sin2_theta_12 - 0.307) / (0.013), 2);
//     // pull[3] = pow((sin2_theta_13 - 2.19e-2) / (0.07e-2), 2);
//     // pull[4] = pow((sin2_theta_23 - 0.553) / (0.02), 2);
//     // pull[5] = pow((delta_CP - 1.19*TMath::Pi()) / (0.22*TMath::Pi()), 2);

//     // pull[0] = pow((delta_m2_21 - 7.53e-5) / (0.18e-5), 2);
//     // pull[1] = pow((delta_m2_32 - 2.455e-3) / (0.028e-3), 2);
//     // pull[2] = pow((sin2_theta_12 - 0.307) / (0.013), 2);
//     // pull[3] = pow((sin2_theta_13 - 2.19e-2) / (0.07e-2), 2);
//     // pull[4] = pow((sin2_theta_23 - 0.558) / (0.018), 2);
//     // pull[5] = pow((delta_CP - 1.19*TMath::Pi()) / (0.22*TMath::Pi()), 2);

//     // for (int i = 0; i < 6; i++) {
//     //     fval += pull[i];
//     // }
// }




// Chi-squared function for minimization
void chi_squared_function() {
    double chi = 0.0;

    // IO
    // const double delta_m2_21   = 7.53e-5;
    // const double delta_m2_32   = -TMath::Abs(2.529e-3);
    // const double sin2_theta_12 = 0.307;
    // const double sin2_theta_13 = 2.19e-2;
    // const double sin2_theta_23 = 0.553;
    // const double delta_CP      = 1.19*TMath::Pi();

    // NO
    const double delta_m2_21   = 7.53e-5;
    const double delta_m2_32   = TMath::Abs(2.455e-3);
    const double sin2_theta_12 = 0.307;
    const double sin2_theta_13 = 2.19e-2;
    const double sin2_theta_23 = 0.558;
    const double delta_CP      = 1.19*TMath::Pi();

    // pull[0] = pow((delta_m2_21 - 7.53e-5) / (0.18e-5), 2);
    // pull[1] = pow((delta_m2_32 + 2.529e-3) / (0.029e-3), 2);
    // pull[2] = pow((sin2_theta_12 - 0.307) / (0.013), 2);
    // pull[3] = pow((sin2_theta_13 - 2.19e-2) / (0.07e-2), 2);
    // pull[4] = pow((sin2_theta_23 - 0.553) / (0.02), 2);
    // pull[5] = pow((delta_CP - 1.19*TMath::Pi()) / (0.22*TMath::Pi()), 2);

    // observed_event_rate->Reset();

    // for (int kNuBar : {1, -1}) { // Loop over neutrino (1) and antineutrino (-1)
        std::vector<TH2D*> osc_probabilities;
        calculate_oscillation_probabilities(delta_m2_21, delta_m2_32, sin2_theta_12, sin2_theta_13, sin2_theta_23, delta_CP, osc_probabilities);

        double E = 0.0;
        double O = 0.0;

        for (size_t i = 0; i < flux_histograms.size(); i += 2) { // Step through numu and nue together
            TH2D* numu_flux = flux_histograms[i];
            TH2D* nue_flux = flux_histograms[i + 1];
            TH2D* numu_prob_dis = osc_probabilities[i*2];
            TH2D* numu_prob_app = osc_probabilities[i*2 + 1];
            TH2D* nue_prob_dis = osc_probabilities[i*2 + 2];
            TH2D* nue_prob_app = osc_probabilities[i*2 + 3];
            TH1D* numu_xsec_C12 = xsec_C12_histograms[i];
            TH1D* numu_xsec_H1 = xsec_H1_histograms[i];
            TH1D* nue_xsec_C12 = xsec_C12_histograms[i + 1];
            TH1D* nue_xsec_H1 = xsec_H1_histograms[i + 1];

            // TH1D* temp_observed_event_rate = (TH1D*)observed_event_rate->Clone("temp_observed_event_rate");
            // temp_observed_event_rate->Reset();

            for (int e_bin = 1; e_bin <= numu_flux->GetNbinsX(); ++e_bin) {
                for (int cos_bin = 1; cos_bin <= numu_flux->GetNbinsY(); ++cos_bin) {
                    double numu_flux_val = numu_flux->GetBinContent(e_bin, cos_bin);
                    double nue_flux_val = nue_flux->GetBinContent(e_bin, cos_bin);

                    double numu_prob_dis_val = numu_prob_dis->GetBinContent(e_bin, cos_bin);
                    double numu_prob_app_val = numu_prob_app->GetBinContent(e_bin, cos_bin);
                    double nue_prob_dis_val = nue_prob_dis->GetBinContent(e_bin, cos_bin);
                    double nue_prob_app_val = nue_prob_app->GetBinContent(e_bin, cos_bin);

                    double numu_event_C12 = (numu_flux_val * numu_prob_dis_val + nue_flux_val * numu_prob_app_val )
                                            * numu_xsec_C12->GetBinContent(e_bin) * C12_factor;
                    double numu_event_H1 = (numu_flux_val * numu_prob_dis_val + nue_flux_val * numu_prob_app_val ) 
                                            * numu_xsec_H1->GetBinContent(e_bin) * H1_factor;

                    double nue_event_C12 = (nue_flux_val * nue_prob_dis_val + numu_flux_val * nue_prob_app_val) 
                                            * nue_xsec_C12->GetBinContent(e_bin) * C12_factor;
                    double nue_event_H1 = (nue_flux_val * nue_prob_dis_val + numu_flux_val * nue_prob_app_val)
                                             * nue_xsec_H1->GetBinContent(e_bin) * H1_factor;

                    double total_event_numu = numu_event_C12 + numu_event_H1;
                    double total_event_nue = nue_event_C12 + nue_event_H1;
                    // temp_observed_event_rate->Fill(numu_flux->GetXaxis()->GetBinCenter(e_bin), total_event);

                    observed_event_hists[i]->SetBinContent(e_bin, cos_bin, total_event_numu);
                    observed_event_hists[i+1]->SetBinContent(e_bin, cos_bin, total_event_nue);
                }
            }

            // observed_event_rate->Add(temp_observed_event_rate);
            // delete temp_observed_event_rate;
        }
        for (auto* hist : osc_probabilities) delete hist; // Clean up memory
    // }

    rebin_wing_histograms();

    for (size_t i = 0; i < flux_histograms.size(); i++) {
        for (int e_bin = 1; e_bin <= rebinned_expected_event_hists[i]->GetNbinsX(); ++e_bin) {
            for (int cos_bin = 1; cos_bin <= rebinned_expected_event_hists[i]->GetNbinsY(); ++cos_bin) {
                double O = rebinned_expected_event_hists[i]->GetBinContent(e_bin, cos_bin);
                double E = rebinned_observed_event_hists[i]->GetBinContent(e_bin, cos_bin); 
                if (E > 0 && O > 0) chi  += 2 * (E - O + O * TMath::Log(O / E));
            }
        }
    }

    std::cout << "Chi-squared value: " << chi << std::endl;

    // add pull terms to chi2
    // double pull[6];
    // pull[0] = pow((delta_m2_21 - 7.53e-5) / (0.18e-5), 2);
    // pull[1] = pow((delta_m2_32 + 2.529e-3) / (0.029e-3), 2);
    // pull[2] = pow((sin2_theta_12 - 0.307) / (0.013), 2);
    // pull[3] = pow((sin2_theta_13 - 2.19e-2) / (0.07e-2), 2);
    // pull[4] = pow((sin2_theta_23 - 0.553) / (0.02), 2);
    // pull[5] = pow((delta_CP - 1.19*TMath::Pi()) / (0.22*TMath::Pi()), 2);

    // pull[0] = pow((delta_m2_21 - 7.53e-5) / (0.18e-5), 2);
    // pull[1] = pow((delta_m2_32 - 2.455e-3) / (0.028e-3), 2);
    // pull[2] = pow((sin2_theta_12 - 0.307) / (0.013), 2);
    // pull[3] = pow((sin2_theta_13 - 2.19e-2) / (0.07e-2), 2);
    // pull[4] = pow((sin2_theta_23 - 0.558) / (0.018), 2);
    // pull[5] = pow((delta_CP - 1.19*TMath::Pi()) / (0.22*TMath::Pi()), 2);

    // for (int i = 0; i < 6; i++) {
    //     fval += pull[i];
    // }
}


// Main function
int main(int argc, char** argv) {
    try {
        // Initialize the fitter
        initialize();
        
        // Run the fit
        chi_squared_function();

        // save_rebinned_histograms("check_event_hists.root");

        // Output success
        std::cout << "Fitter executed successfully." << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unknown error occurred." << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}















// // Main fitter function
// void compute_final_flux_fitter() {
//     // Initialize TMinuit with one parameter
//     TMinuit minuit(6);
//     minuit.SetFCN(chi_squared_function);

//     // // Define the parameter for the fit
//     minuit.DefineParameter(0, "delta_m2_21", 7.65e-5, 1e-7, 5.17e-6, 2.01e-3);
//     minuit.DefineParameter(1, "delta_m2_32", 2.529e-3, 1e-6, 0.0, 2.371e-1);
//     minuit.DefineParameter(2, "sin2_theta_12", 0.307, 0.001, 0.0, 1.0);
//     minuit.DefineParameter(3, "sin2_theta_13", 0.0219, 0.001, 0.0, 1.0);
//     minuit.DefineParameter(4, "sin2_theta_23", 0.553, 0.001, 0.50, 1.0);
//     minuit.DefineParameter(5, "delta_CP", 1.19 * TMath::Pi(), 0.01, 0.75 * TMath::Pi(), 1.63 * TMath::Pi());

//         // Define the parameter for the fit
//     // minuit.DefineParameter(0, "delta_m2_21", 1.531e-4, 1e-7, 0.0, 1.0);
//     // minuit.DefineParameter(1, "delta_m2_32", 4.855e-3, 1e-6, 0.0, 1.0);
//     // minuit.DefineParameter(2, "sin2_theta_12", 0.607, 0.001, 0.0, 1.0);
//     // minuit.DefineParameter(3, "sin2_theta_13", 0.0419, 0.001, 0.0, 1.0);
//     // minuit.DefineParameter(4, "sin2_theta_23", 0.778, 0.001, 0.5, 1.0);
//     // minuit.DefineParameter(5, "delta_CP", 0.59 * TMath::Pi(), 0.01, 0.0, 2.0 * TMath::Pi());

//     // Perform the minimization
//     minuit.Migrad();

//     // Retrieve and output the final chi-squared value
//     double fmin, fedm, errdef;
//     int npari, nparx, istat;
//     minuit.mnstat(fmin, fedm, errdef, npari, nparx, istat);

//     // Output the final chi-squared value
//     std::cout << "Fit completed successfully." << std::endl;
//     std::cout << "Final chi-squared value: " << fmin << std::endl;
//     std::cout << "Estimated distance to minimum (EDM): " << fedm << std::endl;
//     std::cout << "Number of parameters: " << npari << std::endl;
//     std::cout << "Fit status: " << istat << std::endl;
// }


// // Main function
// int main(int argc, char** argv) {
//     try {
//         // Initialize the fitter
//         initialize();
        
//         // Run the fit
//         compute_final_flux_fitter();

//         // save_rebinned_histograms("check_event_hists.root");

//         // Output success
//         std::cout << "Fitter executed successfully." << std::endl;
//     } catch (const std::exception& ex) {
//         std::cerr << "Error: " << ex.what() << std::endl;
//         return EXIT_FAILURE;
//     } catch (...) {
//         std::cerr << "Unknown error occurred." << std::endl;
//         return EXIT_FAILURE;
//     }
//     return EXIT_SUCCESS;
// }

// Function to rebin histograms for expected and observed event rates
void rebin_histograms() {
    const int blockSizeX = 40; // Number of bins to combine along X (energy)
    const int blockSizeY = 40; // Number of bins to combine along Y (cosine)

    // Clear previous rebinned histograms to avoid duplication
    for (auto* hist : rebinned_expected_event_hists) {
        delete hist; // Free memory
    }
    rebinned_expected_event_hists.clear();

    for (auto* hist : rebinned_observed_event_hists) {
        delete hist; // Free memory
    }
    rebinned_observed_event_hists.clear();

    auto createReducedHistogram = [](const std::string& title, const std::string& name) {
        return new TH2D(name.c_str(), title.c_str(), 10, 0, 10, 10, 0, 10); // Use bin numbers as range
    };

    auto fillReducedHistogram = [&](TH2D* original, TH2D* reduced) {
        for (int bx = 0; bx < 10; ++bx) {
            for (int by = 0; by < 10; ++by) {
                double sum = 0.0;
                for (int ix = 1; ix <= blockSizeX; ++ix) {
                    for (int iy = 1; iy <= blockSizeY; ++iy) {
                        int globalX = bx * blockSizeX + ix;
                        int globalY = by * blockSizeY + iy;
                        sum += original->GetBinContent(globalX, globalY);
                    }
                }
                reduced->SetBinContent(bx + 1, by + 1, sum);
            }
        }
    };

    // Rebin expected_event_hists
    for (size_t i = 0; i < expected_event_hists.size(); ++i) {
        TH2D* original_hist = expected_event_hists[i];
        if (!original_hist) {
            std::cerr << "Null histogram in expected_event_hists at index " << i << std::endl;
            continue;
        }

        std::string new_hist_name = std::string(original_hist->GetName()) + "_rebinned";
        TH2D* reduced_hist = createReducedHistogram(original_hist->GetTitle(), new_hist_name);

        fillReducedHistogram(original_hist, reduced_hist);

        rebinned_expected_event_hists.push_back(reduced_hist); // Save rebinned histogram
    }

    // Rebin observed_event_hists
    for (size_t i = 0; i < observed_event_hists.size(); ++i) {
        TH2D* original_hist = observed_event_hists[i];
        if (!original_hist) {
            std::cerr << "Null histogram in observed_event_hists at index " << i << std::endl;
            continue;
        }

        std::string new_hist_name = std::string(original_hist->GetName()) + "_rebinned";
        TH2D* reduced_hist = createReducedHistogram(original_hist->GetTitle(), new_hist_name);

        fillReducedHistogram(original_hist, reduced_hist);

        rebinned_observed_event_hists.push_back(reduced_hist); // Save rebinned histogram
    }
}


// Function to save rebinned expected histograms to a ROOT file
void save_rebinned_expected_histograms(const std::string& output_file_name) {
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (!output_file.IsOpen()) {
        std::cerr << "Error: Could not open file " << output_file_name << " for writing." << std::endl;
        return;
    }

    std::set<std::string> unique_names;
    for (TH2D* hist : rebinned_expected_event_hists) {
        if (hist) {
            std::string hist_name = hist->GetName();
            if (unique_names.insert(hist_name).second) { // Insert only if unique
                hist->Write();
                std::cout << "Saved histogram: " << hist_name << " to file." << std::endl;
            } else {
                std::cerr << "Warning: Duplicate histogram name " << hist_name << " encountered. Skipping." << std::endl;
            }
        }
    }

    output_file.Close();
    std::cout << "All histograms saved to " << output_file_name << "." << std::endl;
}

// Function to save rebinned expected and observed histograms to a ROOT file
void save_rebinned_histograms(const std::string& output_file_name) {
    // Open a ROOT file for writing
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (!output_file.IsOpen()) {
        std::cerr << "Error: Could not open file " << output_file_name << " for writing." << std::endl;
        return;
    }

    std::set<std::string> unique_names; // To track unique histogram names

    // Save rebinned expected histograms
    for (size_t i = 0; i < rebinned_expected_event_hists.size(); ++i) {
        TH2D* hist = rebinned_expected_event_hists[i];
        if (hist) {
            std::string hist_name = hist->GetName();
            if (unique_names.insert(hist_name).second) { // Check if the name is unique
                hist->Write();
                std::cout << "Saved histogram: " << hist_name << " to file." << std::endl;
            } else {
                std::cerr << "Warning: Duplicate histogram name " << hist_name
                          << " encountered. Skipping to avoid overwrite." << std::endl;
            }
        } else {
            std::cerr << "Warning: Null pointer encountered in rebinned_expected_event_hists at index " << i << "." << std::endl;
        }
    }

    // Save rebinned observed histograms
    for (size_t i = 0; i < rebinned_observed_event_hists.size(); ++i) {
        TH2D* hist = rebinned_observed_event_hists[i];
        if (hist) {
            std::string hist_name = hist->GetName();
            if (unique_names.insert(hist_name).second) { // Check if the name is unique
                hist->Write();
                std::cout << "Saved histogram: " << hist_name << " to file." << std::endl;
            } else {
                std::cerr << "Warning: Duplicate histogram name " << hist_name
                          << " encountered. Skipping to avoid overwrite." << std::endl;
            }
        } else {
            std::cerr << "Warning: Null pointer encountered in rebinned_observed_event_hists at index " << i << "." << std::endl;
        }
    }

    // Close the file
    output_file.Close();
    std::cout << "All histograms saved to " << output_file_name << "." << std::endl;
}

void rebin_wing_histograms() {
    // Define new bin edges for the first dimension
    const int newBinsX = 9;
    const int newBinsY = 10;
    const int blockSizeY = 40; // Number of bins to combine along Y (cosine)
    
    // Precomputed bin index mappings for X dimension
    const int new_bin_indices_x[10] = {0, 135, 156, 173, 196, 221, 248, 281, 322, 399}; 

    // Clear previous rebinned histograms to avoid duplication
    for (auto* hist : rebinned_expected_event_hists) {
        delete hist;
    }
    rebinned_expected_event_hists.clear();

    for (auto* hist : rebinned_observed_event_hists) {
        delete hist;
    }
    rebinned_observed_event_hists.clear();

    auto createReducedHistogram = [](const std::string& title, const std::string& name) {
        return new TH2D(name.c_str(), title.c_str(), newBinsX, 0, newBinsX, newBinsY, 0, newBinsY);
    };

    auto fillReducedHistogram = [&](TH2D* original, TH2D* reduced) {
        for (int bx = 0; bx < newBinsX; ++bx) {
            int startX = new_bin_indices_x[bx];
            int endX = new_bin_indices_x[bx + 1]; // Upper bound for this bin

            for (int by = 0; by < newBinsY; ++by) {
                double sum = 0.0;
                
                for (int ix = startX + 1; ix <= endX; ++ix) { // +1 because TH2D bins are 1-based
                    for (int iy = 1; iy <= blockSizeY; ++iy) { // Grouping 40 bins in Y
                        int globalY = by * blockSizeY + iy;
                        sum += original->GetBinContent(ix, globalY);
                    }
                }

                reduced->SetBinContent(bx + 1, by + 1, sum); // Set rebinned value
            }
        }
    };

    // Rebin expected_event_hists
    for (size_t i = 0; i < expected_event_hists.size(); ++i) {
        TH2D* original_hist = expected_event_hists[i];
        if (!original_hist) {
            std::cerr << "Null histogram in expected_event_hists at index " << i << std::endl;
            continue;
        }

        std::string new_hist_name = std::string(original_hist->GetName()) + "_rebinned";
        TH2D* reduced_hist = createReducedHistogram(original_hist->GetTitle(), new_hist_name);

        fillReducedHistogram(original_hist, reduced_hist);
        rebinned_expected_event_hists.push_back(reduced_hist);
    }

    // Rebin observed_event_hists
    for (size_t i = 0; i < observed_event_hists.size(); ++i) {
        TH2D* original_hist = observed_event_hists[i];
        if (!original_hist) {
            std::cerr << "Null histogram in observed_event_hists at index " << i << std::endl;
            continue;
        }

        std::string new_hist_name = std::string(original_hist->GetName()) + "_rebinned";
        TH2D* reduced_hist = createReducedHistogram(original_hist->GetTitle(), new_hist_name);

        fillReducedHistogram(original_hist, reduced_hist);
        rebinned_observed_event_hists.push_back(reduced_hist);
    }
}






