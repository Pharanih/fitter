#ifndef OSCILLATION_FITTER_H
#define OSCILLATION_FITTER_H

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TMinuit.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <string>
#include <set>

// Prob3++ headers
#include "Prob3plusplus/BargerPropagator.h"
#include "Prob3plusplus/NeutrinoPropagator.h"

// Constants
const double solid_angle = 0.2 * TMath::Pi();   // Solid angle factor
const double mass_lab = 20e6;                   // 20 kt in kg
const double unit_conversion = 1e-42;           // Cross-section units from 10^-38 cm² to m²
const double runtime = 6 * 365.25 * 24 * 3600;  // Seconds in 6 years

// Atom counts
const double carbon_ratio = 0.88;
const double molar_mass_carbon = 12.0106;  // g/mol
const double mass_carbon = mass_lab * carbon_ratio * 1e3;  // Convert to grams
const double moles_carbon = mass_carbon / molar_mass_carbon;
const double atoms_C12 = moles_carbon * 6.02214076e23;  // Avogadro's number

const double hydrogen_ratio = 0.12
const double molar_mass_hydrogen = 1.0079750;  // g/mol
const double mass_hydrogen = mass_lab * hydrogen_ratio * 1e3;  // Convert to grams
const double moles_hydrogen = mass_hydrogen / molar_mass_hydrogen;
const double atoms_H1 = moles_hydrogen * 6.02214076e23;  // Avogadro's number

// Factors for Event Rates
const double C12_factor = atoms_C12 * runtime * solid_angle * unit_conversion;
const double H1_factor = atoms_H1 * runtime * solid_angle * unit_conversion;

// Global variables for fitting
extern TH1D* expected_event_rate;  // Expected event rates under normal ordering
extern TH1D* observed_event_rate;  // Observed event rates under inverted ordering (to be calculated)
extern std::vector<TH2D*> flux_histograms;   // Flux histograms for each flavor
extern std::vector<TH1D*> xsec_C12_histograms;  // Cross-section histograms for C12
extern std::vector<TH1D*> xsec_H1_histograms;   // Cross-section histograms for H1
extern std::vector<std::string> flavor_names;   // Flavor names for identification
extern std::vector<TH2D*> normal_event_hists;   // Event rate histograms from normal ordering (from event_normal.root)

// Function declarations
TH2D* loadTH2D(TFile* file, const std::string& name);
TH1D* loadTH1D(TFile* file, const std::string& name);
void initialize();
void calculate_oscillation_probabilities(
    double delta_m2_21, double delta_m2_32, double sin2_theta_12, double sin2_theta_13, 
    double sin2_theta_23, double delta_CP, std::vector<TH2D*>& osc_probabilities);
void chi_squared_function(int& npar, double* grad, double& fval, double* par, int flag);
void compute_final_flux_fitter();
void rebin_single_histogram(TH2D*& hist, const double* e_bins, int num_e_bins,
                            const double* cos_bins, int num_cos_bins);
void rebin_histograms();
void save_rebinned_histograms(const std::string& output_file_name);

#endif // OSCILLATION_FITTER_H

