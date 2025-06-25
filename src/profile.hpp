// profile.hpp
#pragma once

#include <vector>
#include <string>

#include "maths_ops.hpp"
#include "PhaseTransition.hpp"

/*
TO DO:
- move generate_stream_plot() elsewhere? independent of initial conditions so the same for all instances of FluidProfile
- add some safety thing for private vals and prof in FluidProfile to stop them from changing (helper static function to precompute in initialiser list)
- add dflt ctor for FluidProfile that just uses dflt ctor for PTParams
*/

namespace Hydrodynamics {

/**
 * @brief Computes the Lorentz factor between the wall frame and the universe frame.
 *
 * @param xi Fluid velocity in the wall frame.
 * @param v Fluid velocity in the universe frame.
 * 
 * @return Lorentz factor (mu).
 */
double mu(double xi, double v);

double dvdxi(double xi, double v, const double csq);
double dwdxi(double xi, double v, double w, const double csq);

double dxi_dtau(double xi, double v, const double csq);
double dv_dtau(double xi, double v, const double csq);
double dw_dtau(double xi, double v, double w, const double csq);

void generate_streamplot_data(const PhaseTransition::PTParams& params, int xi_pts, int y_pts, const std::string& filename);
void generate_streamplot_data(const PhaseTransition::PTParams& params);

// struct FluidState {
//     double v;   // velocity
//     double w;   // enthalpy
//     // double T;   // temperature

//     FluidState operator+(const FluidState& other) const {
//         // return {v + other.v, w + other.w, T + other.T};
//         return {v + other.v, w + other.w};
//     }
//     FluidState operator*(double scalar) const {
//         // return {v * scalar, w * scalar, T * scalar};
//         return {v * scalar, w * scalar};
//     }
// };

// struct FluidODE {
//   double dvdxi(double xi, double v, const double csq);
//   double dwdxi(double xi, double v, double w, const double csq);

//   double dxi_dtau(double xi, double v, const double csq);
//   double dv_dtau(double xi, double v, const double csq);
//   double dw_dtau(double xi, double v, double w, const double csq);
// };

using state_type = std::vector<double>;

// class FluidProfile_Veff {
//   public:
//     FluidProfile_Veff(const PhaseTransition::PTParams& params, const Thermo& therm);
  
//   private:
//     const state_type p_vals_, e_vals_, w_vals_; // p(T), e(T), w(T)
//     double cpsq_, cmsq_; // speed of sound in broken/symmetric phases
//     double pp_, pm_, ep_, em_, wp_, wm_; // pressure, energy density & enthalpy in b/s phases
//     int mode_; // hydrodynamic mode (deflagration=0, hybrid=1, detonation=2)
// };

/**
 * @class FluidProfile
 * @brief Represents the hydrodynamic profile of a bubble wall in a first-order phase transition.
 */
class FluidProfile {
  public:
    FluidProfile(const double vw, const double cpsq, const double cmsq, const double alN, const std::string method);
    FluidProfile(const PhaseTransition::PTParams& params); // bag EoS
    // FluidProfile(const Thermo& therm); // Veff EoS

    // are these needed?
    PhaseTransition::PTParams params() const { return params_; }; // PT parameters
    std::string get_eos_type() const { return method_; };

    // double xi0() const { return xi0_; };
    // double xif() const { return xif_; };
    // std::vector<double> init_state() const { return y0_; }; // Initial state vector {xi0, v0, w0}
    
    state_type xi_vals() const { return xi_vals_; }; // Vector of xi=r/t
    state_type v_vals() const { return v_vals_; }; // v(xi)
    state_type w_vals() const { return w_vals_; }; // w(xi)
    state_type la_vals() const { return la_vals_; }; // la(xi)

    void write(const std::string& filename = "bubble_prof.csv") const; // write bubble profile to disk
    void plot(const std::string& filename = "bubble_prof.png") const; // Plots bubble profiles

  private:
    const std::string method_;
    const double cpsq_, cmsq_, vw_, alN_;
    int mode_; // hydrodynamic mode (deflagration=0, hybrid=1, detonation=2)

    const PhaseTransition::PTParams params_; // local copy of PT parameters (get rid of this)
    
    state_type xi_vals_, v_vals_, w_vals_, la_vals_; // xi, v(xi), w(xi), la(x)

    int get_mode() const;
    double vJ_det(double alp);

    double calc_vm_bag(double vp, double alpha_p) const;
    double calc_vp_bag(double vm, double alpha_p) const;
    double calc_wm(double wp, double vp, double vm) const;
    double calc_w1wN(double xi_sh) const;

    double xi_shock(double v1UF) const; // position of shock front
    std::vector<double> get_alp_minmax(double vw, double cpsq, double cmsq) const;
    double get_alp_wall(double vpUF, double vw) const;
    double get_alp_shock(double vpUF, double v1UF, double alphaN) const;
    double v1UF_residual_func(double xif, double v1UF, const deriv_func& dydxi);

    // put number of integration points in input file? seems bad to hardcode
    std::vector<state_type> solve_profile(int n=100);
    state_type calc_lambda_vals(state_type w_vals) const;

    // testing purposes ONLY
    std::vector<state_type> read(const std::string& filename) const; // read bubble profile from disk

};

} // namespace Hydrodynamics