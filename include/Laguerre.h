#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <memory>
#include <armadillo>



std::shared_ptr<arma::dmat> laguerre_l(int l, int p, const arma::dmat &x);
std::shared_ptr<arma::cx_mat> laguerre_guassian(const arma::dmat &r2p,
        const std::array<int, 2> &mode,
        const std::array<double, 3> &params);

std::shared_ptr<arma::dmat> get_points(int w, int h, double scale, bool center=true, bool regular=false);
std::shared_ptr<arma::dmat> cartesian_to_r2p (const arma::dmat &xy);
std::shared_ptr<arma::cx_mat> two_laguerre_guassian(const arma::dmat &r2p,
        const std::array<int, 2> modes0,
        const std::array<int, 2> modes1,
        const std::array<double, 3> params0,
        const std::array<double, 3> params1
    );
void blaze_inplace(arma::dmat &pnts, double angle, double scale, double phase_offset,
        double amplitude, double amp_phase,
        size_t x, size_t y, size_t w, size_t h, bool clamped = false);
std::shared_ptr<arma::dmat> clamp(const arma::dmat &arr);
std::shared_ptr<arma::dmat> spherical(double scale, double offset, size_t w, size_t h);

std::shared_ptr<arma::cx_mat> laguerre_guassian_zx(const arma::dmat &r2pz,
        const std::array<int, 2>    &mode,
        const std::array<double, 2> &params);
std::shared_ptr<arma::cx_mat> two_laguerre_guassian_zx(const arma::dmat &r2pz,
        const std::array<int, 2> modes0,
        const std::array<int, 2> modes1,
        const std::array<double, 2> params0,
        const std::array<double, 2> params1
    );
std::shared_ptr<arma::dmat> cartesian_to_r2pz (const arma::dmat &zx, double y, bool ZX);

void clamp_inplace(arma::dmat &arr);
#endif
