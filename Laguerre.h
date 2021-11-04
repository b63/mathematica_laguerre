#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <memory>
#include <armadillo>



std::shared_ptr<arma::dmat> laguerre_l(int l, int p, const arma::dmat &x);
std::shared_ptr<arma::cx_mat> laguerre_guassian(const arma::dmat &r2p,
        const std::array<int, 2> &mode,
        const std::array<double, 3> &params);

std::shared_ptr<arma::dmat> get_points(int w, int h, double wscale, double hscale);
std::shared_ptr<arma::dmat> cartesian_to_r2p (const arma::dmat &xy);
std::shared_ptr<arma::cx_mat> two_laguerre_guassian(const arma::dmat &r2p,
        const std::array<int, 2> modes0,
        const std::array<int, 2> modes1,
        const std::array<double, 3> params0,
        const std::array<double, 3> params1
    );
#endif
