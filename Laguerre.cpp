
#include <armadillo>
#include <memory>
#include <array>
#include <iomanip>
#include <cstdio>

#include "Laguerre.h"


using namespace arma;

std::shared_ptr<dmat> laguerre_l(int n, int a, const dmat &x)
{
    if (n == 0)
    {
        return std::make_shared<dmat>(size(x), fill::ones);
    }

    std::vector<std::shared_ptr<dmat>> table;


    table.reserve(n+1);

    for (size_t i = 0; i <= n; i++)
    {
        if (i == 0)
        {
            table.push_back(std::make_shared<dmat>(size(x), fill::ones));
        }
        else if (i == 1)
        {
            table.push_back(std::make_shared<dmat>(x));
            table[i]->for_each([=](double &v){v = 1 + a-v;});
        }
        else
        {
            // recurrence relation from 
            // https://en.wikipedia.org/wiki/Laguerre_polynomials#Generalized_Laguerre_polynomials
            double c1 = 2*(double)i-1+a;
            double c2 = (double)i-1+a;
            table.push_back(std::make_shared<dmat>(
                        (c1-x) % *table[i-1] - c2 * *table[i-2])
                    );
            *table[i] /= (double)i;
        }
    }

    return table[n];
}


std::shared_ptr<arma::dmat> get_points(int w, int h, double wscale, double hscale)
{
    const size_t N = w*h;
    double hw = w/2.0, hh = h/2.0;

    dcolvec xv {linspace<dcolvec>(-hw, hw, w)};
    dcolvec yv {linspace<dcolvec>(-hh, hh, h)};

    std::shared_ptr<dmat> xys {std::make_shared<dmat>(N, 2)};
    for (size_t i = 0; i < h; ++i)
    {
        double v = yv[i];
        (*xys)(span(i*w, i*w+w-1), 0).for_each([=](double &x){ x = v; });
        (*xys)(span(i*w, i*w+w-1), 1) = xv;
    }

    return xys;
}

std::shared_ptr<cx_mat> laguerre_guassian(const dmat &r2p,
        const std::array<int, 2>    &mode,
        const std::array<double, 3> &params)
{
    int l = mode.at(0);
    int p = mode.at(1);

    double k  = params.at(0);
    double w0 = params.at(1);
    double z  = params.at(2);
    // w(z) = w0 * sqrt(1 + 4*z^2/(k w0^2)^2)
    // R(z) = (z^2+(k w0^2/2)^2 )/z
    // psi(z) = arctan(2 z/(k w0^2))

    // compute w(z), R(z), and psi(z)
    double zR = k*w0*w0/2;
    double z2 = z*z;
    double zR2 = zR*zR;

    double w  = w0 * sqrt(1+ z2/zR2);
    double w2  = w*w;
    double R  = z+zR2/z;
    double psi = atan(z/zR);

    dmat phi     {-l*r2p.col(1) - psi*(2*p+l+1) + k*z};
    cx_mat phi_p {exp(std::complex<double>(0,1)*conv_to<cx_mat>::from(phi))};

    dmat mag { w0/w * exp(-r2p.col(0)/w2) % pow(sqrt(2*r2p.col(0))/w,l)  };
    mag %= *laguerre_l(p, l, 2*r2p.col(0)/w2);

    return std::make_shared<cx_mat>(conv_to<cx_mat>::from(mag) % phi_p);
}

std::shared_ptr<cx_mat> two_laguerre_guassian(const dmat &r2p,
        const std::array<int, 2> modes0,
        const std::array<int, 2> modes1,
        const std::array<double, 3> params0,
        const std::array<double, 3> params1
    )
{
    std::shared_ptr<cx_mat> lg0 {laguerre_guassian(r2p, modes0, params0)};
    std::shared_ptr<cx_mat> lg1 {laguerre_guassian(r2p, modes1, params1)};

    std::complex<double> factor {exp(std::complex<double>(0, 1-modes0[1]-modes1[1]))};
    *lg0 += factor * (*lg1);

    return lg0;
}

std::shared_ptr<dmat> cartesian_to_r2p (const dmat &xy)
{
    std::shared_ptr<dmat> result {std::make_shared<dmat>(size(xy))};
    result->col(0) = sum(pow(xy, 2), 1);
    result->col(1) = atan2(xy.col(1), xy.col(0));

    return result;
}




