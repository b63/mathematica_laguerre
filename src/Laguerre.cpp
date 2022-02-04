
#include <armadillo>
#include <memory>
#include <array>
#include <iomanip>
#include <cstdio>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Laguerre.h"


using namespace arma;

std::shared_ptr<dmat> clamp(const dmat &arr)
{
    double mn = arr.min();
    double mx = arr.max();
    double range = mx - mn;

    return std::make_shared<dmat>( (arr-mn)/range);
}

void clamp_inplace(dmat &arr)
{
    double mn = arr.min();
    double mx = arr.max();
    double range = mx - mn;
    arr.for_each([=](double &x){ x = (x-mn)/range; });
}

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

void print_2dpnts_mat(const dmat &mat)
{
    const size_t h = mat.n_rows;
    const size_t w = mat.n_cols;
    for (auto j = 0; j < h; ++j)
    {
        for (auto i = 0;i < w; ++i)
        {
            std::cout <<  mat(i*h+ j,0) << "," << mat(i*h+ j,1) << "  ";

        }
        std::cout << std::endl;
    }
}

void blaze_inplace(dmat &pnts,
        double angle, double scale, double phase_offset,
        size_t x, size_t y, size_t w, size_t h, bool clamped)
{
    std::shared_ptr<dmat> xy {get_points_reg(w, h, 1, false)};
    dmat offsets {scale*(cos(angle)*xy->col(0) + sin(angle)*xy->col(1))};
    offsets = phase_offset+offsets;

    //dmat offsets {amplitude*cos(2*PI*phase_offset + A*PI*scale*(cos(angle)*xy->col(0) + sin(angle)*xy->col(1))) + amp_offset};
    //dmat offsets {amplitude*(phase_offset + ) + amp_offset};
    offsets.reshape(w, h);

    auto view {pnts(span(x, x+w-1), span(y,y+h-1))};

    //std::ofstream outfile;
    //outfile.open("/tmp/log", ios::out | ios::app);
    //outfile << std::setprecision(22);

    view += offsets;
    if (clamped) {
        view.for_each([](double &x) { 
                if (x < 0)      x = 0;
                else if (x > 1) x = 1.0l;
            });
    } else {
        //view(span(112,150), span(69,69)).raw_print(outfile, "before wrap");
        view -= floor(view);
        //view(span(112,150), span(69,69)).raw_print(outfile, "after wrap");
    }
    //outfile.close();
}

std::unique_ptr<cx_dmat> propagate(const cx_dmat &field, double delta, double z, double k, bool ft)
{
    constexpr const static double TAU = 2*M_PI;
    constexpr const static double TAU2 = TAU*TAU;
    cx_dmat fft {delta*delta*fft2(field)};
    const size_t N = field.n_rows;
    double delta_f = 1/(delta*N); // field should be square matrix
    double k2 = k*k;

    dmat kxky {delta_f * *get_shifted_points(N)};
    dcolvec kx2 {TAU2 * kxky.col(0) % kxky.col(0)};
    dcolvec ky2 {TAU2 * kxky.col(1) % kxky.col(1)};
    cx_dcolvec kz2 {std::complex<double>(1,0)*(k2 -kx2 - ky2)};
    cx_dcolvec kz {sqrt(kz2)};

    cx_dcolvec propagator {exp(std::complex<double>(0,z)*kz)};
    cx_dmat ft_field {reshape(propagator, N, N)};
    //abs(ft_field).raw_print(std::cout, "propagator mag");
    ft_field %= fft;

    if (ft)
    {
        auto pfft {std::make_unique<cx_dmat>(shift(shift(ft_field, N/2, 0), N/2, 1))};
        return pfft;
    }

    std::unique_ptr<cx_dmat> pfield {std::make_unique<cx_dmat>(1.0l/(delta*delta)* ifft2(ft_field))};

    // field.raw_print(std::cout, "field:");
    // pfield->raw_print(std::cout, "pfield:");

    return pfield;
}


std::shared_ptr<dmat> spherical(double scale, double offset, size_t w, size_t h)
{
    std::shared_ptr<dmat> xy {get_points_reg(w, h, 1, true)};

    const static std::function<void(double&)> wrap_fn = [](double &v) {
                if (v < 0.0l) {
                    v = 1 - v + floor(v);
                } else if (v > 1.0l) {
                    v = v - floor(v);
                }
    };

    std::shared_ptr<dmat> offsets {std::make_shared<dmat>(offset + scale * sqrt(pow(xy->col(0), 2) + pow(xy->col(1),2)))};
    offsets->for_each(wrap_fn);
    offsets->reshape(w,h);
    return offsets;
}

std::unique_ptr<dmat> get_shifted_points(size_t n)
{
    // N should be even
    double N = static_cast<double>(n);
    dcolvec cycle_i (N);
    cycle_i(span(0, N/2-1)) = regspace<dcolvec>(0, 1, N/2.0-1);
    cycle_i(span(N/2, N-1)) = regspace<dcolvec>(-N/2.0, -1);

    size_t pnts = N*N;
    std::unique_ptr<dmat> xys {std::make_unique<dmat>(pnts, 2)};
    for (size_t i = 0; i < N; ++i)
    {
        double v = cycle_i[i];
        (*xys)(span(i*N, i*N+N-1), 0) = cycle_i;
        (*xys)(span(i*N, i*N+N-1), 1).for_each([=](double &y){ y = v; });
    }

    return xys;
}

std::shared_ptr<arma::dmat> get_points_reg(int w, int h, double scale, bool center)
{
    // sample evenly spaced (by delta) points
    // sampled grid has width w and heigth h
    return get_points_irreg(w, h, scale, scale, center);
}

std::shared_ptr<arma::dmat> get_points_irreg(int w, int h, double deltaw, double deltah, bool center)
{
    const size_t N = w*h;
    double hw = deltaw*w/2.0;
    double hh = deltah*h/2.0;

    double wstart = -hw, wend = hw;
    double hstart = -hh, hend = hh;

    if (!center) {
        wstart = 0; wend = deltaw*w;
        hstart = 0; hend = deltah*h;
    }

    dcolvec xv {regspace<dcolvec>(wstart, deltaw, wend-deltaw/2)};
    dcolvec yv {regspace<dcolvec>(hstart, deltah, hend-deltah/2)};

    std::shared_ptr<dmat> xys {std::make_shared<dmat>(N, 2)};
    for (size_t i = 0; i < h; ++i)
    {
        double v = yv[i];
        (*xys)(span(i*w, i*w+w-1), 0) = xv;
        (*xys)(span(i*w, i*w+w-1), 1).for_each([=](double &y){ y = v; });
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

std::shared_ptr<cx_mat> laguerre_guassian_zx(const dmat &r2pz,
        const std::array<int, 2>    &mode,
        const std::array<double, 2> &params)
{
    int l = mode.at(0);
    int p = mode.at(1);

    double k  = params.at(0);
    double w0 = params.at(1);
    dmat z {r2pz.col(2)};
    // w(z) = w0 * sqrt(1 + 4*z^2/(k w0^2)^2)
    // R(z) = (z^2+(k w0^2/2)^2 )/z
    // psi(z) = arctan(2 z/(k w0^2))

    // compute w(z), R(z), and psi(z)
    double zR = k*w0*w0/2;
    dmat z2 {pow(z, 2)};
    double zR2 = zR*zR;

    dmat w {w0 * sqrt(1+ z2/zR2)};
    dmat w2inv {pow(w%w, -1)};
    dmat R {z+zR2/z};
    dmat psi{atan(z/zR)};

    dmat phi {-l*r2pz.col(1) - psi*(2*p+l+1) + k*z};
    cx_mat phi_p {exp(std::complex<double>(0,1)*conv_to<cx_mat>::from(phi))};

    dmat mag { w0/w % exp(-r2pz.col(0)%w2inv) % pow(sqrt(2*r2pz.col(0))/w,l)  };
    mag %= *laguerre_l(p, l, 2*r2pz.col(0)%w2inv);

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

    // std::complex<double> factor {exp(std::complex<double>(0, M_PI*(1-modes0[1]-modes1[1])))};
    std::complex<double> factor {exp(std::complex<double>(0, M_PI))};
    *lg0 += factor * (*lg1);

    return lg0;
}

std::shared_ptr<cx_mat> two_laguerre_guassian_zx(const dmat &r2pz,
        const std::array<int, 2> modes0,
        const std::array<int, 2> modes1,
        const std::array<double, 2> params0,
        const std::array<double, 2> params1
    )
{
    std::shared_ptr<cx_mat> lg0 {laguerre_guassian_zx(r2pz, modes0, params0)};
    std::shared_ptr<cx_mat> lg1 {laguerre_guassian_zx(r2pz, modes1, params1)};

    // std::complex<double> factor {exp(std::complex<double>(0, M_PI*(1-modes0[1]-modes1[1])))};
    std::complex<double> factor {exp(std::complex<double>(0, M_PI))};
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

std::shared_ptr<dmat> cartesian_to_r2pz (const dmat &zx, double y, bool ZX)
{
    // zx : z (col 0), x (col 1)
    std::shared_ptr<dmat> result {std::make_shared<dmat>(zx.n_rows, 3)};
    dmat ymat {zx.n_rows, 1, fill::value(y)};
    result->col(0) = sum(pow(zx.col(1), 2), 1) + y*y;
    result->col(2) = zx.col(0);

    if (ZX) {
        // using y = x for all points
        result->col(1) = atan2(zx.col(1), ymat);
    } else {
        // using x = x for all points
        result->col(1) = atan2(ymat, zx.col(1));
    }

    return result;
}




