#include "Laguerre.h"
#include <armadillo>
#include <array>
#include <iomanip>

using namespace arma;
int main(int argc, char *argv[])
{
    dmat test = {
        {1, 20},
        {2, 3},
        {7, -2}
    };

    std::cout << "test l = 2, p = 3\n";
    std::cout << *laguerre_l(2, 3, test);

    //std::cout << "\ntest l = 1, p = 3\n";
    //std::cout << *laguerre_l(1, 3, test);

    //std::cout << "\ntest l = 3, p = 5\n";
    //std::cout << *laguerre_l(3, 5, test) << std::endl;

    std::cout << std::setprecision(10);
    std::array<int, 2> modes {3, 2};
    // k, w0, z
    std::array<double, 3> params {10, 2, 3};

    cx_mat res {*laguerre_guassian(test, modes, params)};

    res.raw_print("");
}
