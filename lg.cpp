/* To launch this program from within Mathematica use:
 *   In[1]:= link = Install["addtwo"]
 *
 * Or, launch this program from a shell and establish a
 * peer-to-peer connection.  When given the prompt Create Link:
 * type a port name. ( On Unix platforms, a port name is a
 * number less than 65536.  On Mac or Windows platforms,
 * it's an arbitrary word.)
 * Then, from within Mathematica use:
 *   In[1]:= link = Install["portname", LinkMode->Connect]
 */

#include <armadillo>
#include <sstream>
#include <array>
#include <cstdio>
#include <wstp.h>
#include "Laguerre.h"

#define ERR_MSG_LEN 512

using namespace arma;

extern int WSMain(int, char **);
//#ifdef WINDOWS_WSTP
//#include <Winuser.h>
//extern HWND WSInitializeIcon( HINSTANCE hinstCurrent, int nCmdShow);
//#endif

extern void laguerre_guassian_mag (
        int w, int h, double wscale, double hscale,
        int l, int p, double k, double w0, double z, 
        int clamp);

extern void laguerre_guassian_phase (
        int w, int h, double wscale, double hscale,
        int l, int p, double k, double w0, double z,
        int clamp);

extern void two_laguerre_guassian_phase (
        int w, int h, double wscale, double hscale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp);

extern void two_laguerre_guassian_mag (
        int w, int h, double wscale, double hscale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp);



void send_error(const char *type, const char *msg)
{
    char ss[ERR_MSG_LEN];
    snprintf(ss, ERR_MSG_LEN, "Message[%s, \"%s\"]", type, msg);

    WSClearError(stdlink);
    WSNewPacket(stdlink);
    WSEvaluate(stdlink, ss);
    WSNextPacket(stdlink);
    WSNewPacket(stdlink);
    WSPutSymbol(stdlink, "$Failed");
}

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

void laguerre_guassian_mag (
        int w, int h, double wscale, double hscale,
        int l, int p, double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, wscale, hscale)};

        dmat rp {*cartesian_to_r2p(xys)};
        cx_mat res_cx { *laguerre_guassian(rp, std::array<int,2>{l, p},
                std::array<double,3>{k, w0, z}) };
        dmat res_mag { abs(res_cx) };

        if (clamp)
        {
            clamp_inplace(res_mag);
        }

        const double *data = res_mag.memptr();
        int dims[]         = {w, h};
        int d              = 2;

        WSPutReal64Array(stdlink, data, dims, nullptr, d);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


void laguerre_guassian_phase (
        int w, int h, double wscale, double hscale,
        int l, int p, double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, wscale, hscale)};

        dmat rp {*cartesian_to_r2p(xys)};
        cx_mat res_cx { *laguerre_guassian(rp, std::array<int,2>{l, p},
                std::array<double,3>{k, w0, z}) };
        dmat res_mag { arg(res_cx) };

        const double *data = res_mag.memptr();
        int dims[]         = {w, h};
        int d              = 2;

        if (clamp)
            clamp_inplace(res_mag);

        WSPutReal64Array(stdlink, data, dims, nullptr, d);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void two_laguerre_guassian_phase (
        int w, int h, double wscale, double hscale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, wscale, hscale)};
        dmat rp {*cartesian_to_r2p(xys)};

        cx_mat res_cx { *two_laguerre_guassian(rp, 
                std::array<int, 2>{l0, p0},
                std::array<int, 2>{l1, p1},
                std::array<double, 3>{k, w0, z},
                std::array<double, 3>{k, w0, z}
                )};

        dmat res_phase { arg(res_cx) };
        if (clamp)
            clamp_inplace(res_phase);

        const double *data = res_phase.memptr();
        int dims[]         = {w, h};
        int d              = 2;

        WSPutReal64Array(stdlink, data, dims, nullptr, d);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


void two_laguerre_guassian_mag (
        int w, int h, double wscale, double hscale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, wscale, hscale)};
        dmat rp {*cartesian_to_r2p(xys)};

        cx_mat res_cx { *two_laguerre_guassian(rp, 
                std::array<int, 2>{l0, p0},
                std::array<int, 2>{l1, p1},
                std::array<double, 3>{k, w0, z},
                std::array<double, 3>{k, w0, z}
                )};

        dmat res_mag { abs(res_cx) };
        if (clamp)
            clamp_inplace(res_mag);

        const double *data = res_mag.memptr();
        int dims[]         = {w, h};
        int d              = 2;

        WSPutReal64Array(stdlink, data, dims, nullptr, d);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


//#if WINDOWS_WSTP
//
//int WINAPI WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, PSTR lpszCmdLine, int nCmdShow)
//{
//    char  buff[512];
//    char FAR * buff_start = buff;
//    char FAR * argv[32];
//    char FAR * FAR * argv_end = argv + 32;
//
//    hinstPrevious = hinstPrevious; /* suppress warning */
//
//    if( !WSInitializeIcon( hinstCurrent, nCmdShow)) return 1;
//    WSScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
//    return WSMain( (int)(argv_end - argv), argv);
//}
//
//#else

int main(int argc, char* argv[])
{
//#if WINDOWS_WSTP
//    ShowWindow(GetConsoleWindow(), SW_HIDE);
//#endif
    return WSMain(argc, argv);
}

//#endif
