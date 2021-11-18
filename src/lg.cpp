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
#include <cstring>
#include <math.h>
#include "Laguerre.h"

#define ERR_MSG_LEN 512

using namespace arma;

extern int WSMain(int, char **);
#ifdef WINDOWS_WSTP
#include <Windows.h>
extern HWND WSInitializeIcon( HINSTANCE hinstCurrent, int nCmdShow);
#endif

extern void laguerre_guassian_mag (
        int w, int h, double scale,
        int l, int p, double k, double w0, double z, 
        int clamp);

extern void laguerre_guassian_mag_zx (
        int w, int h, double scale,
        int l, int p, double k, double w0, double x,
        int usex, int clamp);

extern void laguerre_guassian_phase_zx (
        int w, int h, double scale,
        int l, int p, double k, double w0, double x,
        int usex, int clamp);

extern void laguerre_guassian_phase (
        int w, int h, double scale,
        int l, int p, double k, double w0, double z,
        int clamp);

extern void two_laguerre_guassian_phase (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp);

extern void two_laguerre_guassian_mag (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp);

extern void two_laguerre_guassian_mag_zx (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double x, int usex,
        int clamp);

extern void two_laguerre_guassian_phase_zx (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double x, int usex,
        int clamp);

extern void blaze(const double *arr, double angle);


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

void WS_spherical(double scale, double offset, int w, int h)
{

    try {
        std::shared_ptr<dmat> arr {spherical(scale, offset, w, h)};

        double *d_ptr = arr->memptr();
        long dimensions[] = {h, w};
        int ret = WSPutDoubleArray(stdlink, d_ptr, dimensions, nullptr, 2);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (const char *s) {
        send_error("LG::error", s);
    } catch (const std::exception &e) {
        send_error("LG::error", e.what());
    } catch (const std::string &s) {
        send_error("LG::error", s.c_str());
    }

}

void WS_blaze(double angle, double scale, double offset, double amp, 
        double amp_offset, int x, int y, int w, int h,
        int clamped)
{
    try {
        long *dimensions;
        char **heads;
        long depth;
        double *data;

        int ret = WSGetDoubleArray(stdlink, &data, &dimensions, &heads, &depth);
        if (!ret)
            throw WSErrorMessage(stdlink);

        if (depth != 2) {
            throw std::string("size of array must be 2, size ") 
                    + std::to_string(depth) + std::string(" given");
        } else if (strcmp(*heads, "List") || strcmp(*(heads+1), "List")) {
            throw std::string("head is not type {List,List}, instead ") + std::string(*heads)
                + std::string(", ") + std::string(*(heads+1));
        }

        const int W = dimensions[1];
        const int H = dimensions[0];

        dmat pnts{data, (size_t)W, (size_t)H};
        blaze_inplace(pnts, angle, scale, offset, amp, amp_offset, x, y, w, h, clamped);

        double *d_ptr = pnts.memptr();
        ret = WSPutDoubleArray(stdlink, d_ptr, dimensions, nullptr, 2);
        WSReleaseDoubleArray(stdlink, data, dimensions, heads, depth);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (const char *s) {
        send_error("LG::error", s);
    } catch (const std::exception &e) {
        send_error("LG::error", e.what());
    } catch (const std::string &s) {
        send_error("LG::error", s.c_str());
    }

}

void laguerre_guassian_mag (
        int w, int h, double scale,
        int l, int p, double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, scale)};

        dmat rp {*cartesian_to_r2p(xys)};
        cx_mat res_cx { *laguerre_guassian(rp, std::array<int,2>{l, p},
                std::array<double,3>{k, w0, z}) };
        dmat res_mag { abs(res_cx) };

        if (clamp)
        {
            clamp_inplace(res_mag);
        }

        const double *data = res_mag.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        WSPutReal64Array(stdlink, data, dims, nullptr, d);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void laguerre_guassian_mag_zx (
        int w, int h, double scale,
        int l, int p, double k, double w0, double x, 
        int usex, int clamp)
{
    try {
        dmat zx {*get_points(w,h, scale)};

        dmat r2pz { *cartesian_to_r2pz(zx, x, usex) };
        cx_mat res_cx { *laguerre_guassian_zx(r2pz, std::array<int,2>{l, p},
                std::array<double,2>{k, w0}) };
        dmat res_mag { abs(res_cx) };

        if (clamp) {
            clamp_inplace(res_mag);
        }

        const double *data = res_mag.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void laguerre_guassian_phase_zx (
        int w, int h, double scale,
        int l, int p, double k, double w0, double x, 
        int usex, int clamp)
{
    try {
        dmat zx {*get_points(w,h, scale)};

        dmat r2pz {*cartesian_to_r2pz(zx, x, usex)};
        cx_mat res_cx { *laguerre_guassian_zx(r2pz, std::array<int,2>{l, p},
                std::array<double,2>{k, w0}) };
        dmat res_phase { arg(res_cx) + M_PI };

        if (clamp)
            res_phase /= 2*M_PI;

        const double *data = res_phase.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void laguerre_guassian_phase (
        int w, int h, double scale,
        int l, int p, double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h, scale)};

        dmat rp {*cartesian_to_r2p(xys)};
        cx_mat res_cx { *laguerre_guassian(rp, std::array<int,2>{l, p},
                std::array<double,3>{k, w0, z}) };
        dmat res_phase { (arg(res_cx) + M_PI)};

        const double *data = res_phase.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        if (clamp)
            res_phase /= 2*M_PI;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void two_laguerre_guassian_phase (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h,scale)};
        dmat rp {*cartesian_to_r2p(xys)};

        cx_mat res_cx { *two_laguerre_guassian(rp, 
                std::array<int, 2>{l0, p0},
                std::array<int, 2>{l1, p1},
                std::array<double, 3>{k, w0, z},
                std::array<double, 3>{k, w0, z}
                )};

        dmat res_phase { arg(res_cx) + M_PI };
        if (clamp)
            res_phase /= 2*M_PI;

        const double *data = res_phase.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


void two_laguerre_guassian_phase_zx (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double x, int usex,
        int clamp)
{
    try {
        dmat zx {*get_points(w,h,scale)};
        dmat r2pz {*cartesian_to_r2pz(zx, x, usex)};

        cx_mat res_cx { *two_laguerre_guassian_zx(r2pz, 
                std::array<int, 2>{l0, p0},
                std::array<int, 2>{l1, p1},
                std::array<double, 2>{k, w0},
                std::array<double, 2>{k, w0}
                )};

        dmat res_phase { arg(res_cx) + M_PI };
        if (clamp)
            res_phase /= 2*M_PI;

        const double *data = res_phase.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


void two_laguerre_guassian_mag_zx (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double x, int usex,
        int clamp)
{
    try {
        dmat zx {*get_points(w,h,scale)};
        dmat r2pz {*cartesian_to_r2pz(zx, x, usex)};

        cx_mat res_cx { *two_laguerre_guassian_zx(r2pz, 
                std::array<int, 2>{l0, p0},
                std::array<int, 2>{l1, p1},
                std::array<double, 2>{k, w0},
                std::array<double, 2>{k, w0}
                )};

        dmat res_mag { abs(res_cx) };
        if (clamp)
            clamp_inplace(res_mag);

        const double *data = res_mag.memptr();
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}

void two_laguerre_guassian_mag (
        int w, int h, double scale,
        int l0, int p0, int l1, int p1,
        double k, double w0, double z,
        int clamp)
{
    try {
        dmat xys {*get_points(w,h,scale)};
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
        int dims[]         = {h, w};
        int d              = 2;

        int ret = WSPutReal64Array(stdlink, data, dims, nullptr, d);
        if (!ret)
            throw WSErrorMessage(stdlink);
    } catch (std::exception &e) {
        send_error("LG::error", e.what());
    }
}


#if WINDOWS_WSTP

int WINAPI WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, PSTR lpszCmdLine, int nCmdShow)
{
   char  buff[512];
   char FAR * buff_start = buff;
   char FAR * argv[32];
   char FAR * FAR * argv_end = argv + 32;

   hinstPrevious = hinstPrevious; /* suppress warning */

   if( !WSInitializeIcon( hinstCurrent, nCmdShow)) return 1;
   WSScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
   return WSMain( (int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
//#if WINDOWS_WSTP
//    ShowWindow(GetConsoleWindow(), SW_HIDE);
//#endif
    return WSMain(argc, argv);
}

#endif
