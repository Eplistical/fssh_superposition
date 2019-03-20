#include <cstdlib>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include "misc/fmtstring.hpp"
#include "misc/ioer.hpp"
#include "misc/crasher.hpp"
#include "misc/randomer.hpp"
#include "misc/vector.hpp"
#include "misc/timer.hpp"
#include "misc/fft.hpp"
#include "misc/matrixop.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;

double L = 16;
int M = 256;
double mass = 1000.0;

double A = 0.1;
double B = 3.0;

double xI = -3.0;
double sigmax = 1.0;
double kxI = 30.0;

double xwall_left = -99999999.9;
double xwall_right = 99999999.9;

double init_s = 1.0;
int Nstep = 2000;
double dt = 0.1;

int output_step = 100;

string output_mod = "init_s";

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("L", po::value<double>(&L), "grid para, the grid is [-L/2, L/2]")
        ("M", po::value<int>(&M), "grid number")
        ("mass", po::value<double>(&mass), "mass")
        ("A", po::value<double>(&A), "potential para")
        ("B", po::value<double>(&B), "potential para")
        ("init_x", po::value<double>(&xI), "init x")
        ("sigma_x", po::value<double>(&sigmax), "init sigma x")
        ("init_px", po::value<double>(&kxI), "init px")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("xwall_left", po::value<double>(&xwall_left), " left boundary x value to check end")
        ("xwall_right", po::value<double>(&xwall_right), " right boundary x value to check end")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("dt", po::value<double>(&dt), "time step")
        ("output_step", po::value<int>(&output_step), "output step")
        ("output_mod", po::value<string>(&output_mod), "output mode, can be init_s or init_px")
        ;
    po::variables_map vm; 
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return false;
    }
    return true;
}

vector< complex<double> > cal_H(const double x) {
    vector< complex<double> > H(4);
    double theta = 0.5 * M_PI * (erf(B * x) + 1.0);
    H[0+0*2] = -cos(theta);
    H[0+1*2] = sin(theta);
    H[1+0*2] = conj(H[0+1*2]);
    H[1+1*2] = cos(theta);
    return A * H;
}

vector< complex<double> > myfftshift(const vector< complex<double> >& Z) 
{
    vector< complex<double> > rst(M);
    for (int j = 0; j < M/2; ++j) {
        rst[j] = Z[j+M/2];
        rst[j+M/2] = Z[j];
    }
    return rst;
}

bool check_end( const vector< complex<double> >& psi0, const vector< complex<double> >& psi1, 
                const vector<double>& xarr)
{
    // if the WF is outside the wall, return true to stop the program
    double n_outside = 0.0;
    for (int j = 0; j < M; ++j) {
        double x = xarr[j];
        if (x < xwall_left or x > xwall_right) {
            n_outside += pow(abs(psi0[j]), 2) + pow(abs(psi1[j]), 2);
        }
    }
    return (n_outside > 0.01);
}

void exact() {
    // para
    double dx, dkx;
    vector<double> xarr = linspace(-L/2, L/2, M, dx);
    dkx = 2 * M_PI / M / dx;
    vector<double> kxarr = (arange(M) - M / 2) * dkx;
    xwall_left = max(xwall_left, -L/2 * 0.9);
    xwall_right = min(xwall_right, L/2 * 0.9);
    double c0 = sqrt(1 - init_s);
    double c1 = sqrt(init_s);
    // construct TU on k grid
    vector< complex<double> > TU(M);
    for (int j = 0; j < M; ++j) {
        double kx = kxarr[j];
        TU[j] = exp(-matrixop::IMAGIZ * dt * (kx*kx) / 2.0 / mass);
    }
    TU = myfftshift(TU);
    // construct VU on x grid
    vector< complex<double> > V00(M), V01(M), V10(M), V11(M);
    vector< complex<double> > evts00(M), evts01(M), evts10(M), evts11(M);
    vector< complex<double> > H00(M), H01(M), H10(M), H11(M);
    for (int j = 0; j < M; ++j) {
        double x = xarr[j];
        vector< complex<double> > H = cal_H(x);
        vector<double> eva;
        vector< complex<double> > evt, evamat;
        matrixop::eigh(H, eva, evt);

        evamat.assign(4, 0.0);
        evamat[0+0*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[0]);
        evamat[1+1*2] = exp(-matrixop::IMAGIZ * dt / 2.0 * eva[1]);
        auto tmp = matrixop::matmatmatC(evt, evamat, evt, 2, 2);
        V00[j] = tmp[0+0*2];
        V01[j] = tmp[0+1*2];
        V10[j] = tmp[1+0*2];
        V11[j] = tmp[1+1*2];
        H00[j] = H[0+0*2];
        H01[j] = H[0+1*2];
        H10[j] = H[1+0*2];
        H11[j] = H[1+1*2];

        // adiab evts
        if (j == 0) {
            const complex<double> phase0 = evt[0+0*2] / abs(evt[0+0*2]);
            evt[0+0*2] *= conj(phase0);
            evt[1+0*2] *= conj(phase0);
            const complex<double> phase1 = evt[1+1*2] / abs(evt[1+1*2]);
            evt[0+1*2] *= conj(phase1); 
            evt[1+1*2] *= conj(phase1);
        }
        else {
            const complex<double> phase0 = conj(evts00[j-1]) * evt[0+0*2] + conj(evts10[j-1]) * evt[1+0*2];
            evt[0+0*2] *= conj(phase0);
            evt[1+0*2] *= conj(phase0);
            const complex<double> phase1 = conj(evts01[j-1]) * evt[0+1*2] + conj(evts11[j-1]) * evt[1+1*2];
            evt[0+1*2] *= conj(phase1);
            evt[1+1*2] *= conj(phase1);
        }
        evts00[j] = evt[0+0*2];
        evts01[j] = evt[0+1*2];
        evts10[j] = evt[1+0*2];
        evts11[j] = evt[1+1*2];
    }
    // initialized WF
    /*
    vector< complex<double> > psi0(M, 0.0), psi1(M, 0.0);
    for (int j = 0; j < M; ++j) {
        double x = xarr[j];
        psi0[j] = c0 * exp(matrixop::IMAGIZ * (kxI * x)) * exp(-pow((x - xI) / sigmax, 2));
        psi1[j] = c1 * exp(matrixop::IMAGIZ * (kxI * x)) * exp(-pow((x - xI) / sigmax, 2));
    }
    double nm = norm(psi0 | psi1);
    psi0 /= nm;
    psi1 /= nm;
    */
    // initialized WF on adiab
    vector< complex<double> > psi0(M, 0.0), psi1(M, 0.0);
    vector< complex<double> > psiad0(M, 0.0), psiad1(M, 0.0);
    for (int j = 0; j < M; ++j) {
        double x = xarr[j];
        psiad0[j] = c0 * exp(matrixop::IMAGIZ * (kxI * x)) * exp(-pow((x - xI) / sigmax, 2));
        psiad1[j] = c1 * exp(matrixop::IMAGIZ * (kxI * x)) * exp(-pow((x - xI) / sigmax, 2));
    }
    double nm = norm(psiad0 | psiad1);
    psiad0 /= nm;
    psiad1 /= nm;
    psi0 = evts00 * psiad0 + evts01 * psiad1;
    psi1 = evts10 * psiad0 + evts11 * psiad1;

    // covinience vairables
    vector<int> dim{ M };
    // statistics
    double KE0 = 0.0, KE1 = 0.0, PE = 0.0;
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    // propagate WF
    for (int istep = 0; istep < Nstep; ++istep) {
        // exp(-iVdt/2)
        auto psi_k0 = V00 * psi0 + V01 * psi1;
        auto psi_k1 = V10 * psi0 + V11 * psi1;
        // exp(-iTdt)
        psi_k0 = misc::fftn(psi_k0, dim);
        psi_k1 = misc::fftn(psi_k1, dim);
        psi_k0 *= TU;
        psi_k1 *= TU;
        // exp(-iVdt/2)
        psi_k0 = misc::ifftn(psi_k0, dim);
        psi_k1 = misc::ifftn(psi_k1, dim);
        psi0 = V00 * psi_k0 + V01 * psi_k1;
        psi1 = V10 * psi_k0 + V11 * psi_k1;
        // analysis & output
        if (istep % output_step == 0) {
            // get psi_k
            psi_k0 = myfftshift(misc::fftn(psi0, dim));
            psi_k1 = myfftshift(misc::fftn(psi1, dim));
            double nm = norm(psi_k0 | psi_k1);
            psi_k0 /= nm;
            psi_k1 /= nm;
            KE0 = KE1 = PE = 0.0;
            n0trans = n0refl = n1trans = n1refl = 0.0;
            px0trans = px0refl = px1trans = px1refl = 0.0;
            for (int j = 0; j < M; ++j) {
                if (kxarr[j] >= 0.0) {
                    n0trans += pow(abs(psi_k0[j]), 2);
                    n1trans += pow(abs(psi_k1[j]), 2);
                    px0trans += pow(abs(psi_k0[j]), 2) * kxarr[j];
                    px1trans += pow(abs(psi_k1[j]), 2) * kxarr[j];
                }
                else {
                    n0refl += pow(abs(psi_k0[j]), 2);
                    n1refl += pow(abs(psi_k1[j]), 2);
                    px0refl += pow(abs(psi_k0[j]), 2) * kxarr[j];
                    px1refl += pow(abs(psi_k1[j]), 2) * kxarr[j];
                }
                KE0 += pow(abs(psi_k0[j]), 2) * (kxarr[j]*kxarr[j]) / 2.0 / mass;
                KE1 += pow(abs(psi_k1[j]), 2) * (kxarr[j]*kxarr[j]) / 2.0 / mass;
            }
            px0trans /= (n0trans + 1e-16);
            px1trans /= (n1trans + 1e-16);
            px0refl /= (n0refl + 1e-16);
            px1refl /= (n1refl + 1e-16);
            PE = real(sum(conj(psi0) * H00 * psi0 + conj(psi0) * H01 * psi1 + conj(psi1) * H10 * psi0 + conj(psi1) * H11 * psi1));
            // output
            if (istep == 0) {
                ioer::info("# EXACT 2D DIAB OUTPUT");
                ioer::info("# para: ", " L = ", L, " M = ", M, " mass = ", mass, " A = ", A, " B = ", B, 
                                       " xI = ", xI, " sigmax = ", sigmax, " kxI = ", kxI," init_s = ", init_s, " c0 = ", c0, " c1 = ", c1,
                                       " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step);
                ioer::info("# dx = ", dx, " dkx = ", dkx, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right);
                ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "px0refl", "px1trans", "px1refl", "Etot");
            }
            ioer::tabout('#', istep * dt, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl, KE0 + KE1 + PE);
            // check end
            if (check_end(psi0, psi1, xarr) == true) {
                ioer::info("# check_end returns true");
                break;
            }
        }
    }
    // final output
    if (output_mod == "init_s") {
        ioer::tabout(init_s, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl);
    }
    else if (output_mod == "init_px") {
        ioer::tabout(kxI, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl);
    }
    else {
    }
}

int main(int argc, char** argv) {
    if (argparse(argc, argv) == false) {
        return 0;
    }
    randomer::seed(0);
    timer::tic();
    exact();
    ioer::info("# ", timer::toc());
    return 0;
}
