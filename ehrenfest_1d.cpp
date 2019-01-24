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
#include "misc/matrixop.hpp"
#include "boost/numeric/odeint.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/program_options.hpp"

enum {
    HOP_UP,
    HOP_DN,
    HOP_RJ,
    HOP_FR
};

// BE CAREFUL! PHASE IS IMPORTANT HERE! (DC AND BERRY FORCE)

using namespace std;
namespace po = boost::program_options;
using boost::numeric::odeint::runge_kutta4;
using boost::math::erf;
using state_t = vector< complex<double> >;
const complex<double> zI(0.0, 1.0);
double A = 0.10;
double B = 3.0;
const double mass = 1000.0;
double init_x = -3.0;
double sigma_x = 0.5; 
double init_px = 30.0;
double sigma_px = 1.0; 
double init_s = 1.0;
double xwall_left = -5.0;
double xwall_right = 5.0;
int Nstep = 1000000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 2000;
bool enable_nume = true;
string output_mod = "init_s";

vector< complex<double> > lastevt;
vector<double> eva(2);
vector<double> Fx(4);
vector< complex<double> > dcx(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("A", po::value<double>(&A), "potential para A")
        ("B", po::value<double>(&B), "potential para B")
        ("init_x", po::value<double>(&init_x), "init x")
        ("sigma_x", po::value<double>(&sigma_x), "init sigma x")
        ("init_px", po::value<double>(&init_px), "potential para init_px")
        ("sigma_px", po::value<double>(&sigma_px), "init sigma px")
        ("init_s", po::value<double>(&init_s), "init surface")
        ("xwall_left", po::value<double>(&xwall_left), "wall on x direction to check end")
        ("xwall_right", po::value<double>(&xwall_right), "wall on x direction to check end")
        ("Ntraj", po::value<int>(&Ntraj), "# traj")
        ("Nstep", po::value<int>(&Nstep), "# step")
        ("output_step", po::value<int>(&output_step), "# step for output")
        ("dt", po::value<double>(&dt), "single time step")
        ("output_mod", po::value<string>(&output_mod), "output mode, init_s or init_px")
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

double cal_theta(const double x) {
    return 0.5 * M_PI * (erf(B * x) + 1);
}

double cal_der_theta(const double x) {
    return sqrt(M_PI) * B * exp(-B * B * x * x);
}

vector< complex<double> > cal_H(const double x) {
    const double theta = cal_theta(x);
    vector< complex<double> > H {
        -cos(theta), sin(theta), sin(theta), cos(theta)
    };
    return A * H;
}

vector< complex<double> > cal_nablaHx(const double x) {
    const double theta = cal_theta(x);
    const double der_theta = cal_der_theta(x);
    vector< complex<double> > nablaHx {
        sin(theta), cos(theta), cos(theta), -sin(theta)
    };
    return A * der_theta * nablaHx;
}

void init_state(state_t& state) {
    state.resize(4, matrixop::ZEROZ);
    // state = x, vx, c0, c1
    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_px, sigma_px) / mass); 
    state[2].real(sqrt(1.0 - init_s));
    state[3].real(sqrt(init_s));
}

void cal_info_nume(const double x)
{
    // nume
    vector< complex<double> > evt;
    matrixop::hdiag(cal_H(x), eva, evt);
    // correct phase
    if (not lastevt.empty()) {
        auto tmp = matrixop::matCmat(lastevt, evt, 2);
        for (int j = 0; j < 2; ++j) {
            complex<double> eip = tmp[j+j*2] / abs(tmp[j+j*2]);
            for (int k = 0; k < 2; ++k) {
                evt[k+j*2] /= eip;
            }
        }
    }
    // F, dc
    dcx = matrixop::matCmatmat(evt, cal_nablaHx(x), evt, 2, 2);
    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            Fx[j+k*2] = -dcx[j+k*2].real();
            if (j == k) {
                dcx[j+k*2] = 0.0;
            }
            else {
                dcx[j+k*2] /= (eva[k] - eva[j]);
            }
        }
    }
    lastevt = move(evt);
}

void cal_info(const double x)
{
    cal_info_nume(x);
}

void sys(const state_t& state, state_t& state_dot, const double /* t */) {
    // state = x, v, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    complex<double> a00 = c[0] * conj(c[0]);
    complex<double> a01 = c[0] * conj(c[1]);
    complex<double> a10 = c[1] * conj(c[0]);
    complex<double> a11 = c[1] * conj(c[1]);
    // state_dot
    state_dot[0] = vx;
    state_dot[1] = (a00 * Fx[0+0*2] + a01 * Fx[1+0*2] + a10 * Fx[0+1*2] + a11 * Fx[1+1*2]) / mass;
    state_dot[2] = -zI * c[0] * eva[0] - c[1] * (vx * dcx[0+1*2]) - c[0] * (vx * dcx[0+0*2]);
    state_dot[3] = -zI * c[1] * eva[1] - c[0] * (vx * dcx[1+0*2]) - c[1] * (vx * dcx[1+1*2]);
}

bool check_end(const state_t& state) {
    double x = state[0].real();
    double vx = state[1].real();
    return ((x > xwall_right and vx > 0.0) or (x < xwall_left and vx < 0.0));
}

void ehrenfest() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        init_state(state[itraj]);
    }
    // main loop
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    // statistics
    double pxtrans = 0.0, pxrefl = 0.0;
    double KE = 0.0, PE = 0.0;
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info(state[itraj][0].real());
                // propagate
                rk4.do_step(sys, state[itraj], istep * dt, dt);
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
            // momentum
            pxtrans = pxrefl = 0.0;
            for_each(state.begin(), state.end(), 
                    [&pxtrans, &pxrefl] (const state_t& st) {
                        double vx = st[1].real();
                        if (vx >= 0.0) {
                            pxtrans += mass * vx;
                        }
                        else {
                            pxrefl += mass * vx; 
                        }
                    }
                    );
            // energy
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&KE, &PE] (const state_t& st) {
                        double x = st[0].real();
                        double vx = st[1].real();
                        vector< complex<double> > c { st[2], st[3] };
                        complex<double> a00 = c[0] * conj(c[0]);
                        complex<double> a11 = c[1] * conj(c[1]);

                        vector<double> eva;
                        matrixop::hdiag(cal_H(x), eva);
                        KE += 0.5 * mass * vx * vx;
                        PE += eva[0] * a00.real() + eva[1] * a11.real(); 
                    }
                    );
            // output
            if (istep == 0) {
                // para & header
                ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step, " output_mod = ", output_mod,
                            " mass = ", mass, " A = ", A, " B = ", B, 
                            " init_x = ", init_x, " init_px = ", init_px, 
                            " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                            " xwall_left = ", xwall_left, " xwall_right = ", xwall_right
                        );
                ioer::tabout('#', "t", "pxtrans", "pxrefl", "Etot");
            }
            ioer::tabout('#', istep * dt, pxtrans / Ntraj, pxrefl / Ntraj, (KE + PE) / Ntraj);

            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                break;
            }
        }

    }
    if (output_mod == "init_px") {
        ioer::tabout(init_px, pxtrans, pxrefl);
    }
    else if (output_mod == "init_s") {
        ioer::tabout(init_s,  pxtrans / Ntraj, pxrefl / Ntraj);
    }
    else {
        misc::crasher::confirm(false, "invalid output mode");
    }
}

void check_surf() {
    for (double x = xwall_left - 1; x < xwall_right + 1; x += 0.01) {
        cal_info(x);
        ioer::tabout(x, eva, Fx, real(dcx));
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        ioer::info("use --help for detailed info");
        return -1;
    }
    else if (string(argv[1]) == "check") {
        check_surf();
        return 0;
    }
    else {
        if (argparse(argc, argv) == false) {
            return 0;
        }
        randomer::seed(0);
        timer::tic();
        ehrenfest();
        ioer::info("# ", timer::toc());
    }
    return 0;
}
