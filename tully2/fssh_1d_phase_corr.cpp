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
double B = 0.28;
double C = 0.015;
double D = 0.06;
double E0 = 0.05;
const double mass = 2000.0;
double init_x = -10.0;
double sigma_x = 0.5; 
double init_px = 30.0;
double sigma_px = 1.0; 
double init_s = 0.0;
double xwall_left = -12.0;
double xwall_right = 12.0;
int Nstep = 1000000;
double dt = 0.1;
int output_step = 100;
int Ntraj = 2000;
bool enable_hop = true;
bool enable_deco = false;
bool enable_nodjj = false;
bool enable_nume = true;
string output_mod = "init_s";

vector< complex<double> > lastevt;
vector<double> eva(2);
vector<double> Fx(2);
vector< complex<double> > dcx(4);

inline bool argparse(int argc, char** argv) 
{
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
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
        ("enable_hop", po::value<bool>(&enable_hop), "enable hopping")
        ("enable_deco", po::value<bool>(&enable_deco), "enable decoherence")
        ("enable_nodjj", po::value<bool>(&enable_nodjj), "let djj = 0")
        ("enable_nume", po::value<bool>(&enable_nume), "use numerical way to calc dc")
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

vector< complex<double> > cal_H(const double x) {
    vector< complex<double> > H(4);
    H[0+0*2] = 0.0;
    H[1+1*2] = -A * exp(-B * x * x) + E0;
    H[0+1*2] = C * exp(-D * x * x);
    H[1+0*2] = H[0+1*2];
    return H;
}

vector< complex<double> > cal_nablaHx(const double x) {
    vector< complex<double> > nablaHx(4);
    nablaHx[0+0*2] = 0.0;
    nablaHx[1+1*2] = 2.0 * A * B * x * exp(-B * x * x);
    nablaHx[0+1*2] = -2.0 * C * D * x * exp(-D * x * x);
    nablaHx[1+0*2] = nablaHx[0+1*2];
    return nablaHx;
}

void init_state(state_t& state) {
    state.resize(5, matrixop::ZEROZ);
    // state = x, vx, c0, c1, s
    state[0].real(randomer::normal(init_x, sigma_x)); 
    state[1].real(randomer::normal(init_px, sigma_px) / mass); 
    state[2].real(sqrt(1.0 - init_s));
    state[3].real(sqrt(init_s));
    state[4].real((randomer::rand() < init_s) ? 1.0 : 0.0);
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
    Fx.assign(2, 0.0);
    for (int j = 0; j < 2; ++j) {
        for (int k = 0; k < 2; ++k) {
            if (j == k) {
                Fx[j] = -dcx[j+k*2].real();
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
    int s = static_cast<int>(state[4].real());
    // state_dot
    state_dot[0] = vx;
    state_dot[1] = (Fx[s]) / mass;
    if (s == 0) {
        // vx <-> v0x
        double v0x = vx;
        double v1x = sqrt(max(0.0, v0x * v0x + 2 * mass * eva[0] - 2 * mass * eva[1]));
        state_dot[2] = zI * mass * v0x * v0x * c[0] - c[1] * (v0x * dcx[0+1*2]) - c[0] * (v0x * dcx[0+0*2]);
        state_dot[3] = zI * mass * v0x * v1x * c[1] - c[0] * (v0x * dcx[1+0*2]) - c[1] * (v0x * dcx[1+1*2]);
    }
    else {
        // vx <-> v1x
        double v1x = vx;
        double v0x = sqrt(max(0.0, v1x * v1x + 2 * mass * eva[1] - 2 * mass * eva[0]));
        state_dot[2] = zI * mass * v1x * v0x * c[0] - c[1] * (v1x * dcx[0+1*2]) - c[0] * (v1x * dcx[0+0*2]);
        state_dot[3] = zI * mass * v1x * v1x * c[1] - c[0] * (v1x * dcx[1+0*2]) - c[1] * (v1x * dcx[1+1*2]);
    }
    state_dot[4] = matrixop::ZEROZ;
}

int hopper(state_t& state) {
    // state = x, vx, c0, c1, s
    double x = state[0].real();
    double vx = state[1].real();
    vector< complex<double> > c { state[2], state[3] };
    int s = static_cast<int>(state[4].real());
    // calc hop prob
    double g = -2 * dt * (c[s] * conj(c[1-s]) * (vx * dcx[1-s+s*2])).real() / (c[s] * conj(c[s])).real();
    double dE = eva[1-s] - eva[s];
    // hop
    if (randomer::rand() < g) {
        double tmp = vx * vx - 2 * dE / mass;
        if (tmp > 0.0) {
            double vnew = sqrt(tmp);
            state[1] = vnew * (vx < 0.0 ? -1 : 1); 
            state[4].real(1.0 - s); 
            return (s == 0) ? HOP_UP : HOP_DN;
        }
        else {
            return HOP_FR;
        }
    }   
    return HOP_RJ;
}

bool check_end(const state_t& state) {
    double x = state[0].real();
    double vx = state[1].real();
    return ((x > xwall_right and vx > 0.0) or (x < xwall_left and vx < 0.0));
}

void fssh() {
    // propagation variables
    runge_kutta4<state_t> rk4;
    vector<state_t> state(Ntraj);
    double hopup = 0.0, hopdn = 0.0, hopfr = 0.0, hoprj = 0.0;
    // initialize
    for (int itraj(0); itraj < Ntraj; ++itraj) {
        init_state(state[itraj]);
    }
    // main loop
    vector< vector< complex<double> > > lastevt_save(Ntraj);
    // statistics
    double n0trans = 0.0, n0refl = 0.0, n1trans = 0.0, n1refl = 0.0;
    double px0trans = 0.0, px0refl = 0.0, px1trans = 0.0, px1refl = 0.0;
    double KE = 0.0, PE = 0.0;
    for (int istep(0); istep < Nstep; ++istep) {
        for (int itraj(0); itraj < Ntraj; ++itraj) {
            if (check_end(state[itraj]) == false) {
                // assign last evt
                lastevt = move(lastevt_save[itraj]);
                // calc info
                cal_info(state[itraj][0].real());
                // hopper
                if (enable_hop) {
                    int hopflag = hopper(state[itraj]);
                    switch (hopflag) {
                        case HOP_UP : hopup += 1.0; break;
                        case HOP_DN : hopdn += 1.0; break;
                        case HOP_FR : hopfr += 1.0; break;
                        case HOP_RJ : hoprj += 1.0; break;
                        default : break;
                    }
                }
                // propagate
                rk4.do_step(sys, state[itraj], istep * dt, dt);
                // save last evt
                lastevt_save[itraj] = move(lastevt);
            }
        }
        if (istep % output_step == 0) {
            // data analysis
            // population
            n0trans = n0refl = n1trans = n1refl = 0.0;
            for_each(state.begin(), state.end(), 
                    [&n0trans, &n0refl, &n1trans, &n1refl](const state_t& st) { 
                        int s = static_cast<int>(st[4].real());
                        double vx = st[1].real();
                        if (s == 0) {
                            (vx >= 0.0) ? n0trans += 1.0 : n0refl += 1.0;
                        }
                        else {
                            (vx >= 0.0) ? n1trans += 1.0 : n1refl += 1.0;
                        }
                    });
            // momentum
            px0trans = px0refl = px1trans = px1refl = 0.0;
            for_each(state.begin(), state.end(), 
                    [&px0trans, &px0refl, &px1trans, &px1refl] (const state_t& st) {
                        int s = static_cast<int>(st[4].real());
                        double vx = st[1].real();
                        if (s == 0) {
                            if (vx >= 0.0) {
                                px0trans += mass * vx;
                            }
                            else {
                                px0refl += mass * vx; 
                            }
                        }
                        else {
                            if (vx >= 0.0) {
                                px1trans += mass * vx;
                            }
                            else {
                                px1refl += mass * vx;
                            }
                        }
                    }
                    );
            // energy
            KE = PE = 0.0;
            for_each(state.begin(), state.end(), 
                    [&KE, &PE] (const state_t& st) {
                        double x = st[0].real();
                        double vx = st[1].real();
                        int s = static_cast<int>(st[4].real());
                        vector<double> eva;
                        matrixop::hdiag(cal_H(x), eva);
                        KE += 0.5 * mass * vx * vx;
                        PE += eva[s]; 
                    }
                    );
            // output
            if (istep == 0) {
                // para & header
                ioer::info("# FSSH para: ", " Ntraj = ", Ntraj, " Nstep = ", Nstep, " dt = ", dt, " output_step = ", output_step, " output_mod = ", output_mod,
                            " mass = ", mass, " A = ", A, " B = ", B, " C = ", C, " D = ", D, " E0 = ", E0,
                            " init_x = ", init_x, " init_px = ", init_px, 
                            " sigma_x = ", sigma_x, " sigma_px = ", sigma_px, 
                            " init_s = ", init_s, " xwall_left = ", xwall_left, " xwall_right = ", xwall_right
                        );
                ioer::tabout('#', "t", "n0trans", "n0refl", "n1trans", "n1refl", "px0trans", "px0refl", "px1trans", "px1refl", "Etot");
            }
            if (n0trans > 0.0) { px0trans /= n0trans; }
            if (n0refl > 0.0) { px0refl /= n0refl; }
            if (n1trans > 0.0) { px1trans /= n1trans; }
            if (n1refl > 0.0) { px1refl /= n1refl; }
            n0trans /= Ntraj;
            n0refl /= Ntraj;
            n1trans /= Ntraj;
            n1refl /= Ntraj;
            ioer::tabout('#', istep * dt, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl, (KE + PE) / Ntraj);

            // check end
            bool end_flag = all_of(state.begin(), state.end(), check_end);
            if (end_flag == true) {
                break;
            }
        }

    }
    if (output_mod == "init_px") {
        ioer::tabout(init_px, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl);
    }
    else if (output_mod == "init_s") {
        ioer::tabout(init_s, n0trans, n0refl, n1trans, n1refl, px0trans, px0refl, px1trans, px1refl);
    }
    else {
        misc::crasher::confirm(false, "invalid output mode");
    }
    // hop info
    ioer::info("# hopup = ", hopup, " hopdn = ", hopdn, " hopfr = ", hopfr, " hopfr_rate = ", hopfr / (hopup + hopdn + hopfr));
}

void check_surf() {
    for (double x = xwall_left - 1; x < xwall_right + 1; x += 0.01) {
        cal_info(x);
        ioer::tabout(x, eva, Fx, real(dcx) / 12.0);
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
        fssh();
        ioer::info("# ", timer::toc());
    }
    return 0;
}
