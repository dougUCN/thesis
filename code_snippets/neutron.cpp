#include <vector>
#include <cmath>
#include <iostream>
#include "neutron.hpp"

using namespace std;

void neutron::setState(const vector<double>& ket)
{
    if (ket.size() != NUM_EQ){
        cout << "neutron::setState ket_size != " << NUM_EQ <<endl;
        exit(-1);
    }
    _u = ket;
}

vector<double> neutron::getState()
{
    return _u;
}

void neutron::larmorPrecess(double precTime, double w0)
{
    vector<double> _uEnd(NUM_EQ);
    double x = precTime * w0/2;
    _uEnd[0] = _u[0] * cos(x) + _u[1] * sin(x);
    _uEnd[1] = _u[1] * cos(x) - _u[0] * sin(x);
    _uEnd[2] = _u[2] * cos(x) - _u[3] * sin(x);
    _uEnd[3] = _u[3] * cos(x) + _u[2] * sin(x);
    _u = _uEnd;
}

void neutron::rkStep(const double t, const double dt, const vector<double>& params)
// RK4 integration step
{
    vector<double> f0 = derivs( t, _u, params );
    double t1, t2, t3;
    vector<double> u1(NUM_EQ);
    vector<double> u2(NUM_EQ);
    vector<double> u3(NUM_EQ);
    vector<double> f1, f2, f3;

    t1 = t + dt/2.0;
    for (int i = 0; i < NUM_EQ; i++)
        u1[i] = _u[i] + dt * f0[i] / 2.0;
    f1 = derivs( t1, u1, params );

    t2 = t + dt / 2.0;
    for (int i = 0; i < NUM_EQ; i++)
        u2[i] = _u[i] + dt * f1[i] / 2.0;
    f2 = derivs( t2, u2, params );

    t3 = t + dt;
    for (int i = 0; i < NUM_EQ; i++)
        u3[i] = _u[i] + dt * f2[i];
    f3 = derivs( t3, u3, params );

    for (int i = 0; i < NUM_EQ; i++)
        _u[i] += ( dt / 6.0 ) * (f0[i] + 2*f1[i] + 2*f2[i] + f3[i]);
}

void neutron::integrate(const double time, const double dt, 
                            const vector<double>& params)
{
    int t = 0;
    while( (double)t * dt < time)
    {
        rkStep( (double)t*dt, dt, params );
        t++;
    }
}

void neutron::integrate(const double time, const double dt, const vector<double>& params,
    vector<double>& tOut, vector<double>& xOut, vector<double>& yOut, vector<double>& zOut)
{
    int t = 0;
    tOut.push_back(0);
    xOut.push_back(getXProb(this->getState()));
    yOut.push_back(getYProb(this->getState()));
    zOut.push_back(getZProb(this->getState()));
    while( (double)t * dt < time)
    {
        rkStep( (double)t*dt, dt, params );
        t++;
        tOut.push_back((double)t*dt);
        xOut.push_back(getXProb(this->getState()));
        yOut.push_back(getYProb(this->getState()));
        zOut.push_back(getZProb(this->getState()));
    }
}

vector<double> neutron::derivs( const double t, const vector<double>& u, 
                                const vector<double>& params)
// vector<double> params should be in the form of {w, w0, wRF, phi, INT_ID}
// w is the driving RF frequency in rad/s
// w0 is the strength of the applied B0 field in rad/s
// wRF is the strength of the linear/circular RF field in rad/s
// INT_ID is either LIN_ID (linear RF) or CIRC_ID (circular RF)
//
// Modified right hand side of eq A.1 - A.4 in Daniel May's nEDM thesis
// using eq 3.28, 3.29 as a basis
// u[0] = Re(a), u[1] = Im(a), u[2] = Re(b), u(3) = Im(b)
{
    vector<double> dudt(NUM_EQ);
    double x = params[0] * t + params[3];
    dudt[0] = 0.5*(params[1]*u[1] + params[2]*cos(x)*u[3]) - params[4]/2*params[2]*u[2]*sin(x);
    dudt[1] = 0.5*(-params[1]*u[0] - params[2]*cos(x)*u[2]) - params[4]/2*params[2]*u[3]*sin(x);
    dudt[2] = 0.5*(-params[1]*u[3] + params[2]*cos(x)*u[1]) + params[4]/2*params[2]*u[0]*sin(x);
    dudt[3] = 0.5*(params[1]*u[2] - params[2]*cos(x)*u[0]) + params[4]/2*params[2]*u[1]*sin(x);
    return dudt;
}

double getXProb(const vector<double>& u)  // Odds of measuring spin up along x
{
    if (u.size() != NUM_EQ){
        cout << "getXProb ket_size != " << NUM_EQ <<endl;
        exit(-1);
    }
    return 0.5 + u[0]*u[2] + u[1]*u[3];
}

double getYProb(const vector<double>& u)  // Odds of measuring spin up along y
{
    if (u.size() != NUM_EQ){
        cout << "getYProb ket_size != " << NUM_EQ <<endl;
        exit(-1);
    }
    return 0.5 + u[1]*u[2] - u[3]*u[0];
}

double getZProb(const vector<double>& u)  // Odds of measuring spin up along z
{
    if (u.size() != NUM_EQ){
        cout << "getZProb ket_size != " << NUM_EQ <<endl;
        exit(-1);
    }
    return u[0]*u[0] + u[1]*u[1];
}
