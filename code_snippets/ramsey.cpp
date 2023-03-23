// Sample program that creates a ramsey fringe, either circular or linear
//
// Output: either linRamsey.txt or circRamsey.txt, depending on INT_ID choice
// The text file will contain columns freq, zProb, and params
// {W0_VAL, WL_VAL, PHI_VAL, INT_ID}

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include "neutron.hpp"

using namespace std;

// Some initial parameters
const double W_STEP = 0.001; //[rad s^-1]    Step value of w to make ramsey fringes
const double W_START = 180;  //[rad s^-1]    What w to start with
const double W_END = 186;    //[rad s^-1]    What w to end with

const double W0_VAL = 183.247172;  //[rad s^-1]    B0 field strength
const double WL_VAL = 0.732988688; //[rad s^-1]     RF strength
const double PHI_VAL = 0;          //[rad]          RF pulse inital phase

// Time parameters
// Reminder that pulse time cannot have more sig figs than rk_step
const double PULSE_1_TIME = 4.286; // [seconds]
const double PULSE_2_TIME = 4.286; // [seconds]
const double RK_STEP = 0.001;      // [seconds]
const double PRECESS_TIME = 180;   // [seconds]

// Output precision to stdout and file
const int PRECISION = 12;

const double INT_ID = USE_LINEAR_RF; // Type of RF pulse (USE_CIRCULAR_RF or USE_LINEAR_RF)

int main()
{
    vector<double> wOut, zOut, params;
    string filename;
    neutron ucn;
    ofstream outfile;
    double phiVal2, wVal;
    double progress = 0.1;

    int numSteps = (int)((W_END - W_START) / W_STEP);
    cout << "Building ramsey curve (This may take a while)" << endl;
    cout << "0%..." << flush;

    for (int i = 0; i < numSteps; i++)
    {
        // Reset and calculate for new sequence
        ucn.setState({1, 0, 0, 0});
        wVal = (double)i * W_STEP + W_START;
        params = {wVal, W0_VAL, WL_VAL, PHI_VAL, INT_ID};

        // First pulse
        ucn.integrate(PULSE_1_TIME, RK_STEP, params);

        // Free precession period
        ucn.larmorPrecess(PRECESS_TIME, W0_VAL);

        // Pulse 2 has to stay in phase with Pulse 1 while the larmor precession occurs
        phiVal2 = wVal * PULSE_1_TIME + PHI_VAL + wVal * PRECESS_TIME;
        params = {wVal, W0_VAL, WL_VAL, phiVal2, INT_ID};
        ucn.integrate(PULSE_2_TIME, RK_STEP, params);

        // Store solution
        wOut.push_back(wVal);
        zOut.push_back(getZProb(ucn.getState()));

        // Print progress
        if ((double)i / (double)numSteps >= progress)
        {
            cout << progress * 100 << "%..." << flush;
            progress += 0.1;
        }
    }

    // Save output to text
    if (INT_ID == USE_LINEAR_RF)
    {
        filename = "linRamsey.txt";
    }
    else
    {
        filename = "circRamsey.txt";
    }

    cout << "...100%" << endl
         << "Saving output to " << filename << "...";

    outfile.open(filename);
    outfile << "#W0_VAL=" << W0_VAL << ",WL_VAL=" << WL_VAL
            << ",PHI_VAL=" << PHI_VAL << ",INT_ID=" << INT_ID << "\n";
    outfile.precision(PRECISION);
    outfile << "#w,zProb\n";

    for (int i = 0; i < wOut.size(); i++)
    {
        outfile << wOut[i] << "," << zOut[i] << "\n";
    }
    outfile.close();

    cout << "Done!\n";

    return 0;
}
