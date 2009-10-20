/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */


// =============================== Robot arm sample (HuMAnsPa10) ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, EventDriven, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 1.9;                   // final computation time
    double h = 0.005;                // time step
    double e = 0.9;                  // restit. coef. for impact on the ground.
    double e2 = 0.0;                 // restit. coef for angular stops impacts.

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS

    // --- DS: robot arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See Robot.fig for more details)
    //


    // Initial position (angles in radian)
    SimpleVector q0(nDof), v0(nDof);
    q0(0) = 0.05;
    q0(1) = 0.05;

    SP::LagrangianDS arm(new LagrangianDS(q0, v0));

    // external plug-in
    arm->setComputeMassFunction("RobotPlugin.so", "mass");
    arm->setComputeNNLFunction("RobotPlugin.so", "NNL");
    arm->setComputeJacobianNNLFunction(1, "RobotPlugin.so", "jacobianVNNL");
    arm->setComputeJacobianNNLFunction(0, "RobotPlugin.so", "jacobianQNNL");

    allDS.insert(arm);

    // -------------------
    // --- Interactions---
    // -------------------

    // Two interactions:
    //  - one with Lagrangian non linear relation to define contact with ground
    //  - the other to define angles limitations (articular stops), with lagrangian linear relation
    //  Both with newton impact ns laws.

    InteractionsSet allInteractions; // The set of all interactions.

    // -- relations --

    // => arm-floor relation
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    string G = "RobotPlugin:G2";
    SP::Relation relation(new LagrangianScleronomousR("RobotPlugin:h2", G));
    SP::Interaction inter(new Interaction("floor-arm", allDS, 0, 2, nslaw, relation));

    // => angular stops

    //     SimpleMatrix H(6,3);
    //     SimpleVector b(6);
    //     H.zero();
    //     H(0,0) =-1;
    //     H(1,0) =1;
    //     H(2,1) =-1;
    //     H(3,1) =1;
    //     H(4,2) =-1;
    //     H(5,2) =1;

    //     b(0) = 1.7;
    //     b(1) = 1.7;
    //     b(2) = 0.3;
    //     b(3) = 0.3;
    //     b(4) = 3.14;
    //     b(5) = 3.14;
    double lim0 = 1.6;
    double lim1 = 3.1;  // -lim_i <= q[i] <= lim_i
    SimpleMatrix H(4, 3);
    SimpleVector b(4);
    H.zero();

    H(0, 0) = -1;
    H(1, 0) = 1;
    H(2, 1) = -1;
    H(3, 1) = 1;

    b(0) = lim0;
    b(1) = lim0;
    b(2) = lim1;
    b(3) = lim1;

    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation2(new LagrangianLinearTIR(H, b));
    SP::Interaction inter2(new Interaction("floor-arm2", allDS, 1, 4, nslaw2, relation2));

    allInteractions.insert(inter);
    allInteractions.insert(inter2);

    // -------------
    // --- Model ---
    // -------------

    SP::Model Robot(new Model(t0, T, allDS, allInteractions));

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::EventDriven s(new EventDriven(t));

    // -- OneStepIntegrators --
    SP::Lsodar OSI(new Lsodar(arm));
    s->recordIntegrator(OSI);
    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 20001; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] =  0.005; // Tolerance
    string solverName = "PGS" ;
    SP::NonSmoothSolver mySolver(new NonSmoothSolver(solverName, iparam, dparam));
    SP::OneStepNSProblem impact(new LCP(mySolver, "impact"));
    SP::OneStepNSProblem acceleration(new LCP(mySolver, "acceleration"));
    s->recordNonSmoothProblem(impact);
    s->recordNonSmoothProblem(acceleration);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================  dataPlot(k,7) = (*inter->y(0))(0);


    // ================================= Computation =================================

    // --- Simulation initialization ---
    Robot->initialize(s);
    cout << "End of model initialisation" << endl;

    int k = 0;
    int N = 10630;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    SP::SiconosVector q = arm->q();
    SP::SiconosVector vel = arm->velocity();
    SP::SiconosVector y = inter->y(0);
    SP::SiconosVector yDot = inter->y(1);
    // When a non-smooth event occurs, pre-impact values are saved in memory vectors at pos. 1:
    SP::SiconosVector qMem = arm->getQMemoryPtr()->getSiconosVector(1);
    SP::SiconosVector velMem = arm->getVelocityMemoryPtr()->getSiconosVector(1);
    SP::SiconosVector yMem = inter->getYOldPtr(0);
    SP::SiconosVector yDotMem = inter->getYOldPtr(1);

    dataPlot(k, 0) =  Robot->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*vel)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*vel)(1);
    dataPlot(k, 5) = (*q)(2);
    dataPlot(k, 6) = (*vel)(2);
    dataPlot(k, 7) = (*y)(0);
    dataPlot(k, 8) = (*y)(1);
    dataPlot(k, 9) = (*yDot)(0);
    dataPlot(k, 10) = (*yDot)(1);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    boost::timer boostTimer;
    boostTimer.restart();

    unsigned int numberOfEvent = 0 ;
    SP::EventsManager eventsManager = s->eventsManager();
    bool nonSmooth = false;
    while (s->nextTime() < T)
    {
      // get current time step
      k++;
      s->advanceToEvent();
      if (eventsManager->nextEvent()->getType() == 2)
        nonSmooth = true;

      s->processEvents();
      // If the treated event is non smooth, we get the pre-impact state.
      if (nonSmooth)
      {
        dataPlot(k, 0) =  s->startingTime();
        dataPlot(k, 1) = (*qMem)(0);
        dataPlot(k, 2) = (*velMem)(0);
        dataPlot(k, 3) = (*qMem)(1);
        dataPlot(k, 4) = (*velMem)(1);
        dataPlot(k, 5) = (*qMem)(2);
        dataPlot(k, 6) = (*velMem)(2);
        dataPlot(k, 7) = (*yMem)(0);
        dataPlot(k, 8) = (*yMem)(1);
        dataPlot(k, 9) = (*yDotMem)(0);
        dataPlot(k, 10) = (*yDotMem)(1);


        k++;
      }
      dataPlot(k, 0) =  s->startingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*vel)(0);
      dataPlot(k, 3) = (*q)(1);
      dataPlot(k, 4) = (*vel)(1);
      dataPlot(k, 5) = (*q)(2);
      dataPlot(k, 6) = (*vel)(2);
      dataPlot(k, 7) = (*y)(0);
      dataPlot(k, 8) = (*y)(1);
      dataPlot(k, 9) = (*yDot)(0);
      dataPlot(k, 10) = (*yDot)(1);
      numberOfEvent++;
      //  cout << k << endl;
      //  if (k==N) break;
    }

    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl << endl;
    cout << "Computation Time: " << boostTimer.elapsed()  << endl;

    cout << endl << "Output writing ..." << endl;
    // --- Output files ---
    ioMatrix out("result.dat", "ascii");
    out.write(dataPlot, "noDim");
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
}
