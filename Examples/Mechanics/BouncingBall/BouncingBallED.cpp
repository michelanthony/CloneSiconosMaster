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

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model without XML input.
  Simulation with an Event-Driven scheme.
*/

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 8.5;                   // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SimpleVector q0(new SimpleVector(nDof));
    SP::SimpleVector v0(new SimpleVector(nDof));
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;
    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));
    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SP::SiconosMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 1.0;
    SP::NonSmoothLaw  nslaw0(new NewtonImpactNSL(e));
    SP::Relation relation0(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(1, nslaw0, relation0));

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    // -------------
    // --- Model ---
    // -------------

    SP::Model bouncingBall(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->link(inter, ball);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- (1) OneStepIntegrators --
    SP::OneStepIntegrator OSI(new Lsodar(ball));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) Non smooth problem --
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::EventDriven s(new EventDriven(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_ACCELERATION);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Simulation initialisation ..." << endl << endl;
    s->setPrintStat(true);
    bouncingBall->initialize(s);

    int N = 1854; // Number of saved points: depends on the number of events ...

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);
    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    //   SiconosVector * y = bouncingBall->nonSmoothDynamicalSystem()->interaction(0)->y(0);

    SP::EventsManager eventsManager = s->eventsManager();

    // For the initial time step:
    // time

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    bool nonSmooth = false;
    unsigned int numberOfEvent = 0 ;
    int k = 0;
    boost::progress_display show_progress(N);
    while (s->nextTime() < T && k < N)
    {
      s->advanceToEvent();
      if (eventsManager->nextEvent()->getType() == 2)
        nonSmooth = true;

      s->processEvents();
      // If the treated event is non smooth, the pre-impact state has been solved in memory vectors during process.
      if (nonSmooth)
      {
        dataPlot(k, 0) = s->startingTime();
        dataPlot(k, 1) = (*ball->qMemory()->getSiconosVector(1))(0);
        dataPlot(k, 2) = (*ball->velocityMemory()->getSiconosVector(1))(0);
        k++;
        nonSmooth = false;
        ++show_progress;
      }
      dataPlot(k, 0) = s->startingTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      ++k;
      ++numberOfEvent;
      ++show_progress;
    }

    // --- Output files ---
    cout << endl;
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl << endl;
    cout << "====> Output file writing ..." << endl << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught." << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}
