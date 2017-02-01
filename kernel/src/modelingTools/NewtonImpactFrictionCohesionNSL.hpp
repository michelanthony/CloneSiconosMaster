/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
/*! \file NewtonImpactFrictionCohesionNSL.hpp
  Newton-Impact Non-Smooth Law
*/

#ifndef NEWTONIMPACTFRICTIONCOHESIONNSLAW_H
#define NEWTONIMPACTFRICTIONCOHESIONNSLAW_H

#include "NewtonImpactFrictionNSL.hpp"

/** Newton Impact-Friction Non Smooth Law with Cohesion
 *
 */
class NewtonImpactFrictionCohesionNSL : public NewtonImpactFrictionNSL
{

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonImpactFrictionCohesionNSL);

  /** cohesion coefficient */
  double _c;

  /** default constructor
   */
  NewtonImpactFrictionCohesionNSL();

public:

  /** basic constructor
   *  \param size size of the ns law
   */
  NewtonImpactFrictionCohesionNSL(unsigned int size);

  /** constructor with the value of the NewtonImpactFrictionCohesionNSL attributes
   *  \param en double : normal e coefficient
   *  \param et double : tangent e coefficient
   *  \param mu double : friction coefficient
   *  \param size unsigned int: size of the ns law
   */
  NewtonImpactFrictionCohesionNSL(double en, double et, double mu, double c,
                          unsigned int size);

  /** Destructor */
  ~NewtonImpactFrictionCohesionNSL();

  /** check the ns law to see if it is verified
   *  \return a boolean value whioch determines if the NS Law is verified
   */
  bool isVerified(void) const;

  // GETTERS/SETTERS

  /** getter of c
   * \return the value of c
   */
  inline double c() const
  {
    return _c;
  };

  /** setter of c
   * \param newVal a double to set c
   */
  inline void setCohesion(double newVal)
  {
    _c = newVal;
  };

  // OTHER FUNCTIONS

  /** print the data to the screen
   */
  void display() const;

  /** Visitors hook
   */
  ACCEPT_STD_VISITORS();

};
DEFINE_SPTR(NewtonImpactFrictionCohesionNSL)
#endif // NEWTONIMPACTFRICTIONCOHESIONNSLAW_H
