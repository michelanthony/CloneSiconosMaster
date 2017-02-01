/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "NewtonImpactFrictionCohesionNSL.hpp"

#include <iostream>

// Default (private)
NewtonImpactFrictionCohesionNSL::NewtonImpactFrictionCohesionNSL():
  NewtonImpactFrictionNSL(), _c(0.0)
{}

NewtonImpactFrictionCohesionNSL::NewtonImpactFrictionCohesionNSL(unsigned int size):
  NewtonImpactFrictionNSL(size), _c(0.0)
{}

NewtonImpactFrictionCohesionNSL::NewtonImpactFrictionCohesionNSL(double newEn, double newEt, double newMu, double newC, unsigned int newSize):
  NewtonImpactFrictionNSL(newEn, newEt, newMu, newSize), _c(newC)
{}

NewtonImpactFrictionCohesionNSL::~NewtonImpactFrictionCohesionNSL()
{}

bool NewtonImpactFrictionCohesionNSL::isVerified() const
{
  bool res = false;
  // to do
  RuntimeException::selfThrow("NewtonImpactFrictionCohesionNSL:: isVerified, not yet implemented!");
  return res;
}

void NewtonImpactFrictionCohesionNSL::display() const
{
  std::cout << "=== Newton impact-friction non-smooth law with cohesion data display ===" <<std::endl;
  std::cout << " Normal Newton coefficient of restitution: " << _en <<std::endl;
  std::cout << " Tangential Newton coefficient of restitution: " << _et <<std::endl;
  std::cout << " Friction coefficient: " << _mu <<std::endl;
  std::cout << " Cohesion coefficient: " << _c <<std::endl;
  std::cout << "==========================================================" <<std::endl;
}
