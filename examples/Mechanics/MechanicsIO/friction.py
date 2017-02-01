#!/usr/bin/env python

#
# Several objects on increasing inclined planes
#

from siconos.mechanics.collision.tools import Contactor
from siconos.io.mechanics_io import Hdf5
import siconos.numerics as Numerics
import siconos.kernel as Kernel
import numpy as np

# Creation of the hdf5 file for input/output
with Hdf5() as io:

    # Definition of a cube
    io.addPrimitiveShape('Cube', 'Box', [1,1,1])

    # Definition of the ground shape
    io.addPrimitiveShape('Ground', 'Box', (2, 10, .5))

    # Definition of a non smooth law. As no group ids are specified it
    # is between contactors of group id 0.
    io.addNewtonImpactFrictionNSL('contact', mu=0.5, e=0.0)

    # Define an array of cubes on surfaces of increasing angle
    for i,angle in enumerate(np.linspace(0,np.pi,20)):

        z = 0.74  # Already small amount of penetration to ensure
                  # contacts are generated
        cosa = np.cos(angle)
        sina = np.sin(angle)

        io.addObject('cube%d'%i, [Contactor('Cube')],
                     translation=[-i*3, z*sina, z*cosa],
                     orientation=[(-1,0,0), angle],
                     mass=1)

        io.addObject('ground%d'%i, [Contactor('Ground')],
                     translation=[-i*3, 0, 0],
                     orientation=[(-1,0,0), angle])

# Run the simulation from the inputs previously defined and add
# results to the hdf5 file. The visualisation of the output may be done
# with the vview command.
with Hdf5(mode='r+') as io:

    # By default earth gravity is applied and the units are those
    # of the International System of Units.
    # Because of fixed collision margins used in the collision detection,
    # sizes of small objects may need to be expressed in cm or mm.
    io.run(t0=0,
           T=10,
           h=0.01,
           multipoints_iterations=True,
           theta=0.50001,
           Newton_max_iter=1,
           solver=Numerics.SICONOS_FRICTION_3D_NSGS,
           itermax=100000,
           tolerance=1e-8)
