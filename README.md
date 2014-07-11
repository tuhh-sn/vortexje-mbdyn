vortexje-mbdyn
==============

Sönke Neumann	<soenke.neumann@tuhh.de>

Interface towards MBDyn within Vortexje

This repository offers an example of the integration of MBDyn's kinematics into the panel code Vortexje. Both Open Source codes are accessible under http://www.aero.polimi.it/mbdyn/ and http://vortexje.org/.

The example refers to the examples/vawt.cpp scenario of Vortexje. The model represents a two-bladed vertical-axis wind turbine. The blades are fixed to one rotating point which interferes with MBDyn. Within MBDyn two point masses simulate the kinematic behavior of the rotor. 
The coupling method is strong coupling and can be tuned with an underrelaxation factor.

