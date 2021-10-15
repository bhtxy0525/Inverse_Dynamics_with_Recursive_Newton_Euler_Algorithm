# Inverse Dynamics with recursive Newton Euler Algorithm
Developed by: Xinyang Tian (2021).

Supervisor: Haolin Ma.

Platform: Matlab 2019b.

In this work, the dynamic model of a 7-DOF robot manipulator has been established. Here are some notes about this script:

- Using Newton-Euler algorithm to calculate the inverse dynamic model of serial robot manipulator (only rotational joints) with modified DH parameters.
- No external force/torque.

# Implementation 
The code is implemented in MATLAB R2019b. Also, other versions of Matlab are available. Before running this script, pleas enter the input parameters: q, qd, and qdd. Cubic or quintic polynomials are recommended for planning. As a consequence, the output parameters are seven torques corresponding to seven joints of the robot manipulator.

# References
