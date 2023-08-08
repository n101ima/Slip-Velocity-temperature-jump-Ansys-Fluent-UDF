# Slip-Velocity-temperature-jump-Ansys-Fluent-UDF
•  The first UDF (inlet_vel) defines a velocity profile for an inlet boundary condition using a modified Poiseuille flow equation that accounts for the Knudsen number. The UDF takes three arguments: the name of the UDF (inlet_vel), the thread pointer (t), and the direction index (i). The UDF loops over all the faces in the face zone and sets the velocity component (F_PROFILE) to the value of the equation evaluated at the face centroid (C_CENTROID). The UDF uses the following constants:

•  H: The height of the channel in meters, which is equal to 5e-5 m.

•  x1: The x-coordinate of the inlet in meters, which is not defined in the UDF.

•  y: The y-coordinate of the face centroid in meters, which is obtained from C_CENTROID.

•  inlet_vel: The maximum velocity at the inlet in meters per second, which is equal to 17 m/s.

•  Kn: The Knudsen number, which is a dimensionless parameter that measures the degree of rarefaction of the gas flow. It is equal to 0.02 in this case.

•  The second UDF (maxwell_slip_velocity_x) implements the Maxwell slip boundary condition for the x-component of velocity at a wall boundary. The Maxwell slip boundary condition accounts for the non-equilibrium effects of rarefied gas flows near solid surfaces, such as velocity slip and temperature jump. The UDF takes three arguments: the name of the UDF (maxwell_slip_velocity_x), the thread pointer (t), and the direction index (i). The UDF loops over all the faces in the face zone and sets the velocity component (F_PROFILE) to the value of the slip velocity plus the thermal creep term. The slip velocity is proportional to the tangential and normal gradients of velocity at the wall, while the thermal creep term is proportional to the tangential gradient of temperature at the wall. The UDF uses several macros and functions from Fluent to access various variables and properties, such as F_C0, THREAD_T0, BOUNDARY_FACE_GEOMETRY, C_U_G, C_V_G, C_W_G, C_MU_L, C_R, and C_T. The UDF also uses some user-defined memory locations (F_UDMI) to store some transformation coefficients that are used to convert from global to local coordinates. The UDF uses the following constants:

•  sigma_square: The squared value of Sigma, which is the molecular diameter of the gas in meters squared. It is equal to 1.35305239e-19 m^2.

•  ambpress: The ambient pressure in pascals, which is equal to 101325 Pa.

•  TMAC: The tangential momentum accommodation coefficient, which is a dimensionless parameter that measures how much momentum is transferred between gas molecules and solid surface. It is equal to 1.0 in this case, which means complete accommodation.

•  ThAC: The thermal accommodation coefficient, which is a dimensionless parameter that measures how much energy is transferred between gas molecules and solid surface. It is also equal to 1.0 in this case, which means complete accommodation.

•  SpHR: The specific heat ratio, which is a dimensionless parameter that measures how much heat capacity changes with temperature. It is equal to 1.4 in this case, which is typical for diatomic gases such as air, oxygen, and nitrogen.

•  UDRLXCOEFF: The under-relaxation coefficient, which is a dimensionless parameter that controls how much the new value of velocity is blended with the old value. It is equal to 0.02 in this case, which means a small amount of under-relaxation.

•  Boltzmann: The Boltzmann constant, which is a physical constant that relates temperature and energy. It is equal to 1.3806505e-23 J/K.

•  PI: The number pi, which is a mathematical constant that appears in many formulas involving circles and spheres. It is approximately equal to 3.14159265358979323846.

•  SQRT_2: The square root of 2, which is another mathematical constant that appears in some formulas involving geometry and algebra. It is approximately equal to 1.41421356237309504880.

•  knodson: The Knudsen number, which is a dimensionless parameter that measures the degree of rarefaction of the gas flow. It is equal to 0.02 in this case.

•  hydraulic_diameter: The hydraulic diameter of the channel in meters, which is equal to 5.0e-5 m.

•  Twall: The temperature of the wall in kelvins, which is equal to 320 K.

•  The third UDF (maxwell_slip_velocity_y) implements the Maxwell slip boundary condition for the y-component of velocity at a wall boundary. It is similar to the second UDF, except that it uses a different transformation coefficient (Coeff1[1]) to calculate the slip velocity and thermal creep term for the y-component. The UDF uses the same arguments, macros, functions, and constants as the second UDF.

# mathematical formulations
•  The first UDF (inlet_vel) uses the following equation for the velocity profile:

$$v(y) = 2 v_{inlet} \left( 1 - \left( \frac{y}{H/2} \right)^2 + 4 Kn \right) / \left( 1 + 8 Kn \right)$$

where $v(y)$ is the velocity at the y-coordinate, $v_{inlet}$ is the maximum velocity at the inlet, $H$ is the height of the channel, and $Kn$ is the Knudsen number.

•  The second UDF (maxwell_slip_velocity_x) uses the following equation for the slip velocity plus the thermal creep term:

$$u_s = \frac{2 - TMAC}{TMAC} \lambda \left( \frac{\partial u}{\partial l} + \frac{\partial u}{\partial n} \right) + 0.75 \frac{\mu}{\rho T} \frac{\partial T}{\partial l}$$

where $u_s$ is the slip velocity, $TMAC$ is the tangential momentum accommodation coefficient, $\lambda$ is the mean free path, $u$ is the velocity component, $l$ and $n$ are the tangential and normal directions to the wall, $\mu$ is the dynamic viscosity, $\rho$ is the density, and $T$ is the temperature.

•  The third UDF (maxwell_slip_velocity_y) uses a similar equation for the slip velocity plus the thermal creep term, but with a different coefficient:

$$v_s = \frac{2 - TMAC}{TMAC} \lambda \left( \frac{\partial v}{\partial l} + \frac{\partial v}{\partial n} \right) + 0.75 \frac{\mu}{\rho T} \frac{\partial T}{\partial l}$$

where $v_s$ is the slip velocity, and $v$ is the velocity component.
