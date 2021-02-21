% runge-kutta method for ode
h = hj 
f2a = 0
f2b = 0
%begin loop
if (R + h)*sin(theta) >= R
  f2a = 2*pi - 8*alpha*(sin(theta)^2)*(R + h)*sin(theta)
  f2b = 4*R*cos(theta)*sqrt(((R + h)*sin(theta))^2 - R^2) - 12*alpha*(((R + h)*sin(theta))^2)*cos(theta) + pi*(2*R*h + h*h)*cos(theta)
