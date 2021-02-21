%runge-kutta method for ode
h = hj
f1 = (-4*R*R*L)/(Sl*Sl*rho)
f2a = 0
f2b = 0
%begin loop
if (R+h)*sin(theta) >= R
  f2a = 2*pi - 8*alpha*(sin(theta)^2)*(R+h)*sin(theta)
  f2b = 4*R*cos(theta)*sqrt(((R+h)*sin(theta))^2 - R^2) - 12*alpha*(((R+h)*sin(theta))^2)*cos(theta) + pi*(2*R*h+h*h)*cos(theta)
else
  f2a = 2*pi*sin(theta)*(R+h)
  f2b = pi*(2Rh+(h*h))*cos(theta)
end
f3 = (12*Uavg*deltad)/(25*h) - (12*Uavg*(deltad^3))/(175*(h^3))
f4 = (((3*(deltad^4))/(175*(h^3)) - (6*(deltad^2))/(25*h))*f1*f2a + (6*Uavg*(deltad^2))/(25*h*h) - (9*Uavg*(deltad^4))/(175*(h^4)))*()

