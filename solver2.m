%initial conditions
hi = ((3*mu*Q)/(2*pi*R*rhol*(omega^2)*Rri*(sin(betab))^2));
alphai = acos(R/((R+hi)*sin(betab)));
if ((R+hi)*sin(beta) >= R)
  Sli = (4*tan(alphai)-pi)*R^2 + (pi - 4*alphai*(sin(beta))^2)*(R+hi)^2;
else
  Sli = pi*(2*R*hi+hi^2);
end
Uavgi = 4*(R^2)*L/(rhol*Sli*sin(betab));
Sc = mu/(rhol*Diff);
Xbetab = 0.00001; %assumed the length of the arc corresponding to betab
Ub = (8/5)*Uavgi; %Uavgi comes from file 1, will edit it
deltai = 4.64*(mu*Xbetab/(rhol*Ub))^(1/2);
deltadi = deltai/(Sc^(1/3));
thetaspan = [betab pi-betab];
%ODE solver part
h = ((3*mu*Q)/(2*pi*R*rhol*(omega^2)*Rri*(sin(theta))^2));
alpha = acos(R/((R+h)*sin(theta)));
if ((R+h)*sin(theta) >= R)
  Sl = (4*tan(alpha)-pi)*R^2 + (pi - 4*alpha*(sin(theta))^2)*(R+h)^2;
else
  Sl = pi*(2*R*h+h^2);
end
Uavg = 4*(R^2)*L/(rhol*Sl*sin(theta));
f1 = (-4*R*R*L)/(Sl*Sl*rho);
if (R+h)*sin(theta) >= R
  f2a = 2*pi - 8*alpha*(sin(theta)^2)*(R+h)*sin(theta);
  f2b = 4*R*cos(theta)*sqrt(((R+h)*sin(theta))^2 - R^2) - 12*alpha*(((R+h)*sin(theta))^2)*cos(theta) + pi*(2*R*h+h*h)*cos(theta);
else
  f2a = 2*pi*sin(theta)*(R+h);
  f2b = pi*(2Rh+(h*h))*cos(theta);
end
f3 = (12*Uavg*deltad)/(25*h) - (12*Uavg*(deltad^3))/(175*(h^3));
f4 = ((((3*(deltad^4))/(175*(h^3)) - (6*(deltad^2))/(25*h))*f1*f2a + (6*Uavg*(deltad^2))/(25*h*h) - (9*Uavg*(deltad^4))/(175*(h^4)))*((-2/3)*(((3*mu*Q)/(2*pi*R*rho*omega*omega*Rri))^0.333)*((sin(theta))^(-1.667))*cos(theta))) + ((3*(deltad^4))/(175*h*h*h) - (6*(deltad^4))/(25*h))*f1*f2b + (3*Uavg*(deltad^4))/(175*h*h*h*tan(theta)) - (9*Uavg*deltad*deltad)/(25*h*tan(theta)) + (3*R*Diff)/(2*deltad);
[theta, deltad] = ode45(@(theta, deltad) f4/f3, thetaspan, deltadi);
