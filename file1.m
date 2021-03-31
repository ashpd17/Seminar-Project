clc;
clear;
global mu R rhol Rri L Diff Q omega;
%---------------------------Physical Parameters--------------------------
dp = 3.5e-3; %Particle diameter dp (m)
rhol = 1029; %Liquid density ?l (kg/m3)
sigma = 7.31e-2; %Liquid surface tension ? (N/m)
mu = 1.01e-3; %Liquid viscosity ? (Pa·s)
Diff = 9.25e-10;%Molecular diffusivity D (m2/s)
thetacon = 70*pi/180; %Liquid-solid contact angle con (°)
g = 9.8; %acceleration due to gravity m/s2
%------------------------Operating parameters----------------------------
omega = 176; %Rotational speed ? (rad/s) rad/s 
L = 1.49; %Liquid superficial mass velocity L (kg/(m2·s))
G = 0.03; %Gas superficial mass velocity G (kg/(m2·s))
%------------------------Geometrical parameters--------------------------
R = dp/2; %radius of the catalyst particle, m
Ri = (37/2)*10^-3; %inner radius of the rotor, m
Ro = (79/2)*10^-3; %outer radius of the rotor, m
H = 18*10^-3; %axial height, m
%------------------------Model-------------------------------------------
i = 1;
N = 1; %assumed 
for i=1:N
  Rri = Ri + (2*i-1)*R;
  V(i) = pi*((Rri + R)^2 - (Rri-R)^2)*H;
  epb = 1 - pi/6;
  Eo = rhol*(omega^2)*Rri*(dp^2)*(epb^2)/(sigma*(1-epb));
  esext = 1/(20 + 0.9*Eo);
  Nc = 3.42/epb - 1.18;
  Vpr = pi*(dp^3)*esext/(3*Nc*(1-epb));
  Q = L*(2*pi*R*H)/rhol; %calculate Q assumed 
  specificsolver1 = @(x) solver1(Vpr,R,thetacon,x);
  betab = fzero(specificsolver1,[pi/360,pi/2]); %made 1 change
  Xbetab = R*betab; %assumed the length of the arc corresponding to betab 
  hj = ((3*mu*Q)/(2*pi*R*rhol*(omega^2)*Rri*(sin(betab))^2));
  alpha = acos(R/((R+hj)*sin(betab)));
  if ((R+hj)*sin(betab) >= R)
    Sl = (4*tan(alpha)-pi)*R^2 + (pi - 4*alpha*(sin(betab))^2)*(R+hj)^2;
  else
     Sl = pi*(2*R*hj+hj^2);
  end
  Uavg = 4*(R^2)*L/(rhol*Sl*sin(betab));
  Ub = (8/5)*Uavg;
 
  thetaspan = [betab pi-betab];
  delta = 4.64*(mu*Xbetab/(rhol*Ub))^(1/2);
  Sc = mu/(rhol*Diff);
  deltad = delta/Sc^(1/3);
  [theta, deltad] = ode45(@solver2, thetaspan, deltad);
  %solve trapz using eqn 24
  Di = 2*Rri;
  Re = (L*dp)/mu;
  Ga = (omega^2)*Rri*(dp^3)*(rhol^2)*(mu)^-2;
  Ka = (mu^4)*g/((sigma^3)*rhol);
  f = 0.408*((dp/Di)^(-0.231))*(Re^(0.072))*((Ga)^0.036)*((Ka)^0.022);
  term = (0.75*f*Diff)*sin(theta)./deltad;
  klsi(i) = trapz(term,theta);
  plot(theta,deltad)
end
Vsum = 0;
kls = 0;
for i=1:N 
  Vsum = Vsum + V(i);
end
for i=1:N
  kls = kls + klsi(i)*V(i);
end
kls = kls/Vsum
