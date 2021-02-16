%---------------------------Physical Parameters--------------------------
dp = 3.5e-3; %Particle diameter dp (m)
rhol = 1029; %Liquid density ?l (kg/m3)
sigma = 7.31e-2; %Liquid surface tension ? (N/m)
mu = 1.01e-3; %Liquid viscosity ? (Pa·s)
Diff = 9.25e-10;%Molecular diffusivity D (m2/s)
thetacon = 70; %Liquid-solid contact angle con (°)
%------------------------Operating parameters----------------------------
omega = 59; %Rotational speed ? (rad/s) rad/s 
L = 1.49; %Liquid superficial mass velocity L (kg/(m2·s))
G = 0.03; %Gas superficial mass velocity G (kg/(m2·s))
%------------------------Geometrical parameters--------------------------
R = dp/2; %radius of the catalyst particle, m
Ri = (37/2)*10^-3; %inner radius of the rotor, m
Ro = (79/2)*10^-3; %outer radius of the rotor, m
H = 18*10^-3; %axial height, m
%------------------------Model-------------------------------------------
i = 1;
N = 10; %assumed 
for i=1:N
  Rri = Ri + (2*i-1)*R;
  V(i) = pi*((Rri + R)^2 - (Rri-R)^2)*H;
  epb = 1 - pi/6;
  Eo = rhol*(omega^2)*Rri*(dp^2)*(epb^2)/(sigma*(1-epb));
  esext = 1/(20 + 0.9*Eo);
  Nc = 3.42/epb - 1.18;
  Vpr = pi*(dp^3)*esext/(3*Nc*(1-epb));
  %calculate beta = (Vpr,R,thetacon) separate program Appendix D;
  beta = 20; %assumed 
  theta = 20; %assumed take radian
  Q = 10; %calculate Q assumed 
  hj = ((3*mu*Q)/(2*pi*R*rhol*(omega^2)*Rri*(sin(theta))^2));
  alpha = acos(R/((R+h)*sin(theta)));
  if ((R+hj)*sin(theta) >= R)
    Sl = (4*tan(alpha)-pi)*R^2 + (pi - 4*alpha*(sin(theta))^2)*(R+hj)^2;
  else
     Sl = pi*(2*R*h+h^2);
  end
  Uavg = 4*(R^2)*L/(rhol*Sl*sin(theta));
  %Y assumed ranges from 0 to h
  %U = -(4*Uavg/(5*h^3))*Y^3 + ((12*Uavg)/(5*h))*Y; Not required ig
  Sc = mu/(rhol*Diff);
  Xbeta = 0.00001 %assumed the length of the arc corresponding to beta 
  Ub = (8/5)*Uavg;
  delta = 4.64*(mu*Xbeta/(rhol*Ub))^(1/2);
  deltad = delta/(Sc^(1/3));
  %solve ode in theta deltad eqn C1,C2...
  %solve trapz using eqn 24 
  klsi(i) = 10;
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



