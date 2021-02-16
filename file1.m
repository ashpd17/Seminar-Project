%---------------------------Physical Parameters--------------------------
dp = 3.5e-3; %Particle diameter dp (m)
rhol = 1029; %Liquid density ?l (kg/m3)
sigma = 7.31e?2; %Liquid surface tension ? (N/m)
mu = 1.01e-3; %Liquid viscosity ? (Pa·s)
Diff = 9.25e-10;%Molecular diffusivity D (m2/s)
thetacon = 70; %Liquid-solid contact angle ?con (°)
%------------------------Operating parameters----------------------------
omega = 59; %Rotational speed ? (rad/s) rad/s 
L = 1.49; %Liquid superficial mass velocity L (kg/(m2·s))
G = 0.03; %Gas superficial mass velocity G (kg/(m2·s))
%------------------------Geometrical parameters--------------------------
R = dp/2 %radius of the catalyst particle, m
Ri = (37/2)*10^-3; %inner radius of the rotor, m
Ro = (79/2)*10^-3; %outer radius of the rotor, m
H = 18*10^-3; %axial height, m
%------------------------Model-------------------------------------------

