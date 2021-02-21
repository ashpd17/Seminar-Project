function [Vprtest] = solver3(Vpr,R,thetacon,x)
  r1 = (1 - cos(x))/cos(x + thetacon);
  r2 = -(sin(x) - (1-cos(x))*((cos(x+thetacon))^-1 - tan(x+thetacon)));
  Vprtest = 2*pi*(R^3)*(((r1 - r2)^2 + r1^2)*(1-cos(x)) - r1*(r1-r2)*((1-cos(x))*sqrt(1-(((1-cos(x))/r1))^2)+r1*asin((1-cos(x))/r1))-(1-(cos(x))^2)) - Vpr;
end
