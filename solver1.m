function [x] = solver1(Vpr,R,thetacon)
   x = pi/2;
   i = 90;
   Vprtest = 1;
   while(i>0 & sign(Vprtest) == 1)
     r1 = (1 - cos(x))/cos(x + thetacon);
     r2 = -(sin(x) - (1-cos(x))*((cos(x+thetacon))^-1 - tan(x+thetacon)));
     Vprtest = 2*pi*(R^3)*(((r1 - r2)^2 + r1^2)*(1-cos(x)) - r1*(r1-r2)*((1-cos(x))*sqrt(1-(((1-cos(x))/r1))^2)+r1*asin((1-cos(x))/r1))-(1-(cos(x))^2)) - Vpr;
     x = x - pi/180;
     i = i - 1;
   end
   if(i== 0)
    error("solver1 error")
   end
end
