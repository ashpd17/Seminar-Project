clc;
clear;
L = linspace(1.5,8);
n=length(L);
dp = 5.5e-3;
for i=1:n
  k(i) = file1(L(i),dp)*10^4;
end
dp = 3.5e-3;
for i=1:n
  k2(i) = file1(L(i),dp)*10^4;
end
figure
plot(L,k,'k')
hold on 
plot(L,k2,'b')
xlabel('L (kg/m^2 s)')
ylabel('kls (10^{-4} m/s)')
legend('dp=5.5mm','dp=3.5mm')