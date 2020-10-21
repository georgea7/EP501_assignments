function y=objfun2prime(x)
gamma=5/3;
rho=1.67*10^(-21);  %kg/m^3
p=1.38*10^(-11);    %Pa
B=10^(-9);          %[T]
mu0=4*pi*10^(-7);

Cs=sqrt((gamma*p)/(rho));   %m/s
Ca=sqrt(B^2/(mu0*rho));     %m/s

Cs=Cs*10^(-3);              %km/s
Ca=Ca*10^(-3);              %km/s

A2prime(1) = 6;                                 %x^5
A2prime(2) = 0;                                 %x^4
A2prime(3) = -4*(Ca^2+Cs^2+Ca^2*(cos(pi/4))^2);  %x^3
A2prime(4) = 0;                                 %x^2
A2prime(5) = 2*Ca^2*Cs^2*(cos(pi/4))^2+Ca^4*(cos(pi/4))^2+Ca^2*Cs^2*(cos(pi/4))^2; %x
A2prime(6) = 0;                                 %constant

y=A2prime(1)*x.^5+A2prime(2)*x.^4+A2prime(3)*x.^3+A2prime(4)*x.^2+A2prime(5)*x+A2prime(6);
end %function