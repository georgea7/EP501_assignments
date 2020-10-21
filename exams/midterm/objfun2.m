function y=objfun2(x)
gamma=5/3;
rho=1.67*10^(-21);  %kg/m^3
p=1.38*10^(-11);    %Pa
B=10^(-9);          %[T]
mu0=4*pi*10^(-7);

Cs=sqrt((gamma*p)/(rho));   %m/s
Ca=sqrt(B^2/(mu0*rho));     %m/s

Cs=Cs*10^(-3);              %km/s
Ca=Ca*10^(-3);              %km/s

A2(1)= 1;                                   %x^6
A2(2)= 0;                                   %x^5
A2(3)= -(Ca^2+Cs^2+Ca^2*(cos(pi/4))^2);     %x^4
A2(4)= 0;                                   %x^3
A2(5)= Ca^2*Cs^2*(cos(pi/4))^2+Ca^4*(cos(pi/4))^2+Ca^2*Cs^2*(cos(pi/4))^2; %x^2
A2(6)= 0;                                   %x
A2(7)= -(Ca^4*Cs^2*(cos(pi/4))^4);          %constant

y=A2(1)*x.^6+A2(2)*x.^5+A2(3)*x.^4+A2(4)*x.^3+A2(5)*x.^2+A2(6)*x+A2(7);
end %function