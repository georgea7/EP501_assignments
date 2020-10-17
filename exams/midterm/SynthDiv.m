function [x,R] = SynthDiv(A,b)
%Synthetic Division Algorithm
%x = SynthDiv(A,b)
%A - Polynomial matrix
%b - Divisor
%x - Divided Polynimial
%R - Remainder

n=size(A,2);
x(1)=A(1);
for i=2:n-1
    x(i)=A(i)+b*x(i-1);
end
R=A(n)+b*x(n-1);
end