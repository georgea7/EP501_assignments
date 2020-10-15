%% Introduction
%Aldous George
%EP 501
%Project 3
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/blob/master/nonlinear_eqns
clc
clearvars
close all
%% Problem 1

%setup
x=linspace(0,20);        %radial independent variable
maxit =100;
tol=1e-3;
Bessel=@objfunbessel;
ygrid=Bessel(x);

x0=x(12);
x0i1=x(10);
[r1,it,success]=newton_approx(Bessel,x0,x0i1,maxit,tol,true);

%b
%plot
figure(1)
plot(x,ygrid)

%c


%% Problem 2

%a
f=@objfuna;
fprime=@objfuna_deriv;
ygrid=f(x);
plot(x,ygrid);

for x0=1:5
    [R2(x0),it,success]=newton_exact(f,fprime,x0,maxit,tol,true);
end

%b
f=@objfunb;
fprime=@objfunb_deriv;
ygrid=f(x);
plot(x,ygrid);

for j=1:4
    x0=j+1i;
    [R3(j),it,success]=newton_exact(f,fprime,x0,maxit,tol,true);
end

%% Problem 3

%a
f=@objfun2Df;
g=@objfun2Dg;
gradf=@grad_objfun2Df;
gradg=@grad_objfun2Dg;
x0=0.1;
y0=0.1;
[rootx,rooty,it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,true);

%b



