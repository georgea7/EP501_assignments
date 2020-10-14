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

f=@objfuna;
fprime=@objfuna_deriv;
ygrid=f(x);
plot(x,ygrid);

x0=1;
[R2,it,success]=newton_exact(f,fprime,x0,maxit,tol,true);

