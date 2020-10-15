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
verbose = false;

x0=x(12);
x0i1=x(10);
[r1,it,success]=newton_approx(Bessel,x0,x0i1,maxit,tol,verbose);

%b
%plot
figure(1)
plot(x,ygrid)

%c


%% Problem 2

%a
disp('2-a)');
f=@objfuna;
fprime=@objfuna_deriv;

for x0=1:5
    [R2(x0),it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose);
end

disp('Roots of x^5 - 15x^4 + 85x^3 - 225x^2 + 274x - 120 = 0')
disp(R2);

%b
disp('2-b)');
f=@objfunb;
fprime=@objfunb_deriv;
ygrid=f(x);

for j=1:4
    x0=j+1i;
    [R3(j),it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose);
end

disp('Roots of x^3 - 3x^2 + 4x - 2 = 0')
disp(R3);

%% Problem 3

%a
disp('3-a)');
f=@objfun2Df;
g=@objfun2Dg;
gradf=@grad_objfun2Df;
gradg=@grad_objfun2Dg;
x0=0.1;
y0=0.1;
[rootx(1),rooty(1),it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose);

x0=-0.1;
y0=-0.5;
[rootx(2),rooty(2),it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose);

disp('     x^2 + y^2 = 2x + y');
disp('(1/4)x^2 + y^2 = 1');
disp('Rootx:');
disp(rootx);
disp('Rooty:');
disp(rooty);

%b
disp('3-b)');
f=@objfun3Df;
g=@objfun3Dg;
h=@objfun3Dh;
gradf=@grad_objfun3Df;
gradg=@grad_objfun3Dg;
gradh=@grad_objfun3Dh;
x0=0.1;
y0=0.1;
z0=0.1;
[rootx(1),rooty(1),rootz(1),it,success]=newton3D_exact(f,gradf,g,gradg,h,gradh,x0,y0,z0,maxit,tol,verbose);

x0=-1;
y0=-1;
z0=-1;
[rootx(2),rooty(2),rootz(2),it,success]=newton3D_exact(f,gradf,g,gradg,h,gradh,x0,y0,z0,maxit,tol,verbose);



disp(' x^2 + y^2 +  z^2= 6');
disp(' x^2 - y^2 + 2z^2= 2');
disp('2x^2 + y^2 -  z^2= 3');
disp('X Root:');
disp(rootx);
disp('Y Root:');
disp(rooty);
disp('Z Root:');
disp(rootz);


