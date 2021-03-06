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
maxit =100;              %Max Iterations 
tol=1e-3;                %Tolerance
Bessel=@objfunbessel;    %Bessel function
ygrid=Bessel(x);         %Y for graph
verbose = false;         %True to see steps, False to hide steps


%b
x0=x(12);               %Initial guess x_i
x0i1=x(10-2);           %x_i-1
[r1,it,success]=newton_approx(Bessel,x0,x0i1,maxit,tol,verbose);
disp('1st Root of Bessel function');
disp(r1);
disp('Iterations:');
disp(it);
%plot
figure(1)
plot(x,ygrid)
hold on
title('Roots of Bessel Function of Order Zero');
%c
%Finding the first 6 roots of the Bessel function of order zero
k=0;
for i=1:6
    k=k+3;          %Initial guess adjuster
    x0=k;           %Initial guess x_i
    x0i1=k-0.5;     %x_i-1
    [r(i),it,success]=newton_approx(Bessel,x0,x0i1,maxit,tol,verbose);
    figure(1)
    plot(r(i),Bessel(r1),'o','MarkerEdgeColor','k');
    plot(x0,Bessel(x0),'*','MarkerEdgeColor','r');
    plot(x0i1,Bessel(x0i1),'*','MarkerEdgeColor','g');
end
legend('Bessel function: order zero','root','x_i','x_i_-_1');
disp('Roots of Bessel Function of Order Zero');
disp(r);
r_theory=[2.404826, 5.520078, 8.653728, 11.791534, 14.930918, 18.071064];
disp('Roots of Bessel Function of Order Zero- by Vrahatis et al') 
disp(r_theory);
fprintf('Vrahatis, M N, et al. “On the Localization and Computation of Zeros of Bessel Functions.” \n Zeros of Bessel Function, University of Patras,\n thalis.math.upatras.gr/~vrahatis/papers/journals/VrahatisGRZ97_Z_ANGEW_MATH_MECH_77_pp467-475_1997.pdf.\n\n')
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

for j=1:3
    k=-2*1i+j*1i; %Initial condition adjuster
    x0=(j-1)+k;   %Initial condition x_i
    [R3(j),it,success]=newton_exact(f,fprime,x0,maxit,tol,verbose);
end

disp('Roots of x^3 - 3x^2 + 4x - 2 = 0')
disp(R3);

%% Problem 3

%a
disp('3-a)');
f=@objfun2Df;           %function f
g=@objfun2Dg;           %function g
gradf=@grad_objfun2Df;  %f'
gradg=@grad_objfun2Dg;  %g'
x0=0.1;                 %Initial condition x_i
y0=0.1;                 %Initial condition y_i
[rootx(1),rooty(1),it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose);

x0=-0.1;                %Initial condition x_i
y0=-0.5;                %Initial condition y_i
[rootx(2),rooty(2),it,success]=newton2D_exact(f,gradf,g,gradg,x0,y0,maxit,tol,verbose);

disp('     x^2 + y^2 = 2x + y');
disp('(1/4)x^2 + y^2 = 1');
disp('Rootx:');
disp(rootx);
disp('Rooty:');
disp(rooty);

%b
disp('3-b)');
f=@objfun3Df;           %function f
g=@objfun3Dg;           %function g
h=@objfun3Dh;           %function h
gradf=@grad_objfun3Df;  %f'
gradg=@grad_objfun3Dg;  %g'
gradh=@grad_objfun3Dh;  %h'
x0=0.1;                 %initial condition x_i
y0=0.1;                 %initial condition y_i
z0=0.1;                 %initial condition z_i
[rootx(1),rooty(1),rootz(1),it,success]=newton3D_exact(f,gradf,g,gradg,h,gradh,x0,y0,z0,maxit,tol,verbose);

x0=-1;                 %initial condition x_i
y0=-1;                 %initial condition y_i
z0=-1;                 %initial condition z_i
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
