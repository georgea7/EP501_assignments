%% Introduction
%Aldous George
%EP 501
%Project 6
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/
clc
clearvars
close all
%% Problem 1
%Initiation
a = 0.01;           %(m)
l = a/5;            %(m)
xprime = -9*a/10;   %(m)
xdprime = 9*a/10;   %(m)
Eps0 = 8.854*10^-12;%(F/m)

n=100;
x=linspace(-a,a,n);

%1-a
disp('1-a');
Eps = Eps0*10*(tanh((x-xprime)/l)-tanh((x-xdprime)/l));

%plot
figure(1)
plot(x,Eps);
set(gca,'FontSize',15);
xlabel('x');
ylabel('\epsilon');
title('Variation of Dielectric Function')




%% Problem 2

%Initiation
m = 1.67*10^(-27);            %(kg)
q = 1.6*10^(-19);             %(C)
B = 50000;                    %(nT)
vx0 = 1;                      %(km/s)
vy0 = 1;                      %(km/s)