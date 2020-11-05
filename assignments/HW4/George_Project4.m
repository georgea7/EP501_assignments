%% Introduction
%Aldous George
%EP 501
%Project 4
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/
clc
clearvars
close all
%% Problem 1
%a
load('test_lsq.mat');
%linear least square fit


y=polyval(a,x);

ytrue=a+b*x;
% y=ytrue+ynoisy;

%Plot of noisy data
figure(1);
plot(x,ynoisy,'ko');
xlabel('x');
ylabel('y');
title('Illustrating a Least Square fit')
hold on;
% plot(x,y,'o','MarkerSize',20);

%Setting up Jacobi


%testing with Matlab built-in function
f1=polyfit(x,ynoisy,1);
F1=polyval(f1,x);
f2=polyfit(x,ynoisy,2);
F2=polyval(f2,x);
f3=polyfit(x,ynoisy,3);
F3=polyval(f3,x);
%Plots
figure(1)
hold on
plot(x,F1)
plot(x,F2)
plot(x,F3)
hold off