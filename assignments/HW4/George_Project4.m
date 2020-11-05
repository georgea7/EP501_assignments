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

%Plot of noisy data
figure(1);
plot(x,ynoisy,'ko');
xlabel('x');
ylabel('y');
title('Illustrating a Least Square fit')
hold on;

%a
%I discussed with Kaijus Palm on how to do this part.
for N=1:3                               %Polynomials of varying degree
    %Calculating M matrix for the polynomial
    for i=1:N+1
        M(:,i)=x.^(i-1);                
    end %for
    %Polynomial function y
    a=flipud(inv((M')*M)*(M')*ynoisy);  %Coefficient a
    y=polyval(a,x);                     %Polynomial
    
    %plot
    figure(1)
    hold on
    plot(x,y,'b', 'LineWidth',1.2);
end %for

%Testing with Matlab built-in function
f1=polyfit(x,ynoisy,1);
F1=polyval(f1,x);
f2=polyfit(x,ynoisy,2);
F2=polyval(f2,x);
f3=polyfit(x,ynoisy,3);
F3=polyval(f3,x);
%Plots
figure(1)
hold on
plot(x,F1,'r--','LineWidth',1.2)
plot(x,F2,'r--','LineWidth',1.2)
plot(x,F3,'r--','LineWidth',1.2)
hold off

legend('Data','Linear fit','Quadratic fit','Cubic fit',...
    'MATLAB built-in Linear fit', 'MATLAB built-in Quadratic fit',...
    'MATLAB built-in Cubic fit');

%c
Chi_sq=(1/509)*sum((ynoisy-y).^(2)/(sigmay).^2)
