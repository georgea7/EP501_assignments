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
n=length(x);                            %size of data

%Plot of noisy data
figure(1);
plot(x,ynoisy,'ko');
xlabel('x');
ylabel('y');
title('Illustrating a Least Square fit')
hold on;

%I discussed with Kaijus Palm on how to do this part.
for N=1:3                               %Polynomials of varying degree
    %1-a)
    %Calculating M matrix for the polynomial
    for i=1:N+1
        M(:,i)=x.^(i-1);                
    end %for
    %Polynomial function y
    a=flipud(inv((M')*M)*(M')*ynoisy);  %Coefficient a
    y=polyval(a,x);                     %Polynomial
    
    %1-b)
    %Error vector and Residual
    df=length(a);                       %number of coefficients
    v=n-df;                             
    error(:,N)=abs(y-ynoisy);
    residual(N)=sum(error,'all');
    
    %1-c)
    Chi_sq(N)=(1/v)*sum(((ynoisy-y).^(2)/(sigmay).^2),'all');
    
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
    'MATLAB built-in Cubic fit','Location','northwest');


