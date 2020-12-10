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
    error(:,N)=(y-ynoisy);
    residual(N)=sum(error,'all');
    
    %1-c)
    Chi_sq(N)=(1/v)*sum(((ynoisy-y).^(2)/(sigmay).^2),'all');
    
    %plot
    figure(1)
    hold on
    plot(x,y,'b', 'LineWidth',1.2);
end %for

%Testing with Matlab built-in function
for N=1:3
f=polyfit(x,ynoisy,N);
F=polyval(f,x);

%error and residual
Merror(:,N)=F-ynoisy;
Mresidual(N)=sum(Merror(:,N),'all');

%Plots
figure(1)
hold on
plot(x,F,'r--','LineWidth',1.2)
end
legend('Data','Linear fit','Quadratic fit','Cubic fit',...
    'MATLAB built-in Linear fit', 'MATLAB built-in Quadratic fit',...
    'MATLAB built-in Cubic fit','Location','northwest');

figure(2)
plot(x,error(:,3));
xlabel('x');
ylabel('Error');
title('Error in Cubic Least Square fit')

%Displaying Results
disp('Problem 1-b)');
disp(['                Residual: ', num2str(residual)]);
disp(['MATLAB built-in Residual: ', num2str(Mresidual)]);

fprintf('\nProblem 1-c)\n');
disp('Reduced Chi^2');
disp(['Linear fit   : ', num2str(Chi_sq(1))]);
disp(['Quadratic fit: ', num2str(Chi_sq(2))]);
disp(['Cubic fit    : ', num2str(Chi_sq(3))]);

fprintf('\nProblem 1-d)\n');
fprintf('Since the cubic fit has the closest reduced Chi^2 value to 1 and because its error \n seem to mostly depict noise about zero, it appears to be the best fit')
%% Problem 2
load('test_interp.mat')

%2-a)
x=[-0.83,-0.2,0.4, 0.9, 1.22];  %X values
i=InterpIndex1D(xg,x);          %Indices i for x_i<x<x_i+1

%2-b)
y=[-0.92,-0.3,0.2, 0.7, 1.12];  %X values
[i2D,j2D]=InterpIndex2D(xg,yg,x,y); %Indices i for x_i<x<x_i+1 ...
                                    %    and j for y_j<y<x_j+1

%2-c)
f=f2D(:);
finterpmanualtest=BlInterp(xg,yg,f,x,y); %Manual Interpolation

%2-d)
[X,Y]=meshgrid(xg,yg);

%Manual Interpolation
finterpmanual=BlInterp(xg,yg,f,xgi,ygi);

%Matlab built-in Interpolation
finterp=interp2(X,Y,f2D,xgi,ygi);

[XGi, YGi]=meshgrid(xgi,ygi);
figure(3)
imagesc(finterpmanual);
