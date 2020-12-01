%% Introduction
%Aldous George
%EP 501
%Project 5
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/
clc
clearvars
close all
%% Problem 1

%Initiation
I   = 10;               %(A)
mu0 = 4*pi*10^(-7);     %(H/m)
a   = 0.004;            %(m)
x = linspace(-3*a,3*a,100);
y = linspace(-3*a,3*a,100);
[X,Y] = meshgrid(x,y);

Bx = zeros(100,100);
By = zeros(100,100);

for i=1:100
    for j=1:100
        cond = sqrt(x(i)^2+y(j)^2);
        if cond<a
            Bx(i,j) = (mu0*I/(2*pi*a^2))*cond*(-(y(j)/cond));
            By(i,j) = (mu0*I/(2*pi*a^2))*cond*(x(i)/cond);
        else
            Bx(i,j) = (mu0*I/(2*pi))*(1/cond)*(-(y(j)/cond));
            By(i,j) = (mu0*I/(2*pi))*(1/cond)*(x(i)/cond);
        end
    end
end

%1-a
%plot
%Bx
figure(1);
pcolor(X,Y,Bx);
xlim([-3*a,3*a]);
ylim([-3*a,3*a]);
xlabel('x');
ylabel('y');
title('B_x');
colorbar 
shading flat
hold on 
%By
figure(2)
pcolor(X,Y,By);
xlim([-3*a,3*a]);
ylim([-3*a,3*a]);
xlabel('x');
ylabel('y');
title('B_y');
colorbar 
shading flat

%1-b
figure(3)
quiver(X,Y,Bx',By');
title('B Quiver plot');
