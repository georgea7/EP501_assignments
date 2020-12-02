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

lx=100;                 %size of x
ly=100;                 %size of y
x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
[X,Y] = meshgrid(x,y);

Bx = zeros(100,100);
By = zeros(100,100);

for i=1:lx
    for j=1:ly
        cond = sqrt(x(i)^2+y(j)^2); %condition for piecewise function
        if cond<a
            Bx(i,j) = (mu0*I/(2*pi*a^2))*cond*(-(y(j)/cond));
            By(i,j) = (mu0*I/(2*pi*a^2))*cond*(x(i)/cond);
        else
            Bx(i,j) = (mu0*I/(2*pi))*(1/cond)*(-(y(j)/cond));
            By(i,j) = (mu0*I/(2*pi))*(1/cond)*(x(i)/cond);
        end %if
    end %for
end %for

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

%1-c
dBx=Bx(2)-Bx(1);
dBy=By(2)-By(1);
gradBx=zeros(size(Bx));
gradBy=zeros(size(By));

gradBx(1,1)=(Bx(2,2)-Bx(1,1))/dBx;
for i=2:lx-1
    for j=2:lx-1
    gradBx(j,i)=(Bx(j,i+1)-Bx(j-1,i-1))/2/dBx;    %\partial/\partial x
    end %for
end %for
gradBx(ly,lx)=(Bx(ly,lx)-Bx(ly-1,lx-1))/dBx;

%% Problem 2

%Initiation
Q = 1;              %(C)
a = 1;              %(m)
Eps0 = 8.854*10^-12;%(F/m)

for i=1:lx
    for j=1:ly
        for k=1:lz
            cond = sqrt(x^2+y^2+x^2); %condition for piecewise function
            if cond<a
            phix
            else
            phix
            end %if
        end%for
    end%for
end%for
