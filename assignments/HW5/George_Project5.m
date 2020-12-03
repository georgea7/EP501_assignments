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

Bx = zeros(lx,ly);
By = zeros(lx,ly);

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

gradBx(1,1)=(Bx(2,1)-Bx(1,1))/dBx;
for i=1:lx
    for j=2:ly-1
        gradBx(i,j)=(Bx(i,j+1)-Bx(i,j-1))/2/dBy;    %\partial/\partial x
    end %for
end %for
gradBx(lx,ly)=(Bx(lx,ly)-Bx(lx,ly-1))/dBy;

gradBy(1,1)=(By(1,2)-By(1,1))/dBx;
for i=2:lx-1
    for j=1:ly
        gradBy(i,j)=(Bx(i+1,j)-Bx(i-1,j))/2/dBx;    %\partial/\partial x
    end %for
end %for
gradBy(lx,ly)=(By(lx-1,ly)-By(lx-1,ly))/dBy;

curlB = gradBy-gradBx;
figure(4)
pcolor(X,Y,curlB);
xlabel('x');
ylabel('y');
title('Curl of B');
colorbar
shading flat

%% Problem 2

%Initiation
Q = 1;              %(C)
a = 1;              %(m)
Eps0 = 8.854*10^-12;%(F/m)

lx = 100;                 %size of x
ly = 100;                 %size of y
lz = 100;                 %size of z
x = linspace(-3*a,3*a,lx);
y = linspace(-3*a,3*a,ly);
z = linspace(-3*a,3*a,lz);

phi=zeros(lx,ly,lz);
[X,Y] = meshgrid(x,y);

for i=1:lx
    for j=1:ly
        for k=1:lz
            cond = sqrt(x(i)^2+y(j)^2+z(k)^2); %condition for piecewise function
            if cond<a
                phi(i,j,k) = Q/(4*pi*Eps0*a)-(Q/(8*pi*Eps0*a^3))*(x(i)^2 + y(j)^2 + z(k)^2 -a^2);
            else
                phi(i,j,k) = Q/(4*pi*Eps0*cond);
            end %if
        end%for
    end%for
end%for
figure(5)
pcolor(X,Y,phi(:,:,end/2));
xlim([-3*a,3*a]);
ylim([-3*a,3*a]);
set(gca,'FontSize',15);
xlabel('x');
ylabel('y');
title('Phi');
colorbar
shading flat

%Laplacian
dx = x(2)-x(1);
dy = y(2)-y(1);
dz = z(2)-z(1);

f=phi;
g=phi;
h=phi;

for l=1:2
    %x-derivative part of the divergence
    divx=zeros(size(f));
    divx(1,:,:)=(f(2,:,:)-f(1,:,:))/dx;
    for i=2:lx-1
        divx(i,:,:)=(f(i+1,:,:)-f(i-1,:,:))/2/dx;
    end %for
    divx(lx,:,:)=(f(lx,:,:)-f(lx-1,:,:))/dx;
    
    %y-derivative part of the divergence
    divy=zeros(size(g));
    divy(:,1,:)=(g(:,2,:)-g(:,1,:))/dy;
    for j=2:ly-1
        divy(:,j,:)=(g(:,j+1,:)-g(:,j-1,:))/2/dy;
    end %for
    divy(:,ly,:)=(g(:,ly,:)-g(:,ly-1,:))/dy;
    
    %z-derivative part of the divergence
    divz=zeros(size(h));
    divz(:,:,1)=(h(:,:,2)-h(:,:,1))/dz;
    for k=2:lz-1
        divz(:,:,k)=(h(:,:,k+1)-h(:,:,k-1))/2/dz;
    end %for
    divz(:,:,lz)=(h(:,:,lz)-h(:,:,lz-1))/dz;
    f=divx;
    g=divy;
    h=divz;
end

div=divx+divy+divz;    %this is really laplacian b/c input is gradient

figure(6);
surface(X,Y,div(:,:,end/2));
set(gca,'FontSize',15);
xlabel('x');
ylabel('y');
title('Laplacian(phi)');
colorbar;
shading flat

%% Problem 3

