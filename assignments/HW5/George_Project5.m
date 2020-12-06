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

disp('Problem 1');
%Initiation
I   = 10;               %(A)
mu0 = 4*pi*10^(-7);     %(H/m)
a   = 0.005;            %(m)

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
disp('1-a');
%plot
%Bx
figure(1);
subplot(1,2,1);
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
figure(1)
subplot(1,2,2)
pcolor(X,Y,By);
xlim([-3*a,3*a]);
ylim([-3*a,3*a]);
xlabel('x');
ylabel('y');
title('B_y');
colorbar
shading flat

%1-b
disp('1-b')
figure(2)
quiver(X,Y,Bx',By');
title('B Quiver plot');

%1-c
disp('1-c')
dx=x(2)-x(1);
dy=y(2)-y(1);
gradBx=zeros(size(Bx));
gradBy=zeros(size(By));

Bx=Bx';
for i=1:lx      %Forward Difference
    gradBx(1,i)=(Bx(2,i)-Bx(1,i))/dy;                   %dBx/dy
    gradBy(1,i)=(By(2,i)-By(1,i))/dx;                   %dBy/dx
end


for j=1:ly      %Centered Difference
    for i=2:lx-1
        gradBx(i,j)=(Bx(i+1,j)-Bx(i-1,j))/2/dy;         %dBx/dy
        gradBy(i,j)=(By(i+1,j)-By(i-1,j))/2/dx;         %dBy/dx
    end %for
end %for

for i=1:lx      %Backward difference
    gradBx(lx,i)=(Bx(lx,i)-Bx(lx-1,i))/dy;              %dBx/dy
    gradBy(lx,i)=(By(lx,i)-By(lx-1,i))/dx;              %dBy/dx
end

curlB = gradBy-gradBx;
figure(3)
subplot(1,2,1)
pcolor(X,Y,curlB);
xlabel('x');
ylabel('y');
title('\nabla x B');
colorbar
shading flat

curlM = curl(X,Y,Bx,By');
figure(3)
subplot(1,2,2)
pcolor(X,Y,curlM);
xlabel('x');
ylabel('y');
title('MATLAB built-in \nabla x B');
colorbar
shading flat

%1-d
disp('1-d')
%% Problem 2
disp('Problem 2')
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

%2-a
disp('2-a');
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
title('\Phi');
colorbar
shading flat

%2-b
disp('2-b')
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

Laplace=divx+divy+divz;    %this is really laplacian b/c input is gradient

figure(6);
surface(X,Y,Laplace(:,:,end/2));
set(gca,'FontSize',15);
xlabel('x');
ylabel('y');
title('Laplacian(\Phi)');
colorbar;
shading flat

%% Problem 3
disp('Problem 3')
%3-a
disp('3-a')
int=0;
intg = phi.*Laplace;
for i=1:lx-1
    for j=1:ly-1
        for k=1:lz-1
            int= int+(1/8)*(intg(i,j,k)+intg(i+1,j,k)+intg(i,j+1,k)+...
                intg(i,j,k+1)+intg(i+1,j+1,k)+intg(i+1,j,k+1)+intg(i,j+1,k+1)+...
                intg(i+1,j+1,k+1))*Eps0*dx*dy*dz;
        end %for
    end %for
end %for
W=-0.5*int;
W=W/10^9;

fprintf('W = %.2f GJ\n\n',W);


%% Problem 4
disp('Problem 4')
%4-a
disp('4-a');
r0=2*a;
r=r0*(cos(phi)+sin(phi));
