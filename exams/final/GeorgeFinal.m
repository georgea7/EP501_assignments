%% Introduction
%Aldous George
%EP 501
%Final Exam
%This code contains excerpts from codes provided by Dr. Zettergen.
%https://github.com/Zettergren-Courses/EP501_matlab/
clc
clearvars
close all
%% Problem 1

n = 100;
dt = linspace(0,1,n);
alpha = 4;
G=zeros(1,n);
for i=1:n
        G(i) = 1-(dt(i)*alpha)+(0.5*(dt(i)*alpha)^2)-...
            ((1/6)*(dt(i)*alpha)^3)+((1/24)*(dt(i)*alpha)^4);
end

figure(1);
alphadt=alpha.*dt;
plot(alphadt,G,'o');
set(gca,'FontSize',15);
ylabel('G');
xlabel('\alpha \Delta t');
hold on
plot(0:0.04:4,ones(101),'r','LineWidth',2);
legend('Gain factor','Stability limit line','Location','NorthWest');
title('G vs \alpha dt (\alpha = 4)')

%1-c
x=linspace(0,20,1000);        %radial independent variable
maxit =100;              %Max Iterations 
tol=1e-3;                %Tolerance
verbose =false;
g = @objfun1;
gprime = @objgradfun1;
ygrid=g(x);
g2 = @objfun2;
ygrid2 = ygrid+2;

for j=1:2
    x0=(j-1);    %Initial condition x_i
    [R(j),it,success]=newton_exact(g,gprime,x0,maxit,tol,verbose);
end
disp(['Real Roots: ', num2str(R)]);
figure(2)
subplot(1,2,1)
plot(x,ygrid)
set(gca,'FontSize',15);
ylim([-1,10]);
xlabel('\Delta t');
ylabel('gain factor');
title('G(\alpha ,\Delta t)-1 = 0')
subplot(1,2,2)
plot(x,ygrid2)
set(gca,'FontSize',15);
ylim([-1,10]);
title('G(\alpha ,\Delta t)+1 = 0')
xlabel('\Delta t');
ylabel('gain factor');

%1-d 
%Gridding in time
N=16;
tmin=0;
tmax=10;
t=linspace(tmin,tmax,N);
%t=[tmin:R(2):tmax];
dt=t(2)-t(1);

%true solution
y0=1;
ybar=y0*exp(-alpha*t);

%RK2
yRK2=zeros(1,N);
yRK2(1)=y0;
for n=2:N
    yhalf=yRK2(n-1)+dt/2*(-alpha*yRK2(n-1));
    yRK2(n)=yRK2(n-1)+dt*(-alpha*yhalf);
end %for

%RK4
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

figure(3);
clf;
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',15);
hold on;
plot(t,yRK2,'--')
plot(t,yRK4,'^-')
ylim([0,8]);
legend('Exact','RK2','RK4')
title('f(t,y) = -\alpha y');

%1-e
%Gridding in time
%Recalculating below the limit
N=15;                           %Slightly below the limit
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);

ybar=y0*exp(-alpha*t);

%RK4
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

figure(4);
subplot(1,2,1);
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',15);
hold on
plot(t,yRK4,'^-')
legend('Exact','RK4')
title('\Delta t = 0.7143');

%Recalculating above the limit
N=16;                           %Slightly above the limit
t=linspace(tmin,tmax,N);
dt=t(2)-t(1);

ybar=y0*exp(-alpha*t);

%RK4
yRK4=zeros(1,N);
yRK4(1)=y0;
for n=2:N
    dy1=dt*fRK(t(n-1),yRK4(n-1),alpha);
    dy2=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy1/2,alpha);
    dy3=dt*fRK(t(n-1)+dt/2,yRK4(n-1)+dy2/2,alpha);
    dy4=dt*fRK(t(n-1)+dt,yRK4(n-1)+dy3,alpha);
    
    yRK4(n)=yRK4(n-1)+1/6*(dy1+2*dy2+2*dy3+dy4);
end %for

figure(4);
subplot(1,2,2);
plot(t,ybar,'o-');
xlabel('t');
ylabel('y(t)');
set(gca,'FontSize',15);
hold on
plot(t,yRK4,'^-')
% ylim([0,8]);
legend('Exact','RK4')
title('\Delta t = 0.6667');

%% Problem 2


