%% EXERCISE 1: Fundamentals of parameter estimation - part I %%%%%%%%%
%%%%%MAP estimation, MMSE estimation; MMAE estimation; ML estimation%%%%%%
% Student: Duong Xuan Luan - Student ID: S2236117
% Date: 19 Nov 2020

clear ; close all %start the new code

% SUMMARY 
% The measurement system: an ultrasonic depth gauge
% z= c/2*TOF; c is speed of sound, TOF is time flight 
% Errors:   
% - Second echo: ToF=4*x/c
% - Electrical noise: p(z|x)= P0*p(N(sigma,x))+P1*p(N(sigma,2x))
%                     P0+P1=1;
% Prior knowledge: 
% - p(x)= beta/(2*alpha*gamma(1/beta))*exp(-(|x-mean|/alpha)^beta)
%   mean=xmax+xmin/2; alpha=1/2(xmax-xmin)

%Input given values
global P0 sigma xmin xmax beta;
P0=0.95;
sigma=0.1;
xmin=1;
xmax=3;
beta=20;

%Interested range of x and z
xrange=0:0.05:7; %step width 0.05 of x 
zrange=0:0.05:7; %step width for 0.05 of z 

%% Question 1
%create a graph of px(x) and graphs of pz_x(z|x) against z 

%Considering values of x
x=[1.5,2];

%Calculate function px
px_range=zeros(length(xrange),1); %allocate array to store px 
for it=1:length(xrange)
    px_range(it)=px(xrange(it)); % store data
end

%Calculate function pz_x
pz_x1_range=zeros(length(zrange),1); %allocate array to store pz_x with x1
pz_x2_range=zeros(length(zrange),1); %allocate array to store pz_x with x2
for it=1:length(zrange)
    pz_x1_range(it)=pz_x(x(1),zrange(it)); % store data of pz_x with x1
end
for it=1:length(zrange)
    pz_x2_range(it)=pz_x(x(2),zrange(it)); % store data of pz_x with x1
end

%Plot graph
figure(1) % Graph of px according to x
plot(xrange,px_range);
xlabel('x');ylabel('p(x)');
title('Graph of prior probability p(x) according to x')
print(gcf,'Ass1-fig1.png','-dpng','-r500'); 
figure(2) % Graph of px according to x=1.5 and x=2
hold on
plot(zrange,pz_x1_range);
plot(zrange,pz_x2_range);
legend('x=1.5','x=2');
xlabel('z');ylabel('p(z|x)');
title('Graph of likelihood probability p(z|x) according to x=1.5 and x=2')
hold off
print(gcf,'Ass1-fig2.png','-dpng','-r500'); 


%% Question 2
%create a graph of pz(z) and graphs of px_z(x|z) against x 
%Considering value of z
z=[3.1,4];

%Calculate function pz
pz_range=zeros(length(zrange),1); %allocate array to store pz according to z1, z2
for it=1:length(zrange)
    pz_range(it)=pz(zrange(it),xrange); % store data
end

%Calculate function px_z
px_z1_range=zeros(length(xrange),1); %allocate array to store px_z with z1
px_z2_range=zeros(length(xrange),1); %allocate array to store px_z with z2
for it=1:length(xrange)
    px_z1_range(it)=px_z(xrange(it),z(1),xrange); % store data of px_z with z1
end
for it=1:length(xrange)
    px_z2_range(it)=px_z(xrange(it),z(2),xrange); % store data of px_z with x1
end

%Plot graph
figure(3) % Graph of pz according to z
plot(zrange,pz_range);
xlabel('z');ylabel('p(z)');
title('Graph of measured probability p(z) according to z')
print(gcf,'Ass1-fig3.png','-dpng','-r500'); 
figure(4) % Graph of p(x|z) according to z=3.1 and z=4
hold on
plot(xrange,px_z1_range);
plot(xrange,px_z2_range);
title('Graph of probability p(x|z) according to z=3.1 and z=4')
legend('z=3.1','x=4');
xlabel('x');ylabel('p(x|z)');
hold off
print(gcf,'Ass1-fig4.png','-dpng','-r500'); 
%% Question 3
% MMSE,MAP,MMAE, ML estimator x for z=3.1 and z=4 
% Change range for better acurate estimator 
xrange=0:0.01:5; %step width 0.01 of x 
zrange=0:0.01:5; %step width for 0.01 of z 

%allocate data
x_map=[0 0]; x_mmse=[0 0]; x_mmae=[0 0]; x_ml= [0 0];

x_mmse(1)=MMSE(z(1),xrange); 
x_mmse(2)=MMSE(z(2),xrange);
x_map(1)=MAP(z(1),xrange); 
x_map(2)=MAP(z(2),xrange);
x_mmae(1) = MMAE(z(1),xrange);
x_mmae(2) = MMAE(z(2),xrange);
x_ml(1)=ML(z(1),xrange);
x_ml(2)=ML(z(2),xrange);

fprintf('Solution for question 3');
x_mmse 
x_mmae
x_map 
x_ml


%% Question 4
%Calculate risks for each estimator and each z from Question 3
%Conditional Risk= int(C(x_est|x)*p(x|z)*dx)

%Cost function
C_quad=@(x_est,x) (x_est-x)^2; %quadratic cost functions
C_abs=@(x_est,x) abs(x_est-x); %absolute cost functions
C_uni=@(x_est,x) 1*(abs(x_est-x)>0.05)+0*(abs(x_est-x)<=0.05);%uniform cost functions

%allocate data
R_mmse_quad=[0 0];R_mmse_uni=[0 0];R_mmse_abs=[0 0];
R_mmae_quad=[0 0];R_mmae_uni=[0 0];R_mmae_abs=[0 0];
R_ml_quad=[0 0];R_ml_uni=[0 0];R_ml_abs=[0 0];
R_map_quad=[0 0];R_map_uni=[0 0];R_map_abs=[0 0];

%for z1=3.1
%MMSE
R_mmse_quad(1)=Risk(C_quad,xrange,x_mmse(1),z(1)); 
R_mmse_abs(1)=Risk(C_abs,xrange,x_mmse(1),z(1)); 
R_mmse_uni(1)=Risk(C_uni,xrange,x_mmse(1),z(1)); 
%MAP
R_map_quad(1)=Risk(C_quad,xrange,x_map(1),z(1)); 
R_map_abs(1)=Risk(C_abs,xrange,x_map(1),z(1)); 
R_map_uni(1)=Risk(C_uni,xrange,x_map(1),z(1)); 
%MMAE
R_mmae_quad(1)=Risk(C_quad,xrange,x_mmae(1),z(1)); 
R_mmae_abs(1)=Risk(C_abs,xrange,x_mmae(1),z(1)); 
R_mmae_uni(1)=Risk(C_uni,xrange,x_mmae(1),z(1)); 
%ML
R_ml_quad(1)=Risk(C_quad,xrange,x_ml(1),z(1)); 
R_ml_abs(1)=Risk(C_abs,xrange,x_ml(1),z(1)); 
R_ml_uni(1)=Risk(C_uni,xrange,x_ml(1),z(1)); 

%for z2=4
%MMSE
R_mmse_quad(2)=Risk(C_quad,xrange,x_mmse(2),z(2)); 
R_mmse_abs(2)=Risk(C_abs,xrange,x_mmse(2),z(2)); 
R_mmse_uni(2)=Risk(C_uni,xrange,x_mmse(2),z(2)); 
%MAP
R_map_quad(2)=Risk(C_quad,xrange,x_map(2),z(2)); 
R_map_abs(2)=Risk(C_abs,xrange,x_map(2),z(2)); 
R_map_uni(2)=Risk(C_uni,xrange,x_map(2),z(2)); 
%MMAE
R_mmae_quad(2)=Risk(C_quad,xrange,x_mmae(2),z(2)); 
R_mmae_abs(2)=Risk(C_abs,xrange,x_mmae(2),z(2)); 
R_mmae_uni(2)=Risk(C_uni,xrange,x_mmae(2),z(2)); 
%ML
R_ml_quad(2)=Risk(C_quad,xrange,x_ml(2),z(2)); 
R_ml_abs(2)=Risk(C_abs,xrange,x_ml(2),z(2)); 
R_ml_uni(2)=Risk(C_uni,xrange,x_ml(2),z(2));

fprintf('Solution for question 4- quad-z1/quad-z2/uni-z1/uni-z2/abs-z1/abs-z2');
R_mmse=[R_mmse_quad,R_mmse_uni,R_mmse_abs]
R_mmae=[R_mmae_quad,R_mmae_uni,R_mmae_abs]
R_map=[R_map_quad,R_map_uni,R_map_abs]
R_ml=[R_ml_quad,R_ml_uni,R_ml_abs]

