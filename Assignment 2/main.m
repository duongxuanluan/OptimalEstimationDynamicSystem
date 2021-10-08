%% Exercise 2: Fundamentals of parameter estimation - part II 
%%%%%linear MMSE and unbiased linear MMSE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Student: Duong Xuan Luan - Student ID: S2236117
% Date: 24 Nov 2020
clear; close all
% System model: the ultrasonic depth gauge discussed in exercise 1 
%First data set: design a linMMSE, unbiased linMMSE estimator 
%Second data set: evaluate these estimators - bias and variance



%% Question 1 
% create an xz-scatter diagram using data set 1
load depthgauge_data_set1.mat % load data set 1

figure(1) % Graph of x-z
scatter(x1,z1,'r', '*');
legend('x1 versus z1','Location','southeast');
xlabel('Values of x'); ylabel('Values of z');
title('True Values x against measurements z - Data Set 1')
print(gcf,'Ass2-fig1.png','-dpng','-r500'); 
hold on 

%% Question 2
%Linear MMSE: X_est= alpha*z where alpha = E(xz)/E(z^2) 
[x1_lmmse,alpha_1]= LMMSE(x1,z1);

plot(x1_lmmse,z1, 'b*-'); % Graph of x_lmmse estimated and z1
legend('x1 versus z1','LMMSE x estimated versus z1');
print(gcf,'Ass2-fig2.png','-dpng','-r500'); 

%% Question 3
%unbias Linear MMSE: X_est= alpha*z +beta 
%alpha = cov(x,z)/var(z); beta= E(x)-alpha*E(z)
[x1_ulmmse,alpha_2,beta]= uLMMSE(x1,z1);

plot(x1_ulmmse,z1, 'g*-'); % Graph of x_lmmse estimated and z1
legend('x1 versus z1','LMMSE x estimated versus z1', ...
    'uLMMSE x estimated versus z1');
print(gcf,'Ass2-fig3.png','-dpng','-r500');
hold off 

%% Question 4
% Apply LMMS estimator calculated to data set 2
% Calculate errors, bias and variance of errors

load depthgauge_data_set2.mat % load data set 2
x2_lmmse=alpha_1.*z2;
x2_lmmse_error= x2_lmmse-x2;
bias_lmmse= E(x2_lmmse_error);
variance_lmmse= var(x2_lmmse_error);

%% Question 5
% Apply uLMMS estimator calculated to data set 2
% Calculate errors, bias and variance of errors

x2_ulmmse=alpha_2.*z2+beta.*ones(size(z2));
x2_ulmmse_error= x2_ulmmse-x2;
bias_ulmmse= E(x2_ulmmse_error);
variance_ulmmse= var(x2_ulmmse_error);



