%% EXERCISE 6:  extended Kalman filtering
% Student: Luan Duong - Student ID:S2236117
% Date: 11th Dec 2020
% Objectives of the exercise: 
%   - how to implement extended kalman filter: linearize system during
%   tracking

clear all; close all % start 

% Summary of question: System equation of yacht
%   - Simplified model:m*v'+d(v)v=t 
%       v is linear& angular velocity
%       d(v) is damping factor: d(v)=a*|v|+b*|v|^3   
%       external force t= [phi t]: dt=-alphat(t-t0)+noise
%   - Time discrete model:
%       X=[Xsi; v; acc; t; phi]; X(i+1)=f(x(i),u(i))+w(i); u=[t0; phi0]
%       Xsi(i+1)=Xsi(i)+delta_v(i)+w_Xsi(i)
%       v(i+1)=v(i)+delta_a(i)+w_v(i)
%       a(i+1)=1/m*(t(i)*[cos(phi(i)*pi/180);sin(phi(i)*pi/180)]-d(v(i))*v(i))+w_a(i)
%       t(i+1)=t(i)-delta_alphat(t(i)-t0)+w_t(i)
%       phi(i+1)=phi(i)-delta_alphaphi(phi(i)-phi0)+w_phi(i)
%       Given Case: 
        sigma_wXsi=10^(-2); %standard deviation process noise on position
        sigma_wv=10^(-2); %standard deviation process noise on velocity
        sigma_wa=10^(-2); %standard deviation process noise on acceleration
        t0=400; % intended thrust
        sigma_wt=8; %standard deviation process noise on thrust
        phi0=45; %intended heading
        sigma_wphi=0.5; %standard deviation process noise on heading
        alphat=1/1200; %coefficient differential equation for thrust
        alphaphi=1/200; %coefficient differential equation for heading
        delta=1; % sampling period
        m=3600; % mass
        a=30; % skin friction constant
        b=5; % wave making constant
%   - Measurement model: sampling period 500s, im=0,500,1000...
%       z(im)= h(x(im))+n(im)
%            = [180/pi*arctan(y0-Xsi_y(im)/x0-Xsi_x(im)); |v(im)|;180/pi*arctan(v_y(im)/v_x(im))]
%              + [n1(im);n2(im);n3(im)]
%       Given case: 
        sigma_n1=1; sigma_n2=0.3; sigma_n3=1;
        x0=[5000;10000]; % beacon position
%   - Prior knowledge
        Xsi=[0;0]; sigma_Xsix=10000;sigma_Xsiy=10000;
        v=[0;0]; sigma_vx=2;sigma_vy=2;
        acc=[0;0]; sigma_ax=0.04;sigma_ay=0.04;
        t=400; sigma_t=300;
        phi=0; sigma_phi=10;





%% Question 1: 
% PURPOSE: - Determine Cw for process noise and Cn for measurement noise
% 1. Determine Cw
% As the noise errors are uncorrelated
Cw=diag([sigma_wXsi^2,sigma_wXsi^2,sigma_wv^2,sigma_wv^2,...
    sigma_wa^2,sigma_wa^2, sigma_wt^2, sigma_wphi^2]); 


% 2. Determine Cn
% As the noise errors are uncorrelated
Cn=diag([sigma_n1^2,sigma_n2^2,sigma_n3^2]);

%% Question 2: 
% PURPOSE: load measurement z, apply EKF from i=0 to i=9999
load('z_yacht.mat') % load measurement 

% Prior knowledge x(0|-1) and C(0|-1)
x=[Xsi;v;acc;t;phi];
C=diag([sigma_Xsix^2,sigma_Xsiy^2,sigma_vx^2,sigma_vy^2,...
    sigma_ax^2,sigma_ay^2, sigma_t^2, sigma_phi^2]);

% Saving ininital data
Clist{1}=C;
xlist{1}=x;
e=[]; % Question 5

% EKF loop 
for it=1:delta:10000 
    % update process
    if mod(it-1,500)==0 % update only at im=0,500,1000....
        % Calculate predicted measurement
        z_pred=hmeas(x,x0);
        % Calculate innovation variation 
        S=Hjacobian(x,x0)*C*Hjacobian(x,x0).'+Cn;
        % The Kalman gain matrix: K=Cxz*Cz^-1=(Cx*H^T)*(H*Cx*H^T+Cn)^-1
        K=C*Hjacobian(x,x0).'/S;
        % Update x and C
        x=x+K*(z(:,((it-1)/500)+1)-z_pred);
        C=C-K*S*K.';
        % Question 5: Calculate NIS
        z_inn=z(:,((it-1)/500)+1)-z_pred;%innovation
        e(end+1)=z_inn.'/S*z_inn;
    end
    
    % prediction process
    % propagation of C
    C=Fjacobian(x)*C*Fjacobian(x).'+Cw;
    % prediction x=f(x,u,Cw)
    u=[t0;phi0];
    x=fsys(x,u); 
    
    % Record data
    xlist{it+1}=x;
    Clist{it+1}=C;
end 

%% Question 3: 
% PURPOSE: plot position, velocity, acceleration 
% plot position
figure(1) 
for i=1:length(xlist)
    Xsi(:,i)=xlist{i}(1:2);
end
subplot(2,2,[1 3])
plot(Xsi(1,:),Xsi(2,:),'r-');
xlabel('Xsi x'); ylabel('Xsi y');
subplot(2,2,2)
plot(1:length(xlist),Xsi(1,:),'r-');
xlabel('time'); ylabel('Xsi x');
subplot(2,2,4)
plot(1:length(xlist),Xsi(2,:),'r-');
xlabel('time'); ylabel('Xsi y');
subtitle('Estimate location')
print(gcf,'Ass6-fig1a.png','-dpng','-r500'); 

% plot velocity
figure(2) 
for i=1:length(xlist)
    v(:,i)=xlist{i}(3:4);
end
subplot(2,2,[1 3])
plot(v(1,:),v(2,:),'r-');
xlabel('v x'); ylabel('v y');
subplot(2,2,2)
plot(1:length(xlist),v(1,:),'r-');
xlabel('time'); ylabel('v x');
subplot(2,2,4)
plot(1:length(xlist),v(2,:),'r-');
xlabel('time'); ylabel('v y');
subtitle('Estimate velocity')
print(gcf,'Ass6-fig1b.png','-dpng','-r500'); 

% plot acceleration
figure(3) 
for i=1:length(xlist)
    acc(:,i)=xlist{i}(5:6);
end
subplot(2,2,[1 3])
plot(acc(1,:),acc(2,:),'r-');
xlabel('a x'); ylabel('a y');
subplot(2,2,2)
plot(1:length(xlist),acc(1,:),'r-');
xlabel('time'); ylabel('a x');
subplot(2,2,4)
plot(1:length(xlist),acc(2,:),'r-');
xlabel('time'); ylabel('a y');
subtitle('Estimate acceleration')
print(gcf,'Ass6-fig1c.png','-dpng','-r500'); 


% plot position
figure(4) 
for i=1:length(xlist)
    t_force(:,i)=xlist{i}(7:8);
end

plot(1:length(xlist),t_force(1,:),'r-');
xlabel('time'); ylabel('t');
subtitle('Estimate force')
print(gcf,'Ass6-fig1d.png','-dpng','-r500'); 


%% Question 4
% PURPOSE: Add the uncertainty region for position
figure(5)
hold on
p0=scatter(x0(1),x0(2),'g*');
p1=plot(Xsi(1,:),Xsi(2,:),'r-');
xlabel('\xi x'); ylabel('\xi y');
for i=1:500:10000
    p2=uncertainty(Clist{i}(1:2,1:2),Xsi(:,i));
end
legend([p0 p1 p2],'Beacon','Estimate location',...
    'Uncertainty', 'Location', 'northwest')
axis equal
title('The estimated location with uncertainty ')
hold off
print(gcf,'Ass6-fig2.png','-dpng','-r500'); 

%% Question 5 
% Use NIS to test the algorithm 
figure (6)
hold on
plot(1:length(e),e)
yline(7.81)
hold off
xlabel('Measurement number'); 
title('Normalized Inovation Squared')
print(gcf,'Ass6-fig4png','-dpng','-r500');
% 5% satisfied 
