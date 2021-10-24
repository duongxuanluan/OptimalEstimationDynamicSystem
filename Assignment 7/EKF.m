%% EXERCISE 7: particle filtering
% Student: Luan Duong - Student ID:S2236117
% Objectives of the exercise: 
%   - how to implement particle filter
%   - perform consistency check, compraing result with EKF 

clear all; close all % start 

% Summary of question: System equation of yacht in Ex 6:
%   - Simplified model:m*v'+d(v)v=t 
%       v is linear& angular velocity
%       d(v) is damping factor: d(v)=a*|v|+b*|v|^3   
%       external force t= [phi t]: dt=-alphat(t-t0)+noise
%   - Time discrete model:
%       X=[Xsi; v; acc; t; phi]; X(i+1)=f(x(i),u(i)+w(i); u=[t0; phi0]
%       Xsi(i+1)=Xsi(i)+delta_v(i)+w_Xsi(i)
%       v(i+1)=v(i)+delta_a(i)+w_v(i)
%       a(i+1)=1/m*(t(i)*[cos(phi(i)*pi/180);sin(phi(i)*pi/180)]-d(v(i))*v(i))+w_a(i)
%       t(i+1)=t(i)-delta_alphat(t(i)-t0(i))+w_t(i)
%       phi(i+1)=phi(i)-delta_alphaphi(phi(i)-phi0(i))+w_phi(i)
%       * Differences from Ex6: t0,phi0 depending on time
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
        Xsi=[0;0]; sigma_Xsi=[100;100]; % Changing sigma from Ex6
        v=[0;0]; sigma_v=[2;2];
        acc=[0;0]; sigma_a=[0.04;0.04];
        t=400; sigma_t=300;
        phi=0; sigma_phi=10;
%   - Consitency checks for particle filtering 
%       u_n(i)=F_n(z_n(i)|i) where F is culmulative distribution of z_n(i)
%       u_n uniformly distributed ?
%       Keff(i)= 1/(sum(w_norm(i)^2) => evenly distributed ?

load('x_yacht2.mat'); % Load phi0(t)
load('z_yacht2.mat'); % Load the dataset (measurements in Z)
Z=z;
I = length(x); % Length of the sequence
Cw=diag([sigma_wXsi^2,sigma_wXsi^2,sigma_wv^2,sigma_wv^2,...
    sigma_wa^2,sigma_wa^2, sigma_wt^2, sigma_wphi^2]);
Cn=diag([sigma_n1^2,sigma_n2^2,sigma_n3^2]);

%% Question 6 
% PURPOSE: - Perform EKF and compare.
% EKF
% Prior knowledge x(0|-1) and C(0|-1)
x_EKF=[Xsi;v;acc;t;phi];
C_EKF=diag([sigma_Xsi(1)^2,sigma_Xsi(2)^2,sigma_v(1)^2,sigma_v(2)^2,...
    sigma_a(1)^2,sigma_a(2)^2, sigma_t^2, sigma_phi^2]);

% Saving ininital data
Clist_EKF{1}=C_EKF;
xlist_EKF{1}=x_EKF;
Xi_EKF(:,1)=x_EKF(1:2);
e=[]; 

% EKF loop 
for it=1:I
    % update process
    if mod(it-1,500)==0 % update only at im=0,500,1000....
        % Calculate predicted measurement
        z_pred=hmeas(x_EKF,x0);
        % Calculate innovation variation 
        S=Hjacobian(x_EKF,x0)*C_EKF*Hjacobian(x_EKF,x0).'+Cn;
        % The Kalman gain matrix: K=Cxz*Cz^-1=(Cx*H^T)*(H*Cx*H^T+Cn)^-1
        K=C_EKF*Hjacobian(x_EKF,x0).'/S;
        % Update x and C
        a=z(:,((it-1)/500)+1)-z_pred;
        a(1,:) = wrapTo180(a(1,:));
        a(3,:) = wrapTo180(a(3,:));
        x_EKF=x_EKF+K*a;
        C_EKF=C_EKF-K*S*K.';
        % Question 5: Calculate NIS
        e(end+1)=a.'/S*a;
        
    end
    
    % prediction process
    % propagation of C
    C_EKF=Fjacobian(x_EKF)*C_EKF*Fjacobian(x_EKF).'+Cw;
    Xi_EKF(:,it)=x_EKF(1:2);
    % prediction x=f(x,u,Cw)
    u=[t0;fi0(it)];
    x_EKF=fsys(x_EKF,u); 
    
    % Record data
    xlist_EKF{it+1}=x_EKF;
    Clist_EKF{it+1}=C_EKF;
    
end 
figure(6)
hold on
p0=scatter(x0(1),x0(2),'g*');
p1=plot(Xi_EKF(1,:),Xi_EKF(2,:),'r-');
xlabel('\xi x'); ylabel('\xi y');
for i=1:1000:I
    p2=uncertainty(Clist_EKF{i}(1:2,1:2),Xi_EKF(:,i));
end
legend([p0 p1 p2],'Beacon','Estimate location',...
    'Uncertainty', 'Location', 'northwest')
axis equal
title('The EKF estimated location with uncertainty ')
hold off
print(gcf,'Ass7-fig6.png','-dpng','-r500'); 


figure (7)
hold on
plot(1:length(e),e)
yline(7.81)
hold off
xlabel('Measurement number'); 
title('Normalized Inovation Squared')
print(gcf,'Ass7-fig7.png','-dpng','-r500');
% 5% satisfied 
