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

%% Question 1: 
% PURPOSE: - implement a particle filter. use MMSE, create the estimated
%            path

load('x_yacht2.mat'); % Load phi0(t)
load('z_yacht2.mat'); % Load the dataset (measurements in Z)
Z=z;
I = length(x); % Length of the sequence
M = 8; % Dimension of state vector
N = 3; % Dimension of measurement vector

% Initialisation
Cw=diag([sigma_wXsi^2,sigma_wXsi^2,sigma_wv^2,sigma_wv^2,...
    sigma_wa^2,sigma_wa^2, sigma_wt^2, sigma_wphi^2]);
Cn=diag([sigma_n1^2,sigma_n2^2,sigma_n3^2]);
invCn = inv(Cn); % Inv. cov. of sensor noise
Ncond = 500; % Number of particles

%Generates the samples: Gaussian distribution with mean and std dev
Xs(:,:) = [Xsi;v;acc;t;phi]+[sigma_Xsi;sigma_v;sigma_a;sigma_t;sigma_phi].*randn(8,Ncond);
z_counter=1;
for i = 1:I  
    if mod(i-1,500)==0 % update only at im=0,500,1000....
        % Generate predicted meas. representing p(z(i)|Z(i-1))
        % *NOTE: avoid loop, hmeas can handle all the calculation at one call
        Zs = hmeas(Xs,x0,Cn);
    
        % Question 3: Get uniform distributed rv
        u(1,z_counter)= sum((Zs(1,:) < Z(1,(i-1)/500+1)))/Ncond;
        u(2,z_counter)= sum((Zs(2,:) < Z(2,(i-1)/500+1)))/Ncond;
        u(3,z_counter)= sum((Zs(3,:) < Z(3,(i-1)/500+1)))/Ncond;
        
        % Update
        res = Zs-Z(:,(i-1)/500+1)*ones(1,Ncond); % Residuals
        res(1,:) = wrapTo180(res(1,:));
        res(3,:) = wrapTo180(res(3,:));
        
        W = exp(-0.5*sum(res.*(invCn*res)))'; % Weights
        if (sum(W) == 0), error('process did not converge'); end
        W = W/sum(W); CumW = cumsum(W);
        

        xest(:,z_counter) = Xs(:,:)*W; % Sample mean
        
        % Question 2: Covariance 
        Cx{z_counter}= (Xs(1:2,:) - xest(1:2,z_counter)).*W'*(Xs(1:2,:) ...
            - xest(1:2,z_counter))';
       
        
        % Question 3: Effective number of particles 
        Klist(z_counter)=1/(sum(W.^2));
        
        % Find an index permutation using golden rule root finding
        for j = 1:Ncond
            R = rand; ja = 1; jb = Ncond;
            while (ja < jb-1)
                jx = floor(jb-0.382*(jb-ja));
                fa = R-CumW(ja); fb = R-CumW(jb); fxx = R-CumW(jx);
                if (fb*fxx < 0), ja = jx; else, jb = jx; end
            end
            ind(j) = jb;
        end

        % Resample
        for j = 1:Ncond, Xs(:,j) = Xs(:,ind(j)); end
        % Predict
        ut=[t0;fi0(i)]; % phi0 gained from data 
        Xs = fsys(Xs,ut,Cw); % Update samples
        
        z_counter=z_counter+1;
        
    else % not update
    % Predict
    ut=[t0;fi0(i)]; % phi0 gained from data 
    Xs = fsys(Xs,ut,Cw); % Update samples
    end
end

figure(1)
t = 1:Ncond;
plot(t, CumW)
title('The cumulative value of the weights')
xlabel('Sample number')
ylabel('Cum value of W')
%print(gcf,'Ass7-fig1.png','-dpng','-r500'); 
%print(gcf,'Ass7-fig8c.png','-dpng','-r500');% Changing for different Ncond 

figure(2)
hold on
plot(x(1,:),x(2,:))
plot (xest(1,:),xest(2,:));
scatter(x0(1),x0(2))
xlabel('\xi x'); ylabel('\xi y');
title('Plot of the estimated path') 
legend('Real path','Estimated path', 'Beacon')
%print(gcf,'Ass7-fig2.png','-dpng','-r500'); 

%% Question 2: 
% PURPOSE: - Calculate C(i|i) for each measurement and draw uncertainty
%            area
% Solve: Look at question 1
figure(3)
hold on
p1=plot (xest(1,:),xest(2,:));
p3=scatter(x0(1),x0(2));
for m = 1:length(z(1,:))
    p2=uncertainty(Cx{m}(1:2,1:2),[xest(1,m);xest(2,m)]);
end
hold off
xlabel('Xsi x'); ylabel('Xsi y');
title('Plot of the estimated path') 
legend([p1 p2 p3], 'Estimated path', 'Uncertainty area','Beacon'...
    ,'Location', 'northwest')
%print(gcf,'Ass7-fig3.png','-dpng','-r500');
%print(gcf,'Ass7-fig11c.png','-dpng','-r500');

%% Question 3: 
% PURPOSE: - Calculate u for each measurement and draw the graph
% Solve: Look at question 1
figure (4)
hold on
plot(1:length(u),u(1,:),'r.')
plot(1:length(u),u(2,:),'b.')
plot(1:length(u),u(3,:),'g.')
hold off
legend('u1','u2','u3')
xlabel('time im/500'); ylabel('u');
title('Plot of the test variable u') 
%print(gcf,'Ass7-fig4.png','-dpng','-r500');
%print(gcf,'Ass7-fig9c.png','-dpng','-r500');

%% Question 4: 
% PURPOSE: - Calculate K_eff for each measurement and draw the graph
% Solve: Look at question 1
figure (5)
scatter(1:length(Klist),Klist,'*')
xlabel('Measurement'); ylabel('K');
title('Plot of the effective number of particles')
%print(gcf,'Ass7-fig5.png','-dpng','-r500');
print(gcf,'sahel 3.png','-dpng','-r500');

%% Question 5
% PURPOSE: - Perform experiments to make sure sufficient no of particles.
% print(gcf,'Ass7-fig10c.png','-dpng','-r500'); % Changing for different Ncond
%% Question 6 
% PURPOSE: - Perform EKF and compare.
% EKF file 
