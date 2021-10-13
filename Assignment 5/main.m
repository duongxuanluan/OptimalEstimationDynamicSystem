%% EXERCISE 5:  discrete kalman filtering
% Student: Luan Duong - Student ID:S2236117
% Objectives of the exercise: 
%   - how to implement discrete kalman filter 

% Summary of question: Radar tracking
%   - Xsi(i):2D position, v(i):velocity,a(i):acceleration  
%   - Model: x(i+1)=F*x(i)+w(i), w(i) is white noise 
%   - State vector: x(i)= (Xsi(i),v(i),a(i))
%   - Xsi(i+1)=Xsi(i)+v(i), v(i+1)=v(i)+a(i),a(i+1)=F1*a(i)+w1(i)

clear all; close all % start 
% Given Case: 
F1=[0.97 0;0 0.97]; 
Cw1=[0.0016 0; 0 0.0016];

% Radar system: z(i)=Xsi(i)+n(i), n(i) is measurement noise 
Cn=[49 0; 0 49];

%% Question 1: 
% PURPOSE: - Determine F, Cw, Cx(0) and x_pred(0|-1)
% 1. Determine F
% Xsi(i+1)=Xsi(i)+v(i) => x_xsi(i+1)=x_xsi(i)+x_v(i)+ w
F_Xsi=[1 0 1 0 0 0; 0 1 0 1 0 0];
% v(i+1)=v(i)+a(i) => x_v(i+1)=x_v(i)+x_a(i)+ w
F_v=[0 0 1 0 1 0; 0 0 0 1 0 1];
% a(i+1)=F1*a(i)+w 
F_a=[0 0 0 0 F1(1,:); 0 0 0 0 F1(2,:)];
F=[F_Xsi;F_v;F_a]; % Matrix of F

% 2. Determine Cw
% Xsi(i+1)=Xsi(i)+v(i), v(i+1)=v(i)+a(i) => w_Xsi=w_v=0
Cw=[zeros(4,6);zeros(2,4),Cw1];

% 3. Determine Cx(0)=p(x(0))
% Given: uncertain range of s_Xi(0)=100, s_v(0)=4, s_a(0)= 0.2
% C_pred(0|-1)= Cx=E[(x-x_mean)*(x-x_mean)^T), the variables are incorrelation 
Cpred=diag([100^2 100^2 4^2 4^2 0.2^2 0.2^2]);


% 4. xpred(0|-1)=E(x(0))=0
xpred=[0;0;0;0;0;0];

%% Question 2: 
% PURPOSE: load measurement z, create DKF loop to find x_est
load('zradar.mat') % load measurement 

% Saving ininital data
% Cpredlist{1}=Cpred;
% xpredlist{1}=xpred;

% Define H
H=zeros(2,6);H(1,1)=1; H(2,2)=1;

% DKF loop 
for it=1:length(z) 
    % update process
        % The Kalman gain matrix: K=Cxz*Cz^-1=(Cx*H^T)*(H*Cx*H^T+Cn)^-1
        K=Cpred*H.'/(H*Cpred*H.'+Cn);
        % the uLMMSE estimate of the position: x_est=Kz+b
        b=(eye(length(K))-K*H)*xpred;
        xest=K*z(:,it)+b;
        % the corresponding error cov matrix:
        Cerr=Cpred-K*H*Cpred.';
        % Saving data
        xestlist{it}=xest;
        Cerrlist{it}=Cerr;
        K_list{it}=K;
    % predict process
        % prediction xpred(i+1|i)=F*xest(i|i)//+L*u(i)//
        xpred=F*xest;
        % propagation of covariance matrix
        Cpred=F*Cerr*F.'+Cw;
        % Saving data
        xpredlist{it}=xpred;
        Cpredlist{it}=Cpred;
    
    % Question 5: Comment the above part of the loop, uncomment below
%         [K,Cpred,Cerr]=dlqe(F,eye(6),H,Cw,Cn);
%         b=(eye(length(K))-K*H)*xpred;
%         xest=K*z(:,it)+b;
%         xestlist{it}=xest;
%         Cerrlist{it}=Cerr;
%         K_list{it}=K;
%         xpred=F*xest;
%         xpredlist{it}=xpred;
%         Cpredlist{it}=Cpred;
        
    
end 

%% Question 3: 
% PURPOSE: plot Xipred,Xiest and z
figure(1) 
hold on
p1=plot(z(1,:),z(2,:),'r-.');
for i=1:length(xestlist)
    Xiest(:,i)=xestlist{i}(1:2);
end
p2=plot(Xiest(1,:),Xiest(2,:),'g-.');
for i=1:length(xpredlist)
    Xipred(:,i)=xpredlist{i}(1:2);
end
p3=plot(Xipred(1,:),Xipred(2,:),'b-.');
axis equal
xlabel('z0'), ylabel('z1')
legend([p1 p2 p3],'measurement', 'estimation', 'prediction',...
    'Location', 'southeast')
title('Measurement vs estimation vs prediction location')
%print(gcf,'Ass5-fig1.png','-dpng','-r500'); 
%print(gcf,'Ass5-fig3.png','-dpng','-r500');

%% Question 4
% PURPOSE: Add the uncertainty region 
for i=1:3:100
    p4=uncertainty(Cerrlist{i}(1:2,1:2),Xiest(:,i));
end
legend([p1 p2 p3 p4],'measurement', 'estimation', 'prediction', 'uncertainty')
title('Measurement vs estimation vs prediction location with uncertainty')
%print(gcf,'Ass5-fig2.png','-dpng','-r500'); 
%print(gcf,'Ass5-fig4.png','-dpng','-r500'); 

%% Question 5
% Using the dlqe function 
K_list_before=K_list;
Cerrlist_before=Cerrlist;
Cpredlist_before=Cpredlist;

% Delete clear all at the beginning, Activate part of comment 2
% K_list_after=K_list;
% Cerrlist_after=Cerrlist;
% Cpredlist_after=Cpredlist;

% for i=1:length(K_list_before)
% K_diff(i)=norm(K_list_after{i}-K_list_before{i});
% Cerr_diff(i)=norm(Cerrlist_after{i}-Cerrlist_before{i});
% Cpred_diff(i)=norm(Cpredlist_after{i}-Cpredlist_before{i});
% end
% figure (2)
% hold on
% p5=plot(1:length(K_list_before),K_diff(:));
% p6=plot(1:length(K_list_before),Cerr_diff(:));
% p7=plot(1:length(K_list_before),Cpred_diff(:));
% legend([p5 p6 p7],' Difference of K', 'Different of Cerr', 'Different of Cpred')
% xlabel('Number of measurement')
% title('Measurement vs estimation vs prediction location with uncertainty')
% hold off
% print(gcf,'Ass5-fig5.png','-dpng','-r500'); 

