%% EXERCISE 3:  Fundamentals of parameter estimation - part III
% Student: Luan Duong - Student ID:S2236117
% Date: 27 Nov 2020
% Objectives of the exercise:
%   - the concept of covariance matrices
%   - Unbiased linear MMSE estimation

% Summary of question: Ship location(x,y) navigation problem 
%   - Prior knowledge: prior expectation Ex and a covariance matrix Cx
%   - Measurement theta=phi+theta_err from ship to beacon (x0,y0): 
%     Line of sight (Xi,Eta) plane: x0sin(theta)- y0cos(theta)= Xisin(theta)- Etacos(theta)
%     Real location: x0sin(phi)+ y0cos(phi)= xsin(phi)+ ycos(phi)
%   - Model derived from above eqn x=(x,y): z(theta,x0)=H(theta)*x+v  
%        x0sin(theta)- y0cos(theta)= xsin(theta)- ycos(theta) + d*theta_err
%        d- distance between beacon and ship: d =|x0-Ex|
%        SD_v=d*SD_theta_err (unit: radian)
%   - Uncertainty regions and principal axes:
%        Contour ellipse area: (x-Ex)^T*Cx^-1*(x-Ex)=k^2 (k=1,2,3)
%        Principal axis of ellipse: v0=Cxeigvec0/a0,v1=Cxeigvec1/a1
%        a1,a2 = (Cxeigval)^(1/2)

clear all; close all 
% Given case
Ex=[10;20]; Cx=[ 25 -25; -25 70]; % prior knowledge
x0=[100;100]; theta=deg2rad(35);SD_theta_err=deg2rad(1); % measurement 
% Note- deg2rad is to convert degree to radian

%% Question 1 
% PURPOSE: Find eig(Cx) and the MATLAB code to draw the uncertainty region

%Activate the below comment for question 6
% for i=1:9
%    SD_theta_err=deg2rad(i);
%    figure(1)
%    sgtitle(' Question 6: Comparision when SD theta error is changing')
%    subplot(3,3,i)

%Activate the below comment for question 7
% for i=1:9
%    Cx=[ 25 -25; -25 70]*i;
%    figure(1)
%    sgtitle(' Question 7: Comparision when Cx is changing')
%    subplot(3,3,i)

hold on % continuous drawing in the second question 
p1=uncertainty(Cx,Ex,'b','r*'); %b for blue color, r* for red star 
% xlim([0 30]); %limit for the graph
% legend( 'Uncertainty', 'Mean');
% xlabel('\xi');ylabel('\eta')
% title('Uncertainty circle of covariance Cx')
% print(gcf,'Ass3-fig2.png','-dpng','-r500'); 

%% Question 2 
% PURPOSE: Add the line of sight and uncertainty range of measurement to
% the figure

% Add the line of sight: x0sin(theta)- y0cos(theta)= Xisin(theta)- Etacos(theta)
Xi=-50:1:150; % range of sample for the line 
Eta=(x0(1)*sin(theta)- x0(2)*cos(theta)-sin(theta)*Xi)./(-cos(theta));
p2=plot(Xi,Eta,'y',x0(1),x0(2),'rx');

% Add the uncertainty area of the measurement [theta-SD_theta_err,theta+SD_theta_err]
Eta_err1=(x0(1)*sin(theta-SD_theta_err)- x0(2)*cos(theta-SD_theta_err)-sin(theta-SD_theta_err)*Xi)./(-cos(theta-SD_theta_err));
Eta_err2=(x0(1)*sin(theta+SD_theta_err)- x0(2)*cos(theta+SD_theta_err)-sin(theta+SD_theta_err)*Xi)./(-cos(theta+SD_theta_err));
p3=plot(Xi,Eta_err1,'m--',Xi,Eta_err2,'m--'); % lines of the uncertainty range
% xlim([0 160]); %limit for the graph
% legend( 'Uncertainty range of position', 'Mean','Line of sight',... 
%     'Beacon', 'Uncertainty region','Location','southeast');
% title('Line of sight and uncertainty region due to bearing measurement')
% print(gcf,'Ass3-fig3.png','-dpng','-r500'); 

%% Question 3
% PURPOSE: Draw the linearized measurement function (as a bar) instead of
% the cone above
% Linearized measurement function: x0sin(theta)- y0cos(theta)= xsin(theta)- ycos(theta) + d*theta_err

% Calculate the width between uncertainty range
d=norm(x0-Ex);
SD_v=d*SD_theta_err;

% Add the uncertainty area of the measurement using linearized function
Eta_linear_err1=(x0(1)*sin(theta)- x0(2)*cos(theta)-sin(theta)*Xi+SD_v)./(-cos(theta));
Eta_linear_err2=(x0(1)*sin(theta)- x0(2)*cos(theta)-sin(theta)*Xi-SD_v)./(-cos(theta));
p4=plot(Xi,Eta_linear_err1,'c--',Xi,Eta_linear_err2,'c--'); % lines of the uncertainty range
% xlim([0 120]); %limit for the graph
% legend( [p1(1) p1(2) p2(1) p2(2) p3(1) p4(1)],...
%     'Uncertainty range of position', 'Mean',...
%     'Line of sight','Beacon', 'Uncertainty region cone',...
%     'Uncertainty region bar','Location','southeast');
% title('Uncertainty region using linearized measurement system')
% print(gcf,'Ass3-fig4.png','-dpng','-r500'); 

%% Question 4 
% PURPOSE: - Find z, H and the Kalman gain matrix
%          - the uLMMSE estimate of the position: x_est=Kz+b
%          - the corresponding error cov matrix.

% Find z, H 
z=x0(1)*sin(theta)- x0(2)*cos(theta); % Value of z
H=[sin(theta) -cos(theta)]; % Value of H
% The Kalman gain matrix: K=Cxz*Cz^-1=(Cx*H^T)*(H*Cx*H^T+Cv)^-1
Cv=SD_v.^2;
K=Cx*H.'/(H*Cx*H.'+Cv);
% the uLMMSE estimate of the position: x_est=Kz+b
b=(eye(length(K))-K*H)*Ex;
x_est=K*z+b; 
% the corresponding error cov matrix:
Ce=Cx-K*H*Cx.';


%% Question 5 
% PURPOSE: Draw the posterior mean and uncertainty region 
p5=uncertainty(Ce,x_est,'r','b*');
% Activate the following comment if legend, title and label is needed
% xlabel('Xi'); ylabel('Eta');
% legend([p1(1) p1(2) p2(1) p2(2) p3(1) p4(1) p5(1) p5(2)],...
%     'Prior uncertainty area', 'Prior mean of the ship',  ...
%     'Line of sight', 'Beacon', 'Uncertainty of measurement (cone)',...
%     'Uncertainty of measurement (bar)',...
%     'Posterior uncertainty area', 'Posterior mean of the ship', ...
%     'Location','southwest');
% title('Prior and posterior uncertainty of ship location and the measurement')
 ylim([0 50]);  %limit for the graph
% print(gcf,'Ass3-fig5.png','-dpng','-r500');
hold off    



%% Question 6
% PURPOSE: Draw the graph with different SD_theta_err
% Delete the comment in Question 1 and the comment below
%title(['\sigma_{\delta\theta}=' num2str(SD_theta_err)])
%end
%print(gcf,'Ass3-fig6.png','-dpng','-r500') 



%% Question 7
% PURPOSE: Draw the graph with different Cx
% Delete the comment in Question 1 and the comment below
% title(['Cx=Cx*' num2str(i)]) 
% end
% print(gcf,'Ass3-fig7.png','-dpng','-r500') 