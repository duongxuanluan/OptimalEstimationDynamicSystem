%% EXERCISE 4:   propagation of uncertainty; prediction
% Student: Luan Duong - Student ID:S2236117
% Date: 1st December 2020
% Objectives of the exercise:
%   - the propagation of mean and covariance matrix in dynamic systems
%   - how to use this insight for prediction

% Summary of question: Dead reckoning - location of the ship
%   - Xsi_est(i):location, v(i):velocity =>Xsi_est(i+l) where l: steps ahead
%   - Model: x(i+1)=F*x(i)+w(i), w(i) is white noise 
%   - State vector: x(i)= (Xsi(i),v(i),a(i))

clear all; close all % start 
% Given Case: Using system identification for a(i+1)= F1*a(i)+w1
F1=[-0.0595 -0.1530;-0.0813 -0.176]; 
Cw1=[0.1177 -0.0026; -0.0026 0.0782]*10^(-3);

%% Question 1: 
% PURPOSE: - Determine F, Cw and create state vector x(i)
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

% 3. Create state vector x
load log.mat
Xsi=xsi; 
v=zeros(2,length(Xsi)-1); a=zeros(2,length(v)-1); %Saving places
x=zeros(6,length(a)-1);
for it=1:length(Xsi)-1
    v(:,it)=Xsi(:,it+1)-Xsi(:,it); % calculate v
end
for it=1:length(v)-1
    a(:,it)=v(:,it+1)-v(:,it); % calculate a
end
for it=1:length(a)
    x(:,it)=[Xsi(:,it);v(:,it);a(:,it)]; % calculate x
end

%% Question 2: 
% PURPOSE: i=10, predict position Xi and x at j=i+90
i=10; % current time
%Active the below comment for question 5
% i_list=[30 50 70 90];
% for counter=1:length(i_list)
%     i=i_list(counter);
%     figure (2)
%     subplot (2,2,counter)
    
j=100; % Predict position 
x_pred=zeros(6,j-i);% saving place for x
x_pred(:,1)=x(:,i+1); % Saving in matlab starts from 1
% n-step prediction:
for it=2:length(x_pred)
    x_pred(:,it)=F*x_pred(:,it-1); 
end


%% Question 3
% PURPOSE: Plot the measured position vs predicted position 
p1=plot (x(1,:),x(2,:),'r.',x_pred(1,:),x_pred(2,:),'b.'); 
xlabel('x'); ylabel('y');
title('The measured position vs The predicted position');
legend('the measured postition', 'the predicted position');
hold on % draw the question 4
print(gcf,'Ass4-fig1.png','-dpng','-r500');

%% Question 4
% PURPOSE: Recording all the Cx and draw the uncertainty area
Cx_pred=zeros(6,6,j-i);% saving place for Cx, ininital Cx=0
for it=2:(j-i)
    Cx_pred(:,:,it)=F*Cx_pred(:,:,it-1)*F.'+Cw;
    p2=uncertainty(Cx_pred(1:2,1:2,it),x_pred(:,it),'g-.');
end
legend([p1(1) p1(2) p2],'the measured postition', 'the predicted position', 'uncertainty area');
title('The measured position vs The predicted position and uncertainty area');
print(gcf,'Ass4-fig2.png','-dpng','-r500');

%% Question 5
% Purpose: Test with different i 
% Go to question 2 and activate the comment
%title(['i=' num2str(i)])
%end
%print(gcf,'Ass4-fig3.png','-dpng','-r500');
    




