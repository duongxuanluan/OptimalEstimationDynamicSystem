close all
clear all  
clc

load('x_yacht2.mat')
load('z_yacht2.mat')

%% Constants

%Process noise 
sigw_xsi = 10^-2;
sigw_v = 10^-2;
sigw_a = 10^-2;
sigw_t = 8;
sigw_phi = 0.5;

%% INITIALISATION
Cw = zeros(8,8); 
Cw([1,10,19,28,37,46]) = 0.01^2;
Cw(7,7) = 8^2;
Cw(8,8) = 0.5^2; 

Np = 400; %number of particles
x0 = [0; 0; 0; 0; 0; 0; 400; 0];
xb0 = [5000; 10000];
sigma_n = [1; 0.3; 1];
Cn = eye(3).*sigma_n.^2;

%Samples drawn from prior distribution
Xs(1,:) = x0(1) + 100*randn(1, Np);
Xs(2,:) = x0(2) + 100*randn(1, Np);
Xs(3,:) = x0(3) + 2*randn(1, Np);
Xs(4,:) = x0(4) + 2*randn(1, Np);
Xs(5,:) = x0(5) + 0.04*randn(1, Np);
Xs(6,:) = x0(6) + 0.04*randn(1, Np);
Xs(7,:) = x0(7) + 300*randn(1, Np);
Xs(8,:) = x0(8) + 10*randn(1, Np);


z_counter = 1;

for i = 1:length(fi0)
	
    if mod(i,500)~=1 
       Xs = fsys(Xs,[400; fi0(i)],Cw); 
    else
        
    % Update
    % Generate predicted meas. representing p(z(i)|Z(i-1))
    Zs = hmeas(Xs,xb0,Cn);

    % Get uniform distributed rv
    u(1,z_counter) = sum((Zs(1,:)<z(1,z_counter)))/Np;
    u(2,z_counter) = sum((Zs(2,:)<z(2,z_counter)))/Np;
    u(3,z_counter) = sum((Zs(3,:)<z(3,z_counter)))/Np;

    %Residuals
    res = hmeas(Xs,xb0,Cn)-z(:,z_counter)*ones(1,Np);    % Residuals
    res(1,:) = wrapTo180(res(1,:));
    res(3,:) = wrapTo180(res(3,:));

    W = exp(-0.5*sum(res.*(inv(Cn)*res)))';  % Weights
    if (sum(W)==0)
        error("process did not converge"); 
    end

    %Normalizing
    W = W/sum(W); 
    CumW = cumsum(W);
    K_eff(z_counter) = 1/(sum(W.^2));

    % Sample mean
    xest(:,z_counter) = Xs(:,:)*W;

    %Covariance matrix
    Cx{z_counter}= (Xs(1:2,:) - xest(1:2,z_counter)).*W'*(Xs(1:2,:) - xest(1:2,z_counter))';

    % Find an index permutation using golden rule root finding
    for j=1:Np
        R=rand; 
        ja=1; 
        jb = Np;
        while (ja<jb-1)
            jx = floor(jb-0.382*(jb-ja));
            fa = R-CumW(ja);
            fb = R-CumW(jb);
            fxx = R-CumW(jx);
            if (fb*fxx<0)
                ja = jx; 
            else
                jb=jx; 
            end
        end
        ind(j) = jb;
    end
    
    % Resample                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     q
    for j = 1:Np
        Xs(:,j) = Xs(:,ind(j));
    end
        
    %Plot a live figure
%     figure(5)
%     scatter(Xs(1,:),Xs(2,:),5,'x');
%     hold on
%     plot(xest(1,z_counter),xest(2,z_counter),'o','LineWidth',2);
%     ellipse = transform(Cx{z_counter},[xest(1,z_counter);xest(2,z_counter)]);
%     plot(ellipse(1,:), ellipse(2,:),'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',2)
%     if z_counter > 2
%         plot(x(1,i-1000:i),x(2,i-1000:i),'k');
%     end
%     hold off
%     axis equal
%     drawnow

    % Predict
    Xs = fsys(Xs,[400;fi0(i)], Cw);
    
    
    z_counter = z_counter + 1;
end
end

%% Figures
clc
close all;

figure(1)
t = 1:Np;
plot(t, CumW)
title('The cumulative value of the weights')
xlabel('Sample number')
ylabel('Cum value of W')

figure(2) 
hold on
axis equal
grid on
plot(x(1,:),x(2,:))
plot(xest(1,:),xest(2,:))
scatter(xb0(1),xb0(2))
for m = 1:length(z(1,:))
    ellipse = transform(Cx{m},[xest(1,m);xest(2,m)]);
    plot(ellipse(1,:), ellipse(2,:),'Color',[0.8500, 0.3250, 0.0980])
end
title('Path of the vessel')
legend('Real path','Estimated path', 'Lighthouse')
xlabel('x position [m]')
ylabel('y position [m]')
ylim([-3.4 3.1]*10^4)

figure(3)
hold on 
grid on 
t2 = (1:length(z(1,:)));
sgtitle('Test variable u')
ylabel('Normalized deviation from particles')
subplot(3,1,1)
scatter(t2, u(1,:),20,'o','filled','MarkerEdgeColor',[0, 0.4470, 0.7410],'MarkerFaceColor',[0, 0.4470, 0.7410])
subplot(3,1,2)
scatter(t2, u(2,:),20,'d','filled','MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980])
subplot(3,1,3)
scatter(t2, u(3,:),20,'s','filled','MarkerEdgeColor',[0.9290, 0.6940, 0.1250],'MarkerFaceColor',[0.9290, 0.6940, 0.1250])
xlabel('Measurement iteration')

figure(4)
scatter(t2,K_eff)
title('Effective number of particles at each measurement')
xlabel('Measurement iteration')
ylabel('Number of effective particles')

figure(5)
plot(1:length(xest(2,:)),xest(3,:))

%% FUNCTIONS

function fx = transform(C_x,mu_x)
    t = linspace(0,2*pi, 1000);
    x = [cos(t); sin(t)];
    [V,D] = eig(C_x);
    x_scale = (D^0.5)*x;
    x_trans = V*x_scale;
    x_shift = x_trans + mu_x;
    fx = x_shift;
end

