function J = fpf (y)
load hyddata.mat; % Load the dataset (measurements in Z)
I = length (Z); % Length of the sequence
R1 = 105.78; % Friction constant (cm^5/s^2)
R2 = 84.532; % Friction constant (cm^5/s^2)
qmin = 0; % Minimal input flow
delta = 5; % Sampling period (s)
C = 420; % Capacity of tank (cm^2)
sigma_v = 0.04; % Standard deviation sensor noise (cm)
Ncond = 1000; % Number of particles
M = 3; % Dimension of state vector
N = 2; % Dimension of measurement vector
% Set the design parameters
qmax = y(1); % Maximum input flow
Pup = y(2); % Transition probability up
Pdn = y(3); % Transition probability down
% Initialisation
hmax = 0; hmin = 0; % Margins of levels at i = 0
H = eye(N,M); % Measurement matrix
invCv = inv(sigma_v^2 * eye(N)); % Inv. cov. of sensor noise
% Generate the samples
Xs(1,:) = hmin + (hmax-hmin)*rand(1,Ncond);
Xs(2,:) = hmin + (hmax-hmin)*rand(1,Ncond);
Xs(3,:) = qmin + (qmax-qmin)*(rand(1,Ncond)> 0.5);
for i = 1:I
    % Generate predicted meas. representing p(z(i)|Z(i-1))
    Zs = H*Xs + sigma_v*randn(2,Ncond);
    
    % Get uniform distributed rv
    u(1,i)= sum((Zs(1,:) < Z(1,i)))/Ncond;
    u(2,i)= sum((Zs(2,:) < Z(2,i)))/Ncond;
    
    % Update
    res = H*Xs-Z(:,i)*ones(1,Ncond); % Residuals
    W = exp(-0.5*sum(res.*(invCv*res)))'; % Weights
    if (sum(W) == 0), error('process did not converge'); end
    W = W/sum(W); CumW = cumsum(W);
    
    xest(:,i) = Xs(:,:)*W; % Sample mean
    
    % Find an index permutation using golden rule root finding
    for j = 1:Ncond
        R = rand; ja = 1; jb = Ncond;
        while (ja < jb1)
            jx = floor(jb0:382*(jbja));
            fa = R-CumW(ja); fb = R-CumW(jb); fxx = R-CumW(jx);
            if (fb*fxx < 0), ja = jx; else, jb = jx; end
        end
        ind(j) = jb;
    end
    
    % Resample
    for j = 1:Ncond, Ys(:,j) = Xs(:,ind(j)); end
   
    % Predict
    Tdn = (rand(1,Ncond)<Pdn); % Random transitions
    Tup = (rand(1,Ncond)<Pup); % idem
    kdn = find((Ys(3,:) == qmax) & Tdn); % Samples going down
    kup = find((Ys(3,:) == qmin) & Tup); % Samples going up
    Ys(3,kdn) = qmin; % Turn input flow off
    Ys(1,kdn) = Ys(1,kdn) + ... % Randomize level 1
    (qmax-qmin)*delta*rand(1,length(kdn))/C;
    Ys(3,kup) = qmax; % Turn input flow on
    Ys(1,kup) = Ys(1,kup)- ... % Randomize level 1
    (qmax-qmin)*delta*rand(1,length(kup))/C;
    Xs = f(Ys,R1,R2); % Update samples
end
e = I/10; % Expected number of rv in one bin
% Get histograms (10 bins) and calculate test variables
for i = 1:2, b = hist(u(i,:)); c(i) = sum((be).^2/e); end
J = sum(c); % Full test (chi-square with 19 DoF)
return

