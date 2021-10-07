function [pz_x] = pz_x(x,z)
% Question 1
% Define the likelihood function p(z|x), z given x, using the formular (2) from the reader 
global P0 sigma;
P1=1-P0;
pz_x=P0*1/(sigma*sqrt(2*pi))*exp(-((z-x)^2)/(2*sigma^2))+P1*1/(sigma* ...
 sqrt(2*pi))*exp(-((z-2*x)^2)/(2*sigma^2));
end

