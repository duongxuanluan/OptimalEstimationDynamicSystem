function [px] = px(x)
% Question 1
% Function of prior pdf p(x) using the formula (3) from the reader 
global xmin xmax beta;
alpha=1/2*(xmax-xmin); %define alpha
mean=1/2*(xmax+xmin); %define mean
px=beta/(2*alpha*gamma(1/beta))*exp(-(abs(x-mean)/alpha)^beta); 
end

