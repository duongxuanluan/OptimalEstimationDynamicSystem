function [x_est,alpha, beta] = uLMMSE(x,z)
%Question 3
%unbias Linear MMSE: X_est= alpha*z +beta 
%alpha = cov(x,z)/var(z); beta= E(x)-alpha*E(z)
alpha= (1/(length(x)-1))*sum((x-mean(x)).*(z-mean(z)))./var(z);
beta= E(x)-alpha.*E(z);
x_est=alpha.*z+ beta.*ones(size(z));
end

