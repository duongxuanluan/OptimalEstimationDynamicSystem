function [x_est,alpha] = LMMSE(x,z)
% Question 2
% Linear MMSE: X_est= alpha*z where alpha = E(xz)/E(z^2)
alpha=E(x.*z)./E(z.*z);
x_est=alpha*z;
end

