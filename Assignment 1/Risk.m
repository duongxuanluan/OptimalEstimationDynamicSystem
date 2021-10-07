function R= Risk(C,xrange,x_est,z)
%Question 4 
% Function of risk: R= intergral(C*p(x|z)*dx)
dx=xrange(2)-xrange(1);
R_list=[];
for it=1:length(xrange)
 R_list(it)=C(x_est,xrange(it))*px_z(xrange(it),z,xrange)*dx;
end 
R=sum(R_list);
end

