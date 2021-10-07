function [pz] = pz(z,xrange)
%Question 2
%Define function p(z) from p(z|x) and p(x)
%p(z)= int(p(x,z)dx)=int(p(z|x)*p(x)dx)
dx=xrange(2)-xrange(1); % dx is sample period
pzrange=zeros(length(xrange),1); % allocate array to store pz for numerical intergration
for it=1:length(xrange)
    pzrange(it)= px(xrange(it)).*pz_x(xrange(it),z)*dx; %rectangulars area for each value of x step
end
pz=sum(pzrange);%numerical integration
end


