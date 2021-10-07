function px_z = px_z(x,z,xrange)
%Question 2
%Define function p(x|z) from p(z|x),p(x) and p(z)
%p(x|z)= p(z|x)*p(x)/p(z)
px_z=pz_x(x,z)*px(x)./pz(z,xrange);
end

