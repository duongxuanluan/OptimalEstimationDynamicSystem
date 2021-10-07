function x_ml= ML(z,xrange)
% Question 3
% ML estimator = argmax[p(z|x)]
pz_x_range_z=zeros(length(xrange),1); %allocate array to store pz_x with z1
for it=1:length(xrange)
    pz_x_range_z(it)=pz_x(xrange(it),z); % store data of pz_x with z1
end
[dummy,ind] = max(pz_x_range_z); 
x_ml = xrange(ind);
end

