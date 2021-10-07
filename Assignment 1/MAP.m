function x_map = MAP(z,xrange)
% Question 3
% MAP estimator=argmax p(x|z)
px_z_range=zeros(length(xrange),1); %allocate array to store px_z with z1
for it=1:length(xrange)
    px_z_range(it)=px_z(xrange(it),z,xrange); % store data of px_z with z1
end
[dummy,ind] = max(px_z_range);
x_map = xrange(ind); 
end

