function x_mmae = MMAE(z,xrange)
% Question 3
% MMAE estimator = med[p(x|z)]
px_z_range=zeros(length(xrange),1); %allocate array to store px_z with z1
for it=1:length(xrange)
    px_z_range(it)=px_z(xrange(it),z,xrange); % store data of px_z with z1
end
ind = find((cumsum(px_z_range) ./ sum(px_z_range))>0.5);
x_mmae = xrange(ind(1));
end

