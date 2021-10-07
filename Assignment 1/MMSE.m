function x_mmse= MMSE(z,xrange)
% Question 3
% MMSE estimator= E(x|z) using numerical integration
dx=xrange(2)-xrange(1); % sample period of
dummy=zeros(length(xrange),1); % dummy allocate array to store MMSE values for numerical intergration
for it=1:length(xrange)
   dummy(it)= xrange(it)*px_z(xrange(it),z,xrange)*dx; %rectangular area 
end
x_mmse=sum(dummy); %numerical integration - Euler method
end

