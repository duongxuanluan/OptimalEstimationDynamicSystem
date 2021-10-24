function p1= uncertainty(Cx,Ex)
%Draw the uncertainty circle corresponding to Cx and Ex
%PRE-STEP: Find the Cx eigenvalue and eigenvector
[Cx_eigenvector,Cx_eigenvalue]= eig(Cx);

%STEP 1: Make the unit circle around the origin.
th = 0:pi/50:2*pi; % List of the angle points with the step of pi/50   
xunit = cos(th); % x-coordinate of the circle
yunit = sin(th); % y-coordinate of the circle

%STEP 2: Scaling the circle.
a0=sqrt(Cx_eigenvalue(1,1)); a1=sqrt(Cx_eigenvalue(2,2)); % scaling factor
%Note: eigenvalues of cov matrix are always pos
xunit = a0.*xunit; yunit = a1.*yunit; % Scalling

%STEP 3: Rotate the circle
% Rotation matrix is Cx_eigenvector: (x,y)new=Cx_eigvector*(x,y)old 
for i=1:length(xunit) % for each point  
    a(:,i)= Cx_eigenvector*[xunit(i);yunit(i)];% record for each point
end
xunit=a(1,:);yunit=a(2,:); % update rotation

%STEP 4: Shift to Ex
xunit=xunit+Ex(1)*ones(size(xunit));
yunit=yunit+Ex(2)*ones(size(yunit));

%STEP 5: Plot the graph 
p1=plot(xunit,yunit,'k'); % For calling legend

end