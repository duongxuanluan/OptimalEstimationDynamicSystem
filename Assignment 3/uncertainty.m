function [p1,Cx_eigenvalue,Cx_eigenvector,xunit,yunit]= uncertainty(Cx,Ex,Color1,Color2)
%Exercise 3- Question 1
%Draw the uncertainty circle corresponding to Cx and Ex
%PRE-STEP: Find the Cx eigenvalue and eigenvector
[Cx_eigenvector,Cx_eigenvalue]= eig(Cx);

%STEP 1: Make the unit circle around the origin.
th = 0:pi/50:2*pi; % List of the angle points with the step of pi/50   
xunit = cos(th); % x-coordinate of the circle
yunit = sin(th); % y-coordinate of the circle
% figure (1)
% hold on
% a1=plot(xunit,yunit,'r-');

%STEP 2: Scaling the circle.
a0=sqrt(Cx_eigenvalue(1,1)); a1=sqrt(Cx_eigenvalue(2,2)); % scaling factor
%Note: eigenvalues of cov matrix are always pos
xunit = a0.*xunit; yunit = a1.*yunit; % Scalling
% a2=plot(xunit,yunit,'g-');

%STEP 3: Rotate the circle
% Rotation matrix is Cx_eigenvector: (x,y)new=Cx_eigvector*(x,y)old 
for i=1:length(xunit) % for each point  
    a(:,i)= Cx_eigenvector*[xunit(i);yunit(i)];% record for each point
end
xunit=a(1,:);yunit=a(2,:); % update rotation
% a3=plot(xunit,yunit,'b-');

%STEP 4: Shift to Ex
xunit=xunit+Ex(1)*ones(size(xunit));
yunit=yunit+Ex(2)*ones(size(yunit));
% a4=plot(xunit,yunit,'b-');
% a5=plot(Ex(1),Ex(2), '*');
% legend( 'Step 1: Unit circle', 'Step 2: Scale circle',...
%     'Step 3: Rotate circle', 'Step 4: Circle shifted to Ex', 'Ex', ...
%     'Location', 'southeast');
% xlabel('\xi');ylabel('\eta')
% title('Steps to draw uncertainty circle')
% print(gcf,'Ass3-fig1.png','-dpng','-r500'); 
% axis equal
% hold off


%STEP 5: Plot the graph 
% figure (2)
p1=plot(xunit,yunit, Color1 ,Ex(1),Ex(2),Color2); % For calling legend
axis equal

end

