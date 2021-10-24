function [xfinal,yfinal]=DrawUncertRegion(Cx,meanx)
    [R,Lambda]=eig(Cx);
    A=sqrt(Lambda);
    alpha0=A(1,1);
    alpha1=A(2,2);
    a=1; b=1;          %plotting a unit circle with center=origin 
    x0=0; y0=0;
    t=-pi:0.01:pi;

    ascale=a*alpha0;   %scale horizontal, vertical axises: (with sqrt of cov 
    bscale=b*alpha1;   %matrix eig values)   : convert circle to ellipse
    xscale=x0+ascale*cos(t);
    yscale=y0+bscale*sin(t);

    corRotate=R*[xscale;yscale]; %Rotate the ellipse (with cov matrix eig vectors)
    xrotate=corRotate(1,:);
    yrotate=corRotate(2,:);

    xshift=xrotate+meanx(1);  %shift the centre of ellipse (with the mean value)
    yshift=yrotate+meanx(2);
    xfinal=xshift;
    yfinal=yshift;
    
end