function pos=extinit(sys,data)
%EXTINIT Calculates the external camera parameters using direct decomposition
%of the perspective transformation matrix. The result is used for initialization
%of the nonlinear estimator.
%         
%Usage:
%   pos = extinit(sys,data)  
%
%where           
%   sys  = system configuration information (see configc.m) 
%   data = matrix that contains the 3-D coordinates of the
%          control points (in fixed right-handed frame) and corresponding
%          image observations (in image frame, origo in the upper left
%          corner and the y-axis downwards) 
%          Dimensions: (n x 5) matrix, format: [wx wy wz ix iy]
%          NOTE: in case of a coplanar target wz should be 0 for all points
%   pos  = camera position and orientation [x y z w p r]

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland

NDX=sys(1); NDY=sys(2); Sx=sys(3); Sy=sys(4);
wx=data(:,1); wy=data(:,2); wz=data(:,3);
num=size(data,1);

u0=Sx/2; v0=Sy/2; s=1; f=sys(5);

u = Sx*data(:,4)/NDX;
v = Sy*data(:,5)/NDY;

if any(wz)
   Lu=[wx wy wz 0*u+1 0*u 0*u 0*u 0*u -wx.*u -wy.*u -wz.*u];
   Lv=[0*v 0*v 0*v 0*v wx wy wz 0*v+1 -wx.*v -wy.*v -wz.*v];
   L=reshape([Lu';Lv'],11,2*num)';
   l=reshape([u';v'],2*num,1);

   a=pinv(L)*l;
   a(12)=1;
   F=reshape(a,4,3)';
   lambda=sqrt(sum(F(3,1:3).^2));
   F=F/lambda;

   P=[s*f 0 u0 0;0 f v0 0;0 0 1 0];
   iP=inv(P(:,1:3));
   t=iP*F(:,4);
   R=iP*F(:,1:3);
   [U,S,V]=svd(R);
   Sp=diag([1 1 det(U*V')]);
   R=U*Sp*V';
   pa=asin(R(3,1));
   wa=atan2(-R(3,2)/cos(pa),R(3,3)/cos(pa));
   ra=atan2(-R(2,1)/cos(pa),R(1,1)/cos(pa));
   pos=[t' -[wa pa ra]*180/pi];
else
   Lu=[wx wy 0*u+1 0*u 0*u 0*u -wx.*u -wy.*u];
   Lv=[0*v 0*v 0*v wx wy 0*v+1 -wx.*v -wy.*v];
   L=reshape([Lu';Lv'],8,2*num)';
   l=reshape([u';v'],2*num,1);

   a=pinv(L)*l;
   a(9)=1;
   
   F=reshape(a,3,3)';

   P=[s*f 0 u0 0;0 f v0 0;0 0 1 0];
   iP=inv(P(:,1:3));
   t=iP*F(:,3);
   R=iP*F(:,1:2);
   lambda=sqrt(sum(R(:,1).^2));
   t=t/lambda;
   R=R/lambda;
   R(:,3)=cross(R(:,1),R(:,2));
   [U,S,V]=svd(R);
   Sp=diag([1 1 det(U*V')]);
   R=U*Sp*V';
   pa=asin(R(3,1));
   wa=atan2(-R(3,2)/cos(pa),R(3,3)/cos(pa));
   ra=atan2(-R(2,1)/cos(pa),R(1,1)/cos(pa));
   pos=[t' -[wa pa ra]*180/pi];  
end