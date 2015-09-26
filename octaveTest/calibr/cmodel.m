function ic=cmodel(sys,Bs,pos,par)
%CMODEL Camera model that generates synthetic image coordinates from 
%3-D coordinates and camera parameters. 
%
%Usage:
%   ic = cmodel(sys,Bs,pos,par)
%
%where       
%   sys  = system configuration information (see configc.m) 
%   Bs   = 3-D configuration of the control points (position, orientation, radius)
%   pos  = camera position and orientation
%   par  = camera intrinsic parameters:
%          par(1) = scale factor ~1
%          par(2) = effective focal length
%          par(3:4) = principal point
%          par(5:6) = radial distortion coefficients
%          par(7:8) = tangential distortion coefficients
%   ic   = synthetic image coordinates

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland


NDX=sys(1); NDY=sys(2); Sx=sys(3); Sy=sys(4);

if length(pos)~=6
  error('Position vector should contain [x y z w p r].');
end

if length(par)~=8
  error('Parameter vector should be 1 x 8 matrix.');
end

if size(Bs,1)~=16
  error('Control point matrix should have 16 rows');
end

n=size(Bs,2);
Asp=par(1); Foc=par(2);
Cpx=par(3); Cpy=par(4);
Rad1=par(5); Rad2=par(6);
Tan1=par(7); Tan2=par(8);

M=eye(4);
wa=pos(4)*pi/180;
pa=pos(5)*pi/180;
ra=pos(6)*pi/180;
cw=cos(wa); sw=sin(wa);
cp=cos(pa); sp=sin(pa);
cr=cos(ra); sr=sin(ra);

M(1,:)=[cr*cp -sr*cw+cr*sp*sw sr*sw+cr*sp*cw pos(1)];
M(2,:)=[sr*cp cr*cw+sr*sp*sw -cr*sw+sr*sp*cw pos(2)];
M(3,:)=[-sp cp*sw cp*cw pos(3)];

P=[Foc 0 0 0;0 Foc 0 0;0 0 1 0];

A=P*M;
B1=Bs(1:4,:);
B2=Bs(5:8,:);
B3=Bs(9:12,:);
B4=Bs(13:16,:);

H=B1*A(3,1)+B2*A(3,2)+B3*A(3,3)+B4*A(3,4);
q1=A(1,:)*H;
q2=A(2,:)*H;
q3=A(3,:)*H;
x=(q1./q3)';
y=(q2./q3)';

r2=x.*x+y.*y;
delta=3*Rad1*r2+5*Rad2*r2.*r2;

xn=x.*(1+delta)+6*Tan1*x.*y-Tan2*(r2-6*x.*x); 
yn=y.*(1+delta)-Tan1*(r2-6*y.*y)+6*Tan2*x.*y; 
Q=4*Rad1*r2+6*Rad2*r2.*r2+8*Tan1*y+8*Tan2*x+1;
xn=xn./Q; yn=yn./Q;
  
ic=NDX*Asp*xn/Sx+Cpx;
ic(:,2)=NDY*yn/Sy+Cpy;

