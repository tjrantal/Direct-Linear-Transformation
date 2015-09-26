function a=invmodel(name,par)
%INVMODEL calculates the parameters of the inverse camera model, which
%can be used to correct measured image coordinates to correspond a simple
%pinhole model.
%
%Usage:
%   par2=invmodel(name,par)
%
%where
%   name = string that is specific to the camera and the framegrabber.
%          This string must be defined in configc.m
%   par  = camera intrinsic parameters obtained using cacal.m (cacalw.m)
%   par2 = camera intrinsic parameters for correcting coordinates

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland

sys=configc(name);
NDX=sys(1); NDY=sys(2); Sx=sys(3); Sy=sys(4);
Asp=par(1); Foc=par(2);
Cpx=par(3); Cpy=par(4);
Rad1=par(5); Rad2=par(6);
Tan1=par(7); Tan2=par(8);

[dx,dy]=meshgrid(-NDX/40:NDX/40:NDX+NDX/40,-NDY/40:NDY/40:NDY+NDY/40);
cc=imcorr(name,par,[dx(:) dy(:)]);
cx=(cc(:,1)-Cpx)/NDX*Sx/Asp;
cy=(cc(:,2)-Cpy)/NDY*Sy;

%cpx=Cpx/NDX*Sx/Asp;
%cpy=Cpy/NDY*Sy;

%xs=-cpx*0.95; xe=(Sx-cpx)*0.95;
%ys=-cpy*0.95; ye=(Sy-cpy)*0.95;
%px=Sx/40; py=Sy/40;

%[cx,cy]=meshgrid(xs:px:xe,ye:-py:ys); cx=cx(:); cy=cy(:);

r2=cx.*cx+cy.*cy;
delta=Rad1*r2+Rad2*r2.*r2;

Q=1+(4*Rad1*r2+6*Rad2*r2.*r2+8*Tan1*cy+8*Tan2*cx);

dx=cx-(cx.*delta+2*Tan1*cx.*cy+Tan2*(r2+2*cx.*cx))./Q; 
dy=cy-(cy.*delta+Tan1*(r2+2*cy.*cy)+2*Tan2*cx.*cy)./Q; 


r2=dx.*dx+dy.*dy;
  
Tx=[dx.*r2 dx.*r2.*r2 2*dx.*dy r2+2*dx.*dx];
Ty=[dy.*r2 dy.*r2.*r2 r2+2*dy.*dy 2*dx.*dy];
T=[Tx;Ty];
e=[cx-dx;cy-dy];
a=pinv(T)*e;
par=par(:);
a=[par(1:4);a];
