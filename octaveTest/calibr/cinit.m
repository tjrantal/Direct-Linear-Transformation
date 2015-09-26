function Bs=cinit(sys,data,snorm)
%CINIT forms the matrix containing the control point position, orientation and
%geometry. This is an internal function used by CIRCAL, CIRCALW, and EXTCAL.

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland


n=size(data,1);
r=sys(6);
v3=snorm;
vv=randn(3,1); vv=vv/norm(vv);
v2=ones(n,1)*vv'-v3*[vv vv vv].*v3;
nnn=sqrt(sum(v2'.^2))';
v2=v2./[nnn nnn nnn];
v1=[v2(:,2).*v3(:,3)-v2(:,3).*v3(:,2) ...
    v2(:,3).*v3(:,1)-v2(:,1).*v3(:,3) ...
    v2(:,1).*v3(:,2)-v2(:,2).*v3(:,1)];

x=data(:,1);
y=data(:,2);
z=data(:,3);

B11=x.^2-r^2*(v1(:,1).^2+v2(:,1).^2);
B12=x.*y-r^2*(v1(:,1).*v1(:,2)+v2(:,1).*v2(:,2));
B13=x.*z-r^2*(v1(:,1).*v1(:,3)+v2(:,1).*v2(:,3));
B22=y.^2-r^2*(v1(:,2).^2+v2(:,2).^2);
B23=y.*z-r^2*(v1(:,2).*v1(:,3)+v2(:,2).*v2(:,3));
B33=z.^2-r^2*(v1(:,3).^2+v2(:,3).^2);

Bs=[B11 B12 B13 x B12 B22 B23 y B13 B23 B33 z x y z x*0+1]';
