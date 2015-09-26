function [f,J]=jacobc(sys,Bs,nbr,p,np,cpar)
%JACOBC Evaluates the camera model and approximates the partial derivates
%with respect to the estimated parameters.
%3-D coordinates and camera parameters. 
%
%Usage:
%   [f,J]=jacobc(sys,Bs,nbr,p)
%
%where       
%   sys  = system configuration information (see configc.m) 
%   Bs   = 3-D configuration of the control points (position, orientation, radius)
%   nbr  = number of points per frame
%   p    = parameter vector (8+N*6 x 1)
%          p(1:8) contains the camera intrinsic parameters
%          p(9...) contains the camera position and orientation for
%          each N images. 
%   f    = synthetic image coordinates
%   J    = Jacobian matrix evaluated at p.

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland

delta=1e-10;
if nargin==4
  np=8;
end
if nargin<6
  cpar=[];
end

n=length(p);
t=(n-np)/6;
ind=[1 cumsum(nbr)+1];
f=frames(sys,Bs,nbr,p,np,cpar);
m=length(f)*2;
J=(zeros(m,n));
p=p(:);

for i=1:np
  q=p;
  q(i)=q(i)+delta;
  g=frames(sys,Bs,nbr,q,np,cpar);
  J(:,i)=(g(:)-f(:))/delta;
end
 
for j=1:t
  for i=1:6
    q=[p(1:np);p(np+1+(j-1)*6:np+j*6)];
    q(i+np)=q(i+np)+delta;
    g=frames(sys,Bs(:,ind(j):ind(j+1)-1),nbr(j),q,np,cpar);
    r=f(ind(j):ind(j+1)-1,:);
    J(ind(j):ind(j+1)-1,np+i+(j-1)*6)=(g(:,1)-r(:,1))/delta;
    J(ind(j)+m/2:ind(j+1)-1+m/2,np+i+(j-1)*6)=(g(:,2)-r(:,2))/delta;
  end
end
J=sparse(J);
f=f(:);
