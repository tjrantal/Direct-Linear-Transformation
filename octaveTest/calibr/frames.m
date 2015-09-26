function f=frames(sys,Bs,nbr,p,np,cpar)
%FRAMES Interface between camera model and optimization routine. 
%Allows multiple to be used.
%
%Usage:
%   f = frames(sys,Bs,nbr,p,np,cpar)         
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

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland

nr=size(Bs,2);

npos=(length(p)-np)/6;
if npos~=length(nbr)
  error('Conflict in the number of camera positions');
end

if size(Bs,1)~=16
  error('Data matrix should have 16 rows');
end

index=[1 cumsum(nbr)+1];
if np>0
  par=p(1:np);
else
  par=cpar;
end
f=[];

for i=1:npos
  si=(i-1)*6+np+1; ei=i*6+np;
  pos=p(si:ei);
  cp=cmodel(sys,Bs(:,index(i):index(i+1)-1),p(si:ei),par);
  f=[f;cp];
end

  
