function [pos,iter,res,er,C]=extcal(name,data,cpar)
%EXTCAL Calculates the extrinsic camera parameters for a single image 
%from 3D - 2D correspondences.
%
%Usage:
%   [pos,iter,res,er,C]=extcal(name,data,cpar)         
%
%where
%   name = string that is specific to the camera and the framegrabber.
%          This string must be defined in configc.m
%   data = matrix that contains the 3-D coordinates of the
%          control points (in fixed right-handed frame) and corresponding
%          image observations (in image frame, origin in the upper left
%          corner and the y-axis downwards), and normal vectors of the 
%          object surface. 
%          dimensions: (n x 8) matrices, row format: [wx wy wz ix iy nx ny nz]
%   cpar = intrinsic parameters of the camera (from calibration)
%   pos  = camera position and orientation
%   iter = number of iterations used
%   res  = residual (sum of squared errors)
%   er   = remaining error for each point
%   C    = error covariance matrix of the estimated parameters

%   Version 3.0  10-17-00
%   Janne Heikkila, University of Oulu, Finland

if ~isstr(name)
  error('The first argument should be the camera type');
end
sys=configc(name);
n=length(data);

ipos=extinit(sys,data(:,1:5));

Bs=cinit(sys,data(:,1:3),data(:,6:8));
[pos,iter,res,er,J,succ]=lmoptc(sys,Bs,data(:,4:5),n,ipos(:),1,0,cpar);
C=full(inv(J'*J))*var(er);
disp(sprintf('Standard error in pixels: %f',std(er(:))));
disp(sprintf('Standard deviation of the estimated extrinsic parameters:'));
disp(sprintf('%.5f ',sqrt(diag(C)')));
