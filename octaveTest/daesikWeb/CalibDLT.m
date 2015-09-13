% CalibDLT - computes camera projection matrix from 3D scene points and
% corresponding 2D image points with DLT (direct linear transformation).
%
% Usage:
%           [P, rperr] = CalibDLT(x, X)
%
% Input:
%           x : 2xn, nx2, 3xn or nx3 image points
%           X : 3xn, nx3, 4xn or nx4 object points
%
% Output:
%           P : 3x4 camera projection matrix
%           rperr : re-projection error
%
% cf.:
%           x = P*X
%
% This code follows the algorithm given by
% [1] E. Trucco and A. Verri, "Introductory Techniques for 3-D
%     Computer Vision," pp. 132-134, 1998.
%
% Kim, Daesik, Ph.D
% Graduated from Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea, in 2013.
% E-mail  : daesik80@skku.edu, daesik80@gmail.com
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% June 2008  - Original version.
% Mar. 2014  - A re-projection error is added to the output.


function [P, rperr] = CalibDLT(x, X)


%% The Number of Points
xSize  = size(x);
XSize  = size(X);
noPnts = length(x);


%% image points
if (xSize(2) == 2)
    x = x';
elseif (xSize(1) == 3)
    x = x(1:2,:);
elseif (xSize(2) == 3)
    x = x(:,1:2)';
end


%% object points
if (XSize(2) == 3)
    X = X';
elseif (XSize(1) == 4)
    X = X(1:3,:);
elseif (XSize(2) == 4)
    X = X(:,1:3)';
end


%% Compute the Camera Projection Matrix
A = [X(1,:)'            X(2,:)'             X(3,:)'             ones(noPnt,1) ...
     zeros(noPnt,1)     zeros(noPnt,1)      zeros(noPnt,1)      zeros(noPnt,1) ...
     -x(1,:)'.*X(1,:)'  -x(1,:)'.*X(2,:)'   -x(1,:)'.*X(3,:)'   -x(1,:)';
     zeros(noPnt,1)     zeros(noPnt,1)      zeros(noPnt,1)      zeros(noPnt,1) ...
     X(1,:)'            X(2,:)'             X(3,:)'             ones(noPnt,1) ...
     -x(2,:)'.*X(1,:)'  -x(2,:)'.*X(2,:)'   -x(2,:)'.*X(3,:)'   -x(2,:)'];
     
[U,D,V] = svd(A);
P = reshape(V(:,12),4,3)';


%% Re-Projection
x_rp = P*[X; ones(1,length(X))];
x_rp(1,:) = x_rp(1,:)./x_rp(3,:);
x_rp(2,:) = x_rp(2,:)./x_rp(3,:);


%% Distance between the re-projected points and the measured points
dist_u = x_rp(1,:) - x(1,:);
dist_v = x_rp(2,:) - x(2,:);
    
    
%% Re-projection Error
rperr_u = sqrt(dot(dist_u,dist_u)/noPnts);
rperr_v = sqrt(dot(dist_v,dist_v)/noPnts);
rperr = sqrt(rperr_u*rperr_u + rperr_v*rperr_v);
