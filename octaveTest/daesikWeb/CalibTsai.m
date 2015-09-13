% CalibTsai - computes the intrinsic and extrinsic parameters with Tsai's method.
%
% Usage:    
%           [K, R, t, rperr] = CalibTsai(x, X, u0, v0)
%
% Input:
%           x : 3xn homogeneous image points  (n: # of points)
%           X : 3xn homogeneous object points (n: # of points)
%
% Output:
%           K : Camera Intrinsic Matrix
%           R : 3x3 Rotation Matrix
%           t : 3x1 Translation Vector
%
%
% This code follows the algorithm given by 
% [1] R.Y. Tsai, "A Versatile Camera Calibration Technique for 
%                 High-Accuracy 3D Machine Vision Metrology 
%                 Using Off-the-Shelf TV Cameras and Lenses,"
%                 IEEE Journal of Robotics Automation, 
%                 vol. RA-3, no. 4, pp. 323-344, 1987. 
% [2] E. Trucco and A. Verri, "Introductory Techniques for 3-D
%     Computer Vision," pp. 127-130, 1998.
%
% Kim, Daesik, Ph.D
% Gradudated from Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea, in 2013.
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% May 2011 - Original version.
% Mar 2014 - Bugs are fixed, and output (re-projection error) is added.


function [K, R, t, rperr] = CalibTsai(x, X, u0, v0)

noPnts = length(x);
x_bar = zeros(size(x));

x_bar(:,1) = x(:,1) - u0;
x_bar(:,2) = x(:,2) - v0;

M = [ x_bar(:,1).*X(:,1)  x_bar(:,1).*X(:,2)  x_bar(:,1).*X(:,3)  x_bar(:,1) ...
     -x_bar(:,2).*X(:,1) -x_bar(:,2).*X(:,2) -x_bar(:,2).*X(:,3) -x_bar(:,2)];
  
[U,D,V] = svd(M);


%% Compute scale and aspect ratio
scale = 1/sqrt(V(1,8)^2 + V(2,8)^2 + V(3,8)^2)
alpha = scale*sqrt(V(5,8)^2 + V(6,8)^2 + V(7,8)^2);

if (scale*x_bar(1,2)*(V(1,8)*X(1,1) + V(2,8)*X(1,2) + V(3,8)*X(1,3) + V(4,8)) < 0)
    scale = -scale;
end


%% Compute r1, r2, tx, and ty
r1(1,1) = scale*V(5,8)/alpha;
r1(1,2) = scale*V(6,8)/alpha;
r1(1,3) = scale*V(7,8)/alpha;
r2(1,1) = scale*V(1,8);
r2(1,2) = scale*V(2,8);
r2(1,3) = scale*V(3,8);
tx      = scale*V(8,8)/alpha;
ty      = scale*V(4,8);


%% Compute a rotation matrix
r3 = cross(r1, r2);

R = [r1; r2; r3];
[U, D, V] = svd(R);
R = U*V';


%% Compute fx, fy, and tz
Xc = R(1,1).*X(:,1) + R(1,2).*X(:,2) + R(1,3).*X(:,3) + tx;
Yc = R(2,1).*X(:,1) + R(2,2).*X(:,2) + R(2,3).*X(:,3) + ty;
Zc_tz = R(3,1).*X(:,1) + R(3,2).*X(:,2) + R(3,3).*X(:,3); % Zc - tz

A = zeros(length(x_bar)*2, 2);
B = zeros(length(x_bar)*2, 1);

A(1:2:2*length(x_bar), :) = [Xc, -x_bar(:,1)];
A(2:2:2*length(x_bar), :) = [Yc/alpha, -x_bar(:,2)];

B(1:2:2*length(x_bar), 1) = [x_bar(:,1).*Zc_tz];
B(2:2:2*length(x_bar), 1) = [x_bar(:,2).*Zc_tz];

Y = pinv(A)*B;
fx = Y(1);
tz = Y(2);
fy = fx/alpha;

K = [fx 0 u0; 0 fy v0; 0 0 1];
t = [tx; ty; tz];


%% Re-Projection
P = K*[R t];
x_rp = P*[X'; ones(1,length(X))];
x_rp(1,:) = x_rp(1,:)./x_rp(3,:);
x_rp(2,:) = x_rp(2,:)./x_rp(3,:);


%% Distance between the re-projected points and the measured points
dist_u = x_rp(1,:) - x(:,1)';
dist_v = x_rp(2,:) - x(:,2)';
    
    
%% Re-projection Error
rperr_u = sqrt(dot(dist_u,dist_u)/noPnts);
rperr_v = sqrt(dot(dist_v,dist_v)/noPnts);
rperr = sqrt(rperr_u*rperr_u + rperr_v*rperr_v);
