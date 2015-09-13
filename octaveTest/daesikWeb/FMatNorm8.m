% FMatNorm8 - computes fundamental matrix with normalized 8-point algorithm
%
% Usage:
%           [F, e1, e2] = FMatNorm8(x1, x2)
%
% Input:
%           x1 : 3xn homogeneous image points of the first image
%           x2 : 3xn homogeneous image points of the second image
%
% Output:
%           F  : 3x3 fundamental matrix
%           e1 : epipole of the first image
%           e2 : epipole of the second image
%
% cf.:
%           x2'*F*x1 = 0
%           F*e1 = 0
%           F'*e2 = 0
%
% This code follows the algorithm given by
% [1] Hartley and Zisserman "Multiple View Geometry in Computer Vision,"
%     pp.281-282, 2003.
% [2] E. Trucco and A. Verri "Introductory Techniques for 3-D Computer Vision,"
%     pp.155-157, 1998.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% Jun. 2008 - Original version.
% Jul. 2008 - centroid calculation bug fixed.

function [F,e1,e2] = FMatNorm8(x1, x2)

%% Number of Points
noPnt = length(x1);

%% centroids of the points
centroid1 = mean(x1(1:2,:)')';
centroid2 = mean(x2(1:2,:)')';

%% Shift the origin of the points to the centroid
x1(1,:) = x1(1,:) - centroid1(1);
x1(2,:) = x1(2,:) - centroid1(2);

x2(1,:) = x2(1,:) - centroid2(1);
x2(2,:) = x2(2,:) - centroid2(2);

%% Normalize the points so that the average distance from the origin is equal to sqrt(2).
averagedist1 = mean(sqrt(x1(1,:).^2 + x1(2,:).^2));
averagedist2 = mean(sqrt(x2(1,:).^2 + x2(2,:).^2));

scale1 = sqrt(2)/averagedist1;
scale2 = sqrt(2)/averagedist2;

x1(1:2,:) = scale1*x1(1:2,:);
x2(1:2,:) = scale2*x2(1:2,:);

%% similarity transform 1
T1 = [scale1    0       -scale1*centroid1(1)
      0         scale1  -scale1*centroid1(2)
      0         0       1      ];

%% similarity transform 2
T2 = [scale2    0       -scale2*centroid2(1)
      0         scale2  -scale2*centroid2(2)
      0         0       1      ];

%% Compute the Fundamental Matrix
A = [x2(1,:)'.*x1(1,:)'  x2(1,:)'.*x1(2,:)'  x2(1,:)' ...
     x2(2,:)'.*x1(1,:)'  x2(2,:)'.*x1(2,:)'  x2(2,:)' ...
     x1(1,:)'            x1(2,:)'            ones(noPnt,1) ];       

[U,D,V] = svd(A,0);
F = reshape(V(:,9),3,3)';

%% Rank 2 Constraint Enforcement
[U,D,V] = svd(F,0);
D(3,3) = 0;
F = U*D*V';

%% Denormalization
F = T2'*F*T1;

%% Epipoles
[U,D,V] = svd(F,0);
e1 = V(:,3)/V(3,3);
e2 = U(:,3)/U(3,3);