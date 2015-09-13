% HomoNormDLT - computes the 2D homography
%
% Usage:   
%           H = HomoNormDLT(x1, x2)
%
% Input:
%          x1   : 2xn or 3xn homogeneous points set 1
%          x2   : 2xn or 3xn homogeneous points set 2
%         
% Output:
%          H    : 3x3 homography
%
% cf.:
%          x2 ~ H*x1
%
% This code follows the algorithm given by
% [1] Hartley and Zisserman "Multiple View Geometry in Computer Vision,"
%     pp.88-91, pp107-109, 2003.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% July 2008  - Original version.
% May  2011  - Updated

function H = HomoNormDLT(x1, x2)

[r1,c1] = size(x1);
[r2,c2] = size(x2);

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

if (c1 == c2)
    A = zeros(2*c1,9);
    
    for n = 1:c1
        x1_x = x1(1,n); x1_y = x1(2,n);
        x2_x = x2(1,n); x2_y = x2(2,n);
        
        A(2*n-1,:) = [x1_x x1_y 1 0 0 0 -x2_x*x1_x -x2_x*x1_y -x2_x];
        A(2*n  ,:) = [0 0 0 x1_x x1_y 1 -x2_y*x1_x -x2_y*x1_y -x2_y];
    end
    
    [U,D,V] = svd(A);
    
    H = reshape(V(:,9),3,3)';
end

%% Denormalization
H = inv(T2)*H*T1;