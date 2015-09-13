% HOMODLT - computes the 2D homography
%
% Usage:   
%           H = HomoDLT(x1, x2)
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
%     pp.88-89, 2003.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% May 2008  - Original version.
% May 2011  - Updated

function H = HomoDLT(x1, x2)

[r1,c1] = size(x1);
[r2,c2] = size(x2);

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