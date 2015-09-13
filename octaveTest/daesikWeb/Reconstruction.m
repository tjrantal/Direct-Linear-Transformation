% RECONSTRUCTION - calculates the point in space with two projection
% matrices and the corresponding points.
%
% Usage:    
%           X = Reconstruction(P1, P2, x1, x2)
%
% Input:
%           P1  : 3x4 projection matrix 1
%           P2  : 3x4 projection matrix 2
%           x1  : 2xn (or 3xn) image point set 1
%           x2  : 2xn (or 3xn) image point set 2
%
% Output:
%           X   : 4xn point set in space
%
% This code follows the algorithm given by 
% [1] "Summary Note 2006-SN-004-EN 3D Reconstruction, ISRC, SKKU"
%     available at http://www.daesik80.com
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% May 2008  - Original version.

function X = Reconstruction(P1, P2, x1, x2)

[r1,c1] = size(x1);
[r2,c2] = size(x2);

if (r1 == 3)
    x1(1,:) = x1(1,:)./x1(3,:);
    x1(2,:) = x1(2,:)./x1(3,:);
    x2(1,:) = x2(1,:)./x2(3,:);
    x2(2,:) = x2(2,:)./x2(3,:);
end

if (c1 == c2)
    for (n=1:c1)
        M(1,1) = P1(1,1) - x1(1,n)*P1(3,1);
        M(1,2) = P1(1,2) - x1(1,n)*P1(3,2);
        M(1,3) = P1(1,3) - x1(1,n)*P1(3,3);
        M(1,4) = P1(1,4) - x1(1,n)*P1(3,4);
        M(2,1) = P1(2,1) - x1(2,n)*P1(3,1);
        M(2,2) = P1(2,2) - x1(2,n)*P1(3,2);
        M(2,3) = P1(2,3) - x1(2,n)*P1(3,3);
        M(2,4) = P1(2,4) - x1(2,n)*P1(3,4);
        M(3,1) = P2(1,1) - x2(1,n)*P2(3,1);
        M(3,2) = P2(1,2) - x2(1,n)*P2(3,2);
        M(3,3) = P2(1,3) - x2(1,n)*P2(3,3);
        M(3,4) = P2(1,4) - x2(1,n)*P2(3,4);
        M(4,1) = P2(2,1) - x2(2,n)*P2(3,1);
        M(4,2) = P2(2,2) - x2(2,n)*P2(3,2);
        M(4,3) = P2(2,3) - x2(2,n)*P2(3,3);
        M(4,4) = P2(2,4) - x2(2,n)*P2(3,4);
        
        [svdU,svdD,svdV] = svd(M);
      
        X(1,n) = svdV(1,4)/svdV(4,4);
        X(2,n) = svdV(2,4)/svdV(4,4);
        X(3,n) = svdV(3,4)/svdV(4,4);
        X(4,n) = svdV(4,4)/svdV(4,4);
    end
end