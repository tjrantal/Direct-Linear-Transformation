% CALIBPLANE - computes the 3x3 intrinsic matrix from a set of 2D
% homography matrices.
%
% Usage:    
%           K = CalibPlane(H)
%
% Input:
%           H    : 3nx3 set of 3x3 homography matrix
%
% Output:
%           K    : Camera Intrinsic Matrix
%
% This code follows the algorithm given by 
% [1] Hartley and Zisserman "Multiple View Geometry in Computer Vision,"
%     pp.211-212, 2003.
% [2] Zhengyou Zhang "A Flexible New Technique for Camera Calibration," 
%     IEEE Transactionson Pattern Analysis and Machine Intelligence, 
%     Vol. 22, No. 11, pp. 1330-1334, 2000.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% Apr 2008 - Original version.
% Jun 2010 - The matrix B is checked whether it is positive definite
%            for Cholesky Factorization.

function K = CalibPlane(H)

[sizeR, sizeC] = size(H);

V = [];
for (n=1:3:sizeR)
    V = [V;
        H(n  ,1)*H(n  ,2) ...
        H(n,  1)*H(n+1,2)+H(n+1,1)*H(n  ,2) ...
        H(n+1,1)*H(n+1,2) ...
        H(n+2,1)*H(n  ,2)+H(n  ,1)*H(n+2,2) ...
        H(n+2,1)*H(n+1,2)+H(n+1,1)*H(n+2,2) ...
        H(n+2,1)*H(n+2,2)

        H(n  ,1)*H(n  ,1)-H(n  ,2)*H(n  ,2) ...
        H(n  ,1)*H(n+1,1)+H(n+1,1)*H(n  ,1)-(H(n  ,2)*H(n+1,2)+H(n+1,2)*H(n  ,2)) ...
        H(n+1,1)*H(n+1,1)-H(n+1,2)*H(n+1,2) ...
        H(n+2,1)*H(n  ,1)+H(n  ,1)*H(n+2,1)-(H(n+2,2)*H(n  ,2)+H(n  ,2)*H(n+2,2)) ...
        H(n+2,1)*H(n+1,1)+H(n+1,1)*H(n+2,1)-(H(n+2,2)*H(n+1,2)+H(n+1,2)*H(n+2,2)) ...
        H(n+2,1)*H(n+2,1)-H(n+2,2)*H(n+2,2)];
end

%% conic
%   B = [B11, B12, B13
%        B12, B22, B23
%        B13, B23, B33]

[svdU,svdD,svdV] = svd(V);
B_temp = svdV(:,6);
B(1,1) = B_temp(1); B(1,2) = B_temp(2); B(1,3) = B_temp(4);
B(2,1) = B(1,2);    B(2,2) = B_temp(3); B(2,3) = B_temp(5);
B(3,1) = B(1,3);    B(3,2) = B(2,3);    B(3,3) = B_temp(6);

%% solution 1 (Cholesky Factorization)
if (det(B) < 0) 
    B = -B; % B must be positive definite.
end

K = chol(B);
K = inv(K);
K = K/K(3,3);

%% solution 2 (Closed-Form Solution)
% v0 = (B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
% lambda = B(3,3)-(B(1,3)^2 + v0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
% fx = sqrt(lambda/B(1,1));
% fy = sqrt(lambda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
% skew = -B(1,2)*fx^2*fy/lambda;
% u0 = skew*v0/fy - B(1,3)*fx^2/lambda;
% K = [fx skew u0; 0 fy v0; 0 0 1];
