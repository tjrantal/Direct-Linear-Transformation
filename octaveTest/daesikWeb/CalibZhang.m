% CALIBZHANG - computes the intrinsic and extrinsic parameters 
%              from 3D scene points and corresponding 2D image points 
%              with Zhang's method.
%
% Usage:    
%           [K, R, t] = CalibZhang(x, X, m)
%
% Input:
%           x : 2xn, nx2, 3xn or nx3 image points  (n: number of points)
%           X : 2xn, n2x, 3xn or nx3 object points (n: number of points)
%           m : number of planes
%
% Output:
%           K : Camera Intrinsic Matrix
%           R : 3xm matrix, A set of 3x1 Rodrigues's rotation vector
%           t : 3xm matrix, A set of 3x1 translation vector
%
% See also HomoNormDLT, RefineHomo, Rodrigues
%
%
% This code follows the algorithm given by 
% [1] Zhengyou Zhang "A Flexible New Technique for Camera Calibration," 
%     IEEE Transactionson Pattern Analysis and Machine Intelligence, 
%     Vol. 22, No. 11, pp. 1330-1334, 2000.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% Jun. 2010 - Original version.
% Jun. 2011 - Input and output argument type is modified.
% Feb. 2012 - The homography is refined by the LM.


function [K, R, t] = CalibZhang(x, X, m)

%% The Number of Points
xSize = size(x);
XSize = size(X);


%% image points
if (xSize(2) == 2)
    x = x';
elseif (xSize(1) == 3)
    x = x(1:2,:);
elseif (xSize(2) == 3)
    x = x(1:2,:)';
end

%% object points
if (xSize(2) == 2)
    X = X';
elseif (xSize(1) == 3)
    X = X(1:2,:);
elseif (xSize(2) == 3)
    X = X(1:2,:)';
end

if (length(x) ~= length(X))
    error('The length of x and X are not same.');
end


%% Number of points at each plane.
noPnts = length(x)/m;


%% Compute and Stack the Homography
for (n=1:m)
    H_init = HomoNormDLT(X(:,(n-1)*noPnts+1:n*noPnts), x(:,(n-1)*noPnts+1:n*noPnts));
    H(3*n-2:3*n,:) = RefineHomo(X(:,(n-1)*noPnts+1:n*noPnts), x(:,(n-1)*noPnts+1:n*noPnts), H_init);
end


%% Calibration with the Zhang's Method
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


%% Intrinsic Parameters Estimation
v0 = (B(1,2)*B(1,3)-B(1,1)*B(2,3))/(B(1,1)*B(2,2)-B(1,2)^2);
lambda = B(3,3)-(B(1,3)^2 + v0*(B(1,2)*B(1,3)-B(1,1)*B(2,3)))/B(1,1);
fx = sqrt(lambda/B(1,1));
fy = sqrt(lambda*B(1,1)/(B(1,1)*B(2,2)-B(1,2)^2));
skew = -B(1,2)*fx^2*fy/lambda;
u0 = skew*v0/fy - B(1,3)*fx^2/lambda;
K = [fx skew u0; 0 fy v0; 0 0 1];


%% Extrinsic Parameters Estimation
m = 1;
for (n=1:3:sizeR)
    % Tz Sign fixing
    if (H(n+2,3) < 0)
        H(n:n+2,:) = -H(n:n+2,:);
    end
        
    % Rotation Matirx Estimation
    lambda1 = 1/norm(inv(K)*H(n:n+2,1));
    lambda2 = 1/norm(inv(K)*H(n:n+2,2));
    lambda3 = (lambda1 + lambda2)/2;
    
    r1 = lambda1*inv(K)*H(n:n+2,1);
    r2 = lambda2*inv(K)*H(n:n+2,2);
    r3 = cross(r1,r2);

    % Orthogonality Enforcement
    R_temp = [r1 r2 r3];
    [svdU,svdD,svdV] = svd(R_temp);
    R(:,m) = Rodrigues(svdU*eye(3,3)*svdV');
    
    % Translation Vector Estimation
    t(:,m) = lambda3*inv(K)*H(n:n+2,3);
    
    m = m + 1;
end
