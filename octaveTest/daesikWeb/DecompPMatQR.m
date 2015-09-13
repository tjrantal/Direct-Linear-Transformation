% DecompPMatQR - decompose the camera projection matrix P into intrinsic
%                matrix K, rotation matrix R, and translation vector t.
%
% Usage:
%           [K,R,t] = DecompPMatQR(P)
%
% Input:
%           P : 3x4 camera projection matrix
%
% Output:
%           K : 3x3 intrinsic matrix
%           R : 3x3 rotation matrix
%           t : 3x1 translation vector
%
% cf.:
%           P = K*[R t]
%
% This code follows the algorithm given by
% [1] "Summary Note 2006-SN-003-EN_Camera Calibration, ISRC, SKKU"
%     available at http://www.daesik80.com
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% June 2008  - Original version.
% Aug  2008  - checks the sign of the intrinsic parameters.

function [K,R,t] = DecompPMatQR(P)

Q = inv(P(1:3, 1:3));
[R,K] = qr(Q);

%% Sign Checking
if (K(1,1) < 0)
    S = [-1 0 0;0 1 0;0 0 1];
    R = R*S;
    K = S*K;
end

if (K(2,2) < 0)
    S = [1 0 0;0 -1 0;0 0 1];
    R = R*S;
    K = S*K;
end

if (K(3,3) < 0)
    S = [1 0 0;0 1 0;0 0 -1];
    R = R*S;
    K = S*K;
end

%% Translation Vector
t = K*P(1:3,4);

%% if R is a rotation matrix, then det(R)=1.
if det(R)< 0 
    t = -t;
    R = -R;
end

%% Rotation Matrix
R = inv(R);

%% Intrinsic Matrix
K = inv(K);
K = K/K(3,3);