% DecompPMat - decompose the camera projection matrix P into intrinsic
%              matrix K, rotation matrix R, and translation vector t.
%
% Usage:
%           [K,R,t] = DecompPMat(P)
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
% [2] E. Trucco and A. Verri, "Introductory Techniques for 3-D Computer Vision,"
%     Prentice Hall, pp.134-135, 1998.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% June 2008  - Original version.

function [K,R,t] = DecompPMat(P)

%% Normalized Projection Matrix
scale = sqrt(P(3,1)^2 + P(3,2)^2 + P(3,3)^2);
P = P/scale;

%% Principal Point
K(1,3) = P(1,1:3)*P(3,1:3)';
K(2,3) = P(2,1:3)*P(3,1:3)';

%% Focal Length
K(1,1)   = sqrt(P(1,1:3)*P(1,1:3)' - K(1,3)^2);
K(2,2)   = sqrt(P(2,1:3)*P(2,1:3)' - K(2,3)^2);
K(3,1:3) = [0 0 1];

%% Translation Vector
t(3,1)   = P(3,4);
t(1,1)   = (P(1,4)-K(1,3)*t(3,1))/K(1,1);
t(2,1)   = (P(2,4)-K(2,3)*t(3,1))/K(2,2);

%% Rotation Matrix
R(3,1)   = P(3,1);
R(3,2)   = P(3,2);
R(3,3)   = P(3,3);
R(1,1:3) = (P(1,1:3)-K(1,3)*P(3,1:3))/K(1,1);
R(2,1:3) = (P(2,1:3)-K(2,3)*P(3,1:3))/K(2,2);

%% Orthogonality Enforcement
[U,D,V] = svd(R);
D = eye(3,3);
R = U*D*V';

%% Tz Sign fixing
if (t(3,1) < 0)
    t = -t;
    R = -R;
end