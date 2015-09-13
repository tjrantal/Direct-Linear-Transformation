% PoseQuaternion - estimate the pose with quaternion.
%
% Usage:
%           [R,t] = PoseQuaternion(x1, x2)
%
% Input:
%           x1 : 4xn homogeneous 3D point set 1
%           x2 : 4xn homogeneous 3D point set 2
%
% Output:
%           R : 3x3 rotation matrix
%           t : 3x1 translation vector
%
% cf.:
%           x2 = R*x1 + t
%
% This code follows the algorithm given by
% [1] B.K.P. Horn,¡°Closed-Form Solution of Absolute Orientation Using
%     Unit Quaternions,¡±Journal of the Optical Society of America A, 
%     vol.4, no.4, pp.629-642,s 1987.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% Jul. 2008  - Original version.

function [R,t] = PoseQuaternion(x1, x2)

%% Centroids
c1 = mean(x1(1:3,:)')';
c2 = mean(x2(1:3,:)')';

%% Shift the origin of the points to the centroid
x1(1,:) = x1(1,:) - c1(1);
x1(2,:) = x1(2,:) - c1(2);
x1(3,:) = x1(3,:) - c1(3);

x2(1,:) = x2(1,:) - c2(1);
x2(2,:) = x2(2,:) - c2(2);
x2(3,:) = x2(3,:) - c2(3);

%% Cross Covariance Matrix
H = zeros(3,3);

noPoints = length(x1);
for (n=1:noPoints)
    H_temp = [x1(1,n)*x2(1,n)    x1(1,n)*x2(2,n)    x1(1,n)*x2(3,n);
              x1(2,n)*x2(1,n)    x1(2,n)*x2(2,n)    x1(2,n)*x2(3,n);
              x1(3,n)*x2(1,n)    x1(3,n)*x2(2,n)    x1(3,n)*x2(3,n)];
    H = H + H_temp;
end

%% Symmetric 4x4 Matrix
A = [H(2,3)-H(3,2)   H(3,1)-H(1,3)   H(1,2)-H(2,1)]';
G = [trace(H)   A';
     A          H+H'-trace(H)*eye(3,3)];

%% Eigenvalues and Eigenvectors of G
[V,D] = eig(G);

q0 = V(1,4);
q1 = V(2,4);
q2 = V(3,4);
q3 = V(4,4);

%% Rotation Matrix from unit Quaternion
R(1,1) = q0^2 + q1^2 - q2^2 - q3^2;
R(1,2) = 2*(q1*q2 - q0*q3);
R(1,3) = 2*(q1*q3 + q0*q2);
R(2,1) = 2*(q1*q2 + q0*q3);
R(2,2) = q0^2 - q1^2 + q2^2 - q3^2;
R(2,3) = 2*(q2*q3 - q0*q1);
R(3,1) = 2*(q1*q3 - q0*q2);
R(3,2) = 2*(q2*q3 + q0*q1);
R(3,3) = q0^2 - q1^2 - q2^2 + q3^2;

%% Translation Vector
t = c2(1:3,:) - R*c1(1:3,:);