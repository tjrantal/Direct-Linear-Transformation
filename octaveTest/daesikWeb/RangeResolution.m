% RangeResolution - calculates the range resolution between two cameras.
%
% Usage:    
%           [r1, r2, X1, X2] = RangeResolution(P1,P2,u,v,z)
%
% Input:
%           P1 : 3x4 camera projection matrix of the first camera
%           P2 : 3x4 camera projection matrix of the second camera
%           u  : 2D image point coordinate (horizontal direction)
%           v  : 2D image point coordinate (vertical direction)
%           z  : distance from the camera center
%
% Output:
%           r1  : range resolution when u2 of the second camera becomes u2+delta_u2
%           r2  : range resolution when v2 of the second camera becomes v2+delta_v2
%           X1  : 4x1 homogeneous scene point
%           X2  : 4x1 homogeneous scene point
%
% cf.
%           The world coordinates must coincide with the first camera's
%           coordinates.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% June 2008  - Original version.

function [r1, r2, X1, X2] = RangeResolution(P1,P2,u,v,z)

%% 3D Scene Point
X = inv(P1(:,1:3))*[u;v;1]*z;

%% Corresponding 2D image point of the second camera
U2 = P2*[X;1];
scale = U2(3);

%% Delta
% s(u + delta_u) = su + s*delta_u
delta_u = scale;
delta_v = scale;

%% Newly reconstructed 3D scene point
X1 = inv(P2(:,1:3))*(U2-[delta_u; 0; 0] - P2(:,4));
X2 = inv(P2(:,1:3))*(U2-[0; delta_v; 0] - P2(:,4));

%% Range Resolution
r1 = X - X1;
r2 = X - X2;