% CrossRatioPlane - computes the cross-ratio given five coplanar points.
%
% Usage:   
%           [c1, c2] = CrossRatioPlane(x)
%
% Input:
%          x  : nx3 homogeneous points set
%         
% Output:
%          c1 : first cross-ratio
%          c2 : second cross-ratio
%
% This code follows the algorithm given by
% [1] E. Trucco and A. Verri, "Introductory Techniques for 3-D
%     Computer Vision," pp. 258, 1998.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% July 2008  - Original version.

function [c1, c2] = CrossRatioPlane(x)

m431(1,:) = [x(4,1),x(4,2),1];
m431(2,:) = [x(3,1),x(3,2),1];
m431(3,:) = [x(1,1),x(1,2),1];

m432(1,:) = [x(4,1),x(4,2),1];
m432(2,:) = [x(3,1),x(3,2),1];
m432(3,:) = [x(2,1),x(2,2),1];

m521(1,:) = [x(5,1),x(5,2),1];
m521(2,:) = [x(2,1),x(2,2),1];
m521(3,:) = [x(1,1),x(1,2),1];

m421(1,:) = [x(4,1),x(4,2),1];
m421(2,:) = [x(2,1),x(2,2),1];
m421(3,:) = [x(1,1),x(1,2),1];

m531(1,:) = [x(5,1),x(5,2),1];
m531(2,:) = [x(3,1),x(3,2),1];
m531(3,:) = [x(1,1),x(1,2),1];

m532(1,:) = [x(5,1),x(5,2),1];
m532(2,:) = [x(3,1),x(3,2),1];
m532(3,:) = [x(2,1),x(2,2),1];

%% First Cross-Ratio
c1 = det(m431)*det(m521)/det(m421)/det(m531);

%% Second Cross-Ratio
c2 = det(m421)*det(m532)/det(m432)/det(m521);