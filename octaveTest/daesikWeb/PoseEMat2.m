% PoseEMat2 - estimate the pose from essential matrix.
%
% Usage:    
%           [R1, R2, t1, t2] = PoseEMat2(E)
%
% Input:
%           E   : essential matrix
%
% Output:
%           R1  : 3x3 rotation matrix 1
%           R2  : 3x3 rotation matrix 2
%           t1   : 3x1 translation vector 1
%           t2   : 3x1 translation vector 2
%
% Note:
%           An essential matrix is defined as 'E = RS' in the original paper,
%           but the essential matrix is implmented as 'E = SR' in this algorithm.
%
% This code follows the algorithm given by 
% [1] H. Longuet-Higgins, "A Computer Algorithm for Reconstructing a Scene
%     from Two Projections," Nature, vol.293, no.10, pp.133-135, 1981.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% June 2008  - Original version.

function [R1, R2, t1, t2] = PoseEMat2(E)

scale = sqrt(trace(E*E')/2);
E = E/scale;
EE = E*E';

t1(1,1) = sqrt(1-EE(1,1));
t1(2,1) = sqrt(1-EE(2,2));
t1(3,1) = sqrt(1-EE(3,3));

t2(1,1) = -sqrt(1-EE(1,1));
t2(2,1) = -sqrt(1-EE(2,2));
t2(3,1) = -sqrt(1-EE(3,3));

W(:,1) = cross(E(:,1),t1);
W(:,2) = cross(E(:,2),t1);
W(:,3) = cross(E(:,3),t1);

R1(:,1) = W(:,1) + cross(W(:,2),W(:,3));
R1(:,2) = W(:,2) + cross(W(:,3),W(:,1));
R1(:,3) = W(:,3) + cross(W(:,1),W(:,2));

W(:,1) = cross(-E(:,1),t1);
W(:,2) = cross(-E(:,2),t1);
W(:,3) = cross(-E(:,3),t1);

R2(:,1) = W(:,1) + cross(W(:,2),W(:,3));
R2(:,2) = W(:,2) + cross(W(:,3),W(:,1));
R2(:,3) = W(:,3) + cross(W(:,1),W(:,2));