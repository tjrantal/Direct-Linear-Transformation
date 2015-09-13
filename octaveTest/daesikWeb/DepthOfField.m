% DEPTHOFFIELD - computes the depth of field of a camera
%
% Usage:    
%           [dof, z1, z2, z3] = DepthOfField(f, b, Fstop, z)
%
% Input:
%           f       : focal length
%           b       : circle of confusion (acceptable blur diameter)
%           Fstop   : F-stop (F-number) = focal length / the diameter of the lens aperture
%           z       : object distance
%
% Output:
%           dof    : Depth of Field of a camera
%           z1     : near object distance
%           z2     : far object distance
%           z3     : hyperfocal distance
%
% This code follows the algorithm given by 
% [1] "Summary Note 2006-SN-001-EN_Basic Optics, ISRC, SKKU"
%     available at http://www.daesik80.com
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% April 2008  - Original version.

function [dof, z1, z2, z3] = DepthOfField(f, b, Fstop, z)

%% the diameter of the lens aperture
d = f/Fstop;

%% near object distance
z1 = (f*z*(d-b))/(d*f+b*z);

%% far object distance
z2 = (f*z*(d+b))/(d*f-b*z);

%% hyperfocal distance
z3 = d*f/b;

%% depth of field
dof = (2*b*d*f*z*(f+z))/(d^2*f^2 - b^2*z^2);