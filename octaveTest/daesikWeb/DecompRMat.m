% DecompRMat - decompose the rotation matrix into three angles.
%
% Usage:
%           [angz, angy, angx] = DecompRMat(R)
%
% Input:
%           R: 3x3 Rotation Matrix
%
% Output:
%           angz : angle about Z axis (in degree)
%           angy : angle about Y axis (in degree)
%           angx : angle about X axis (in degree)
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Jun. 2011 - Original version.


function [angz, angy, angx] = DecompRMat(R)

angy = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
angz = atan2(R(2,1)/cos(angy), R(1,1)/cos(angy));
angx = atan2(R(3,2)/cos(angy), R(3,3)/cos(angy));

angz = angz/pi*180;
angy = angy/pi*180;
angx = angx/pi*180;