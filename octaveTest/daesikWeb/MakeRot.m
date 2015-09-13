% MakeRot - makes the rotation matrix.
%
% Usage:
%           R = MakeRot(alpha, beta, gamma, mode)
%
% Input:
%           angleZ: angle about Z axis (in degree)
%           angleY: angle about Y axis (in degree)
%           angleX: angle about X axis (in degree)
%           mode : selects the representation of orientation
%                  'Euler_ZYX', 'Fixed_XYZ'
%                  (currently only above two out of 24 angle set is available)
%
% Output:
%           R: 3x3 Rotation Matrix
%
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Jun. 2011 - Original version.


function R = MakeRot(angleZ, angleY, angleX, mode)

angleZ = angleZ/180*pi;
angleY = angleY/180*pi;
angleX = angleX/180*pi;

if ((mode == 'Euler_ZYX') | (mode == 'Fixed_XYZ'))
    Rz = [cos(angleZ) -sin(angleZ) 0
          sin(angleZ)  cos(angleZ) 0
                   0            0  1];

    Ry = [ cos(angleY) 0 sin(angleY)
                    0  1          0
          -sin(angleY) 0 cos(angleY)];

    Rx = [1          0            0 
          0 cos(angleX) -sin(angleX)
          0 sin(angleX)  cos(angleX)];

    R = Rz*Ry*Rx;
end