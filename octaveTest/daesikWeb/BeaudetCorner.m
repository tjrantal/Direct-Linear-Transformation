% BeaudetCorner - detects corners by P.R. Beaudet's Method.
%
% Usage:
%           gc = BeaudetCorner(img)
%
% Input:
%           img : Gray Level Image
%
% Output:
%           gc : Gaussian Curvature (cornerness)
%
% This code follows the algorithm given by
% [1] P.R. Beaudet, "Rotationally Invariant Image Operators," Int¡¯l Joint
%     Conf. Pattern Recognition, pp. 579-583, 1978.
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Nov. 2008  - Original version.


function gc = BeaudetCorner(img)

%% Sobel Mask
gy = fspecial('sobel');
gy = gy/4;
gx = gy';

%% First Derivative
jx = conv2(img, gx, 'same');
jy = conv2(img, gy, 'same');

%% Second Derivative
jxx = conv2(jx, gx, 'same');
jyy = conv2(jy, gy, 'same');
jxy = conv2(jx, gy, 'same');

%% Gaussian Smoothing
g = fspecial('gaussian', [5,5], 0.5);

jx = conv2(jx, g, 'same');
jy = conv2(jy, g, 'same');

jxx = conv2(jxx, g, 'same');
jyy = conv2(jyy, g, 'same');
jxy = conv2(jxy, g, 'same');

%% Gaussian Curvature
gc = (jxx.*jyy - jxy.^2)/255;