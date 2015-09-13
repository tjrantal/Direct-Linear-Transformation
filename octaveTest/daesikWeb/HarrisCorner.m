% HarrisCorner - detects corners by Harris and Stephens's Method.
%
% Usage:
%           res = HarrisCorner(img, sigma)
%
% Input:
%           img : Gray Level Image
%           sigma : parameter for cornerness (default: 0.04)
%
% Output:
%           res : Corner Response (cornerness)
%
% This code follows the algorithm given by
% [1] C. Harris and M. Stephens. A Combined Corner and Edge Detector. 
%     Proc. Alvey Vision Conf., Univ. Manchester, pp. 147-151, 1988. 
%
% Kim, Daesik
% Intelligent Systems Research Center
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.daesik80.com
%
% Nov 2008  - Original version.

function res = HarrisCorner(img, sigma)

img = double(img);

gx = fspecial('sobel');
gx = gx/4;
gy = gx';

jx = conv2(img, gx, 'same');
jy = conv2(img, gy, 'same');

g = fspecial('gaussian');

jxx = conv2(jx.*jx, g, 'same');
jyy = conv2(jy.*jy, g, 'same');
jxy = conv2(jx.*jy, g, 'same');

res = (jxx.*jyy - jxy.^2) - sigma*(jxx + jyy).^2;