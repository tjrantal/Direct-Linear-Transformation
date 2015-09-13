% UndistImage - undistorts the image by inverse mapping.
%
% Usage:
%        image = UndistImage(image, K, d)
%
% Input:
%       image: distorted image
%       K: 3x3 intrinsic parameter matrix
%       d: lens distortion vector
%          (1x1          : k1, 
%           2x1 (or 1x2) : k1, k2,
%           4x1 (or 4x1) : k1, k2, p1, p2, 
%           5x1 (or 5x1) : k1, k2, p1, p2, k3,
%           where k1, k2, k3 are the radial distortion parameters and
%           p1, p2 are the tangential distortion parameters.) 
%
% Output:
%       image: undistorted image
%
%
% Kim, Daesik, Ph. D
% E-mail  : daesik80[at]gmail.com
% Homepage: http://www.daesik80.com
%
% Apr. 2015 - Original version.



function image = UndistImage(image, K, d)

% Get the size of the image
[height, width, channel] = size(image);
wh = width*height;


% Get the lens distortion parameters
if (length(d) == 1)
    k1 = d(1);  k2 = 0;     k3 = 0;
    p1 = 0;     p2 = 0;
elseif (length(d) == 2)
    k1 = d(1);  k2 = d(2);  k3 = 0;
    p1 = 0;     p2 = 0;
elseif (length(d) == 4)
    k1 = d(1);  k2 = d(2);  k3 = 0;
    p1 = d(3);  p2 = d(4);
elseif (length(d) == 5)
    k1 = d(1);  k2 = d(2);  k3 = d(5);
    p1 = d(3);  p2 = d(4);
end


% Generate the undistorted image point
[x, y] = meshgrid([1:width], [1:height]);
xx_uim(1, :) = reshape(x', [1, wh]);
xx_uim(2, :) = reshape(y', [1, wh]);
xx_uim(3, :) = ones(1, wh);


% Compute the undistorted normalized point
xx_u = inv(K)*xx_uim;


% Compute the distorted normalized point
x_u = xx_u(1, :);
y_u = xx_u(2, :);
r = sqrt(x_u.^2 + y_u.^2);
radial = (1 + k1*r.^2 + k2*r.^4  + k3*r.^6);
xx_d(1, :) = radial.*x_u + 2*p1*x_u.*y_u + p2*(r.^2 + 2*x_u.^2);
xx_d(2, :) = radial.*y_u + p1*(r.^2 + 2*y_u.^2) + 2*p2*x_u.*y_u;
xx_d(3, :) = ones(1, wh);
 

% Compute the distorted image point
xx_dim = K*xx_d;


% Preparation for bilinear interpolation
xf = floor(xx_dim(1,:));
yf = floor(xx_dim(2,:));
xc = ceil(xx_dim(1,:));
yc = ceil(xx_dim(2,:));
dx = xx_dim(1,:) - xf;
dy = xx_dim(2,:) - yf;


% Find the valid point within the image size
valid = find((xf >= 1) & (xf < width ) & (xc >= 1) & (xc < width ) & ...
             (yf >= 1) & (yf < height) & (yc >= 1) & (yc < height));

         
% Get the index number corresponding to image coordinates
yfxf = yf + xf*height - height;
yfxc = yf + xc*height - height;
ycxf = yc + xf*height - height;
ycxc = yc + xc*height - height;


% Compute the interpolated intensity
image = double(image);

for (n = 1:channel)
    intensity = zeros(1, wh);
    intensity(valid) = ...
        (1-dy(valid)).*(1-dx(valid)).*image(yfxf(valid) + (n-1)*wh) + ...
        (1-dy(valid)).*(  dx(valid)).*image(yfxc(valid) + (n-1)*wh) + ...
        (  dy(valid)).*(1-dx(valid)).*image(ycxf(valid) + (n-1)*wh) + ...
        (  dy(valid)).*(  dx(valid)).*image(ycxc(valid) + (n-1)*wh);
    
    image(:, :, n) = reshape(intensity, [width, height])';
end

image = uint8(image);