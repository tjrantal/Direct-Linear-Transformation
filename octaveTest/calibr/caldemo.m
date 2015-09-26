clc
load cademo
more off
echo on
%
% ABOUT CAMERA CALIBRATION
%
% In geometric camera calibration the intrinsic and extrinsic
% camera parameters are computed based on measured image coordinates
% of a know calibration target. These parameters may be used to
% correct distorted images or image coordinates.
%
% What we need is the coordinates of the control points and
% corresponding image observations. In CACAL routine you must also
% specify the configuration of the imaging system. The configuration
% data is stored in the file called configc.m and it is acquired
% based on the name of the setup like 'pulnix'. The actual calibration
% data is stored in data matrices. The first three columns of the data 
% matrix contain the 3-D coordinates, the next two columns contain 
% the observations of the control points, and the last three contain
% the normal vectors of the calibration object surface around the control
% points. CACAL supports up to six data matrices, each containing 
% information from distict images.
%
% Press enter
pause
clc
%
% COPLANAR TARGET
%
% The control points of the calibration target can locate in 3-D
% or they can also be coplanar. In case of coplanar target, multiple
% images captured from different positions and angles are required.
% In the following example three images are used.
%
% Let us plot the data set. The first data matrix looks like this.
%
% Press enter
pause
clf
subplot(1,2,1)
plot3(data1(:,1),data1(:,2),data1(:,3),'x')
title('3-D data')
axis([-400 500 -400 500 0 500])
grid  
subplot(1,2,2)
plot(data1(:,4),data1(:,5),'r+')
axis([0 500 0 500])
axis('ij')
title('Image data')
%
% Press enter
pause
clc
% the second data matrix:
clf
subplot(1,2,1)
plot3(data2(:,1),data2(:,2),data2(:,3),'x')
title('3-D data')
axis([-400 500 -400 500 0 500])
grid  
subplot(1,2,2)
plot(data2(:,4),data2(:,5),'r+')
axis([0 500 0 500])
axis('ij')
title('Image data')
%
% Press enter
pause
clc
% the third data matrix:
clf
subplot(1,2,1)
plot3(data3(:,1),data3(:,2),data3(:,3),'x')
title('3-D data')
axis([-400 500 -400 500 0 500])
grid  
subplot(1,2,2)
plot(data3(:,4),data3(:,5),'r+')
axis([0 500 0 500])
axis('ij')
title('Image data')
%
% Press enter
pause
clc
% OK, let's calibrate
[par,pos,iter,res,er,C]=cacal('pulnix',data1,data2,data3);
% Press enter
pause
clc
%
% 3-D TARGET
%
% The control point structure can be also three-dimensional. The
% advantage which is gained is that only one image is required,
% although several images are also supported by CACAL. This demo
% uses the following data matrix:
%
% Press enter
pause
clf
subplot(1,2,1)
plot3(data3d(:,1),data3d(:,2),data3d(:,3),'x')
axis('equal')
title('3-D data')  
grid
subplot(1,2,2)
plot(data3d(:,4),data3d(:,5),'r+')
axis([0 700 20 550])
axis('ij')
title('Image data')
%
% Press enter
pause
clc
% OK, let's calibrate again
[par,pos,iter,res,er,C]=cacal('sonyz',data3d);
% Press enter
pause
previous=std(er);
clc
%
% TREE-STEP PROCEDURE
%
% The control points are often circular, because they are easy to
% make and accurate to measure in subpixel precision from digital
% images. However, measuring the centroid of the ellipse (that is
% a projection of a circle) introduces a systematic error caused by
% perspective projection which is not an affine transformation. 
% Compensating for this error component requires a three-step calibration
% procedure. The Matlab function CACAL performs this procedure. The
% additional information is the radius of the points which is embedded
% in the configuration file and the normal vector of the surface around
% the points. 
%
% Press enter
pause
clc
% Let's see what happens:
[par,pos,iter,res,er,C]=cacal('sony',data3d);
% Press enter
pause
%
% If we compare this result with the previous one, we notice that
% the standard deviation of the residual is slightly reduced:
current=std(er);
[previous current]
%
% Press enter
pause
clc
%
% THE THIRD STEP
%
% The distorted image coordinates can be corrected by using a simple
% inverse model. The parameters of the inverse model are computed
% with the routine called INVMODEL

par2=invmodel('sony',par);

% Press enter
pause
%
% The process of correcting image coordinates is demonstrated by randomly
% generating one thousand image coordinate pairs all over the image area.

r=[rand(1000,1)*768 rand(1000,1)*576];
clf
plot(r(:,1),r(:,2),'rx');
% Press enter
pause
clc
% Now, we may corrupt these coordinates with radial and tangential
% distortion

d=imdist('sony',par,r);

hold
plot(d(:,1),d(:,2),'go');
% Press enter
pause
clc
% The correction is performed with the Matlab function IMCORR

c=imcorr('sony',par2,d);

% the difference between the original and the corrected coordinates
% may be represented by using histograms

clf
subplot(1,2,1)
hist(r(:,1)-c(:,1))
title('Error in x direction')
xlabel('pixels')
subplot(1,2,2)
hist(r(:,2)-c(:,2))
title('Error in y direction')
xlabel('pixels')

% As we can see, the error is smaller than 0.01 pixels.
echo off
