%Downloaded camera calibration functions from http://www.daesik80.com/matlabfns/matlabfns.htm

close all;
clear all;
%clc;

addpath('daesikWeb');	%Daesik's matlab functions
calibrationDigitized = 1; %=0 if you haven't digitized yet, =1 if you have
dataSaveName = 'digitizedGoPro.mat';
%Corners of the rubic cube squares are used as the calibration object
%Origin is back lower corner of white side
%X-axis extends along the white side
%Y-axis along the orange side
%Z-axis is pependicular to the blue side plane
%the frame is a 3 column matrix, column 1 = X, 2 =  Y, 3 = Z
calibrationFrame = [0,0,0;1,0,0;2,0,0;3,0,0; ... %White bottom
					0,0,1;1,0,1;2,0,1;3,0,1; ... %White lower mid
					0,0,2;1,0,2;2,0,2;3,0,2; ... %White upper mid
					0,0,3;1,0,3;2,0,3;3,0,3; ... %White top
					3,1,0;3,2,0;3,3,0; ... %Orange bottom
					3,1,1;3,2,1;3,3,1; ... %Orange lower mid
					3,1,2;3,2,2;3,3,2; ... %Orange upper mid
					3,1,3;3,2,3;3,3,3; ... %Orange top
					0,1,3;0,2,3;0,3,3; ... %Blue back
					1,1,3;1,2,3;1,3,3; ... %Blue back mid
					2,1,3;2,2,3;2,3,3; ... %Blue front mid
										];

imageNames = {"GOPR0093.JPG","GOPR0099.JPG"};
coordNames = {"GoPro0093.xls","GoPro0099.xls"};

if calibrationDigitized == 1
	load(dataSaveName);	%Load digitization
else
	%Digitize calibration
	%for i = 1:length(imageNames)
	%	[cam(i).image, discard ,discard] = imread(imageNames{i});
	%	cam(i).digitizedCoordinates = digitizeCalibration(cam(i).image,calibrationFrame);
	%end
	%READ digitized points from text files. Digitized the images, and created the text files with imageJ 
	for i = 1:length(imageNames)
		[cam(i).image, discard ,discard] = imread(imageNames{i});
		tempCoords = dlmread(coordNames{i},"\t",1,1);
		cam(i).digitizedCoordinates = dlmread(coordNames{i},"\t",1,1);
	end
	save('-mat', dataSaveName, 'cam');%Save the digitization results
end

%Subsample the image to not run out of memory...
figure('position',[10 10 1000 500])
for i = 1:length(cam)
	digitizedCoords = cam(i).digitizedCoordinates./5;	%Scale the coordinates down by 5
	tempImage = cam(i).image;
	%clear cam;
	scaledImage = imresize(tempImage(:,:,1),0.2,'nearest');
	scaledImage(:,:,2) = imresize(tempImage(:,:,2),0.2,'nearest');
	scaledImage(:,:,3) = imresize(tempImage(:,:,3),0.2,'nearest');
	
	subplot(length(cam),2,i*2-1)
	imshow(scaledImage,[]);
	hold on;
	plot(digitizedCoords(:,1),digitizedCoords(:,2),'r.')
	imSize = size(scaledImage);
	[K, R, t, rperr] = CalibTsai(digitizedCoords, calibrationFrame, imSize(2)/2, imSize(1)/2);
	disp(['cam ' num2str(i) ' rperr no correction ' num2str(rperr)]);	
	%Get distortion coefficients
	%[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0], R, t);
	[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0,0], R, t);
	%[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0,0,0,0,0], R, t);
	%[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0,0,0,0,0], R, t);
	disp(['cam ' num2str(i) ' rperr with correction ' num2str(rperr)]);
	%Test calibration
	undistorted = UndistImage(scaledImage, K, d);
	subplot(length(cam),2,i*2)
	imshow(undistorted,[])
	%Plot undistorted digitized points
	hold on;
	% Compute the undistorted normalized point
	wh = size(digitizedCoords,1);
	k1 = d(1);
  k2 = d(2);
  if length(d) == 5
    
    k3 = d(5);
    p1 = d(3); 
    p2 = d(4);
  end
	xx_u = inv(K)*[digitizedCoords(:,1)';digitizedCoords(:,2)';ones(1,wh)];
	% Compute the distorted normalized point
	x_u = xx_u(1, :);
	y_u = xx_u(2, :);
	r = sqrt(x_u.^2 + y_u.^2);
  if length(d) ==5
    radial = (1 + k1*r.^2 + k2*r.^4  + k3*r.^6);
    xx_d(1, :) = radial.*x_u + 2*p1*x_u.*y_u + p2*(r.^2 + 2*x_u.^2);
    xx_d(2, :) = radial.*y_u + p1*(r.^2 + 2*y_u.^2) + 2*p2*x_u.*y_u;
  else
    %radial = (1 + k1*r.^2 + k2*r.^4);
    %Use the inverse to visualise the result
    radial = 1./(1 + k1*r.^2 + k2*r.^4);
    xx_d(1, :) = radial.*x_u;
    xx_d(2, :) = radial.*y_u;
    

  end
	xx_d(3, :) = ones(1, wh);
	% Compute the distorted image point
	xx_dim = K*xx_d;
	plot(xx_dim(1,:),xx_dim(2,:),'g.')
end


