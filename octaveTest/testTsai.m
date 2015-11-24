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
for i = 1:length(cam)
	imSize = size(cam(i).image);	
	[K, R, t, rperr] = CalibTsai(cam(i).digitizedCoordinates, calibrationFrame, imSize(2)/2, imSize(1)/2);
	disp(['cam ' num2str(i) ' K']);
	K
	disp('R');
	R
	disp('t');
	t	
	
end



