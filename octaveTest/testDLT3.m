close all;
clear all;
clc;

calibrationDigitized = 0; %=0 if you haven't digitized yet, =1 if you have
dataSaveName = 'digitizedData2.mat';
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
imageNames = {"barrel1.jpg","barrel2.jpg"};

if calibrationDigitized == 1
	load(dataSaveName);	%Load digitization
else
	%Digitize calibration
	for i = 1:length(imageNames)
		[cam(i).image, discard ,discard] = imread(imageNames{i});
		cam(i).digitizedCoordinates = digitizeCalibration(cam(i).image,calibrationFrame);
	end
	save('-mat', dataSaveName, 'cam');%Save the digitization results
end

%Calculate DLT-coefficients
for i = 1:length(cam)
	cam(i).coeffs = getDLTcoeffs(calibrationFrame,cam(i).digitizedCoordinates);
end

%disp('camera position')
for i = 1:length(cam)
	cam(i).position = getCamPosition(cam(i).coeffs);
	cam(i).principalPoint = getPrincipalPoint(cam(i).coeffs);
	cam(i).principalPoint
end

%Backproject calibrationpoints

for i = 1:length(cam)
	figure
	temp = zeros(size(calibrationFrame,1),2);
	for j = 1:size(calibrationFrame,1)
		temp(j,:) = backproject(cam(i).coeffs, calibrationFrame(j,:));
 	end
 	cam(i).backprojectedPoints = temp;
 	imshow(cam(i).image)
 	hold on;
 	plot(cam(i).digitizedCoordinates(:,1),cam(i).digitizedCoordinates(:,2),'r.')
 	plot(cam(i).backprojectedPoints(:,1),cam(i).backprojectedPoints(:,2),'go')
 	
end

