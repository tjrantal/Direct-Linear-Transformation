close all;
clear all;
clc;

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
    figure
    for i = 1:length(cam)
    subplot(1,2,i)
    imshow(cam(i).image);
    hold on;
    plot(cam(i).digitizedCoordinates(:,1),cam(i).digitizedCoordinates(:,2),'r*','linestyle','none');
    end
    
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

%Calculate DLT-coefficients
for i = 1:length(cam)
	%11 coeffs
	cam(i).coeffs11 = getDLTcoeffs(calibrationFrame,cam(i).digitizedCoordinates);   %Use regular DLT as the initial guess for optimisation
    %16 coeff
	cam(i).coeffs = get16DLTcoeffs(calibrationFrame,cam(i).digitizedCoordinates,cam(i).coeffs11,getPrincipalPoint(cam(i).coeffs11));    %Optimisation for camera distortion
    %Calculate error-corrected coordinates
    cam(i).deltas = getDeltas(cam(i).digitizedCoordinates,cam(i).coeffs);
    cam(i).corrected = cam(i).digitizedCoordinates-cam(i).deltas;
    cam(i).coeffsCorr = getDLTcoeffs(calibrationFrame,cam(i).corrected);   %Use regular DLT after distortion correction
    
end
% keyboard;

%Calculate calibrated 3D coordinates
globalCoordinatesBasedOnDigitisation = zeros(size(calibrationFrame,1),size(calibrationFrame,2));
global11 = zeros(size(calibrationFrame,1),size(calibrationFrame,2));
globalCorr = zeros(size(calibrationFrame,1),size(calibrationFrame,2));
for c = 1:size(cam(1).corrected)
    coefficients = [];
    coordinates = [];
    coefficients11 = [];
    coordinates11 = [];
    coefficientsCorr = [];
    coordinatesCorr = [];
    for r = 1:length(cam)
        coefficients(r,:) = cam(r).coeffs;
        coordinates(r,:) = cam(r).corrected(c,:);
        coefficients11(r,:) = cam(r).coeffs11;
        coordinates11(r,:) = cam(r).digitizedCoordinates(c,:);
        coefficientsCorr(r,:) = cam(r).coeffsCorr;
        coordinatesCorr(r,:) = cam(r).corrected(c,:);
    end
    
    globalCoordinatesBasedOnDigitisation(c,:) = getGlobalCoordinates(coefficients,coordinates);
    global11(c,:) = getGlobalCoordinates(coefficients11,coordinates11);
    globalCorr(c,:) = getGlobalCoordinates(coefficientsCorr,coordinatesCorr);
end

% 
% %disp('camera position')
% for i = 1:length(cam)
% 	cam(i).position = getCamPosition(cam(i).coeffs);
% 	cam(i).principalPoint = getPrincipalPoint(cam(i).coeffs);
% end
% 
% %Backproject calibrationpoints
% 
% for i = 1:length(cam)
% 	figure
% 	temp = zeros(size(calibrationFrame,1),2);
% 	for j = 1:size(calibrationFrame,1)
% 		temp(j,:) = backproject(cam(i).coeffs, calibrationFrame(j,:));
%  	end
%  	cam(i).backprojectedPoints = temp;
%  	imshow(cam(i).image)
%  	hold on;
%  	plot(cam(i).digitizedCoordinates(:,1),cam(i).digitizedCoordinates(:,2),'r.')
%  	plot(cam(i).backprojectedPoints(:,1),cam(i).backprojectedPoints(:,2),'go')
%  	
% end
% 
% 
%Plot 3D positions of the markers and the cameras
figure
%Markers
ind =1:16;
plot3(calibrationFrame(ind,1),calibrationFrame(ind,2),calibrationFrame(ind,3),'ro','linewidth',5)
hold on
ind = 17:28;
plot3(calibrationFrame(ind,1),calibrationFrame(ind,2),calibrationFrame(ind,3),'b*','linewidth',5)
ind = 29:37;
plot3(calibrationFrame(ind,1),calibrationFrame(ind,2),calibrationFrame(ind,3),'k+','linewidth',5)

% plot3(globalCoordinatesBasedOnDigitisation(:,1),globalCoordinatesBasedOnDigitisation(:,2),globalCoordinatesBasedOnDigitisation(:,3),'co','linewidth',5);
plot3(global11(:,1),global11(:,2),global11(:,3),'m*','linewidth',5);
plot3(globalCorr(:,1),globalCorr(:,2),globalCorr(:,3),'g+','linewidth',5);

disp(sprintf("Squared error dlt 11 %.2f",sum(sum((global11-calibrationFrame).^2),2)));
disp(sprintf("Squared error dlt 16 %.2f",sum(sum((globalCorr-calibrationFrame).^2),2)));



% %Cameras
% colorNames = {'r','b','k','g','y','m','c'};
% for i =1:length(cam)
% 	%plot3(cam(i).position(1),cam(i).position(2),cam(i).position(3),colorNames{i},'linewidth',5)
% 	ua = 2*(-cam(i).position/norm(cam(i).position));
% 	quiver3(cam(i).position(1),cam(i).position(2),cam(i).position(3),ua(1),ua(2),ua(3),colorNames{i},'linewidth',5)
% end

