close all;
clear all;
clc;

addpath('functions');

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

imageNames = {'sampleFigs/GOPR0093.JPG','sampleFigs/GOPR0099.JPG'};

ah = [];
fh = [];
if calibrationDigitized == 1
	load(dataSaveName);	%Load digitization
    fh = figure;
    cnt = 0;
    for i = 1:length(cam)
        subplot(1,2,i)
        imshow(cam(i).image);
        hold on;
        plot(cam(i).digitizedCoordinates(:,1),cam(i).digitizedCoordinates(:,2),'r*','linestyle','none');
        cnt = cnt+1;
        ah(cnt) = gca();
    end
    
else
	%Digitize calibration
    for i = 1:length(imageNames)
		cam(i).image = imread(imageNames{i});
        if isfield(cam(i),'digitizedCoordinates')
            cam(i).digitizedCoordinates = digitizeCalibration(cam(i).image,calibrationFrame,cam(i).digitizedCoordinates);
        else
            cam(i).digitizedCoordinates = digitizeCalibration(cam(i).image,calibrationFrame);
        end
	end
	save(dataSaveName, 'cam');%Save the digitization results
    for i = 1:length(imageNames)
        writetable(array2table(cam(i).digitizedCoordinates,'VariableNames',{'X','Y'}),sprintf('%s_coords.csv',imageNames{i}(1:end-4)));
    end
end

%sampleIndices = datasample(1:size(calibrationFrame,1),11); %Select random points to use for calib.

sampleIndices = 1:size(calibrationFrame,1); %Use all points for calib


%Calculate DLT-coefficients
for i = 1:length(cam)
		%11 coeffs
	cam(i).coeffs11 = getDLTcoeffs(calibrationFrame(sampleIndices,:),cam(i).digitizedCoordinates(sampleIndices,:));   %Use regular DLT as the initial guess for optimisation
    %16 coeff
	cam(i).coeffs = get16DLTcoeffs(calibrationFrame(sampleIndices,:),cam(i).digitizedCoordinates(sampleIndices,:),cam(i).coeffs11,getPrincipalPoint(cam(i).coeffs11));    %Optimisation for camera distortion
    %Calculate error-corrected coordinates
    cam(i).deltas = getDeltas(cam(i).digitizedCoordinates,cam(i).coeffs);
    cam(i).corrected = cam(i).digitizedCoordinates-cam(i).deltas;
    cam(i).coeffsCorr = getDLTcoeffs(calibrationFrame(sampleIndices,:),cam(i).corrected(sampleIndices,:));   %Use regular DLT after distortion correction
    
    %Visualise camera calibration
    rectified = undistort(cam(i).image,cam(i).coeffs);
    
    %Backproject calibration object!
    set(0, 'CurrentFigure', fh)
    set(fh,'currentaxes',ah(i));
    for r = 1:size(calibrationFrame,1)
        bProjected11 = backproject(cam(i).coeffs11,calibrationFrame(r,:));
        plot(bProjected11(1),bProjected11(2),'m+','linewidth',1);
        bProjected = backproject16(cam(i).coeffsCorr,calibrationFrame(r,:),cam(i).coeffs);
        bProj = backproject(cam(i).coeffsCorr,calibrationFrame(r,:));
        deltaBProj = getDeltas(bProjected,cam(i).coeffs);
        plot(bProjected(1),bProjected(2),'c*','linewidth',1);
        disp(sprintf('cam %d marker %d diff x %.03f diff y %.03f',i,r, bProj(1)-(bProjected(1)-deltaBProj(1)),bProj(2)-(bProjected(2)-deltaBProj(2))));
%         keyboard;
    end
    
    
    
    
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
plot3(calibrationFrame(:,1),calibrationFrame(:,2),calibrationFrame(:,3),'ro','linewidth',1)
hold on


% plot3(globalCoordinatesBasedOnDigitisation(:,1),globalCoordinatesBasedOnDigitisation(:,2),globalCoordinatesBasedOnDigitisation(:,3),'co','linewidth',5);
plot3(global11(:,1),global11(:,2),global11(:,3),'m*','linewidth',1);
plot3(globalCorr(:,1),globalCorr(:,2),globalCorr(:,3),'g+','linewidth',1);

disp(sprintf("Variance dlt 11 %.2f object units",mean(sqrt(sum((global11-calibrationFrame).^2,2)))));
disp(sprintf("Variance dlt 16 %.2f object units",mean(sqrt(sum((globalCorr-calibrationFrame).^2,2)))));


% %Cameras
% colorNames = {'r','b','k','g','y','m','c'};
% for i =1:length(cam)
% 	%plot3(cam(i).position(1),cam(i).position(2),cam(i).position(3),colorNames{i},'linewidth',5)
% 	ua = 2*(-cam(i).position/norm(cam(i).position));
% 	quiver3(cam(i).position(1),cam(i).position(2),cam(i).position(3),ua(1),ua(2),ua(3),colorNames{i},'linewidth',5)
% end

