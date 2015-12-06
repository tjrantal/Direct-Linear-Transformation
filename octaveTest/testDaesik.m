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
notes = struct();
if 0
notes.fh = figure('position',[10 10 1000 500])
end
for i = 1:length(cam)
	digitizedCoords = cam(i).digitizedCoordinates./5;	%Scale the coordinates down by 5
	tempImage = cam(i).image;
	%clear cam;
	scaledImage = imresize(tempImage(:,:,1),0.2,'nearest');
	scaledImage(:,:,2) = imresize(tempImage(:,:,2),0.2,'nearest');
	scaledImage(:,:,3) = imresize(tempImage(:,:,3),0.2,'nearest');
	 if 0
	notes.sp(i*2-1) = subplot(length(cam),2,i*2-1);
	imshow(scaledImage,[]);
	hold on;
	plot(digitizedCoords(:,1),digitizedCoords(:,2),'r.')
  end
  %Calculate DLT and backprojection without lens distortion correction
  notes.ocoeffs(i).coeff = getDLTcoeffs(calibrationFrame,digitizedCoords);
  notes.obp(i).bp = backproject(notes.ocoeffs(i).coeff,[1.5,0,1.5]);
 	%plot(notes.obp(i).bp(1),notes.obp(i).bp(2),'k*')
  %Backproject all calibrationframe coordinates
  for cc = 1:size(calibrationFrame,1)
    bbc = backproject(notes.ocoeffs(i).coeff,calibrationFrame(cc,:));
   	if 0
    plot(bbc(1),bbc(2),'go');
    end
  end
  
  
	imSize = size(scaledImage);
  %imSize = size(cam(i).image);
  %digitizedCoords = cam(i).digitizedCoordinates;
	[K, R, t, rperr] = CalibTsai(digitizedCoords, calibrationFrame, imSize(2)/2, imSize(1)/2);
	disp(['cam ' num2str(i) ' rperr no correction ' num2str(rperr)]);	
	
  if 1
    %Test plotting normalised points...
    noPnts = size(calibrationFrame,1);
    XXc = [R t]*[calibrationFrame'; ones(1,noPnts)];
    xu = XXc(1,:)./XXc(3,:);
    yu = XXc(2,:)./XXc(3,:);
    figure('position',[10 10 1000 500])
    subplot(1,2,1)  
    plot(xu,yu,'ro','linestyle','none');
    %axis('equal');
    set(gca,'ydir','reverse');
    %hold on;
    %plot(xu,yu,'ro','linestyle','none');
    uv = K*[xu;yu;ones(1,noPnts)];
    subplot(1,2,2)
    plot(uv(1,:),uv(2,:),'ro','linestyle','none');
    hold on;
    plot(digitizedCoords(:,1),digitizedCoords(:,2),'ko','linestyle','none');
    axis('equal');
    set(gca,'ydir','reverse');
    
  end
  %Get distortion coefficients
	%[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0], R, t);
	%[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0,0], R, t);
	[K, d, R, t, rperr] = RefineCamParam(digitizedCoords, calibrationFrame, K, [0,0,0,0,0], R, t,eps,30);
	
  %Add lens distortion to projected calibration object coordinates
   noPnts = size(calibrationFrame,1);
    XXc = [R t]*[calibrationFrame'; ones(1,noPnts)];
    xu = XXc(1,:)./XXc(3,:);
    yu = XXc(2,:)./XXc(3,:);  
    r = sqrt(xu.^2 + yu.^2);
        %% Lens Distortion Parameters
        k1 =0;k2=0;k3=0;p1=0;p2=0;
        if length(d) > 0
            k1 = d(1);
            
        end
        if length(d) > 1
            k2 = d(2); 
        end
        if length(d) > 3
             p1 = d(3); p2 = d(4);  
        end
        if length(d) ==5 
           
           
            k3 = d(5);
        end
        xd = xu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1*xu.*yu + p2*(r.^2 + 2*xu.^2);
        yd = yu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*yu.^2) + 2*p2*xu.*yu;
    uv = K*[xd;yd;ones(1,noPnts)];
    plot(uv(1,:),uv(2,:),'g+','linestyle','none');
  
  
  %K
	%d
	%R
	%t
	%keyboard;
  disp(['cam ' num2str(i) ' rperr with correction ' num2str(rperr)]);
	
  %Test calibration
	undistorted = UndistImage(scaledImage, K, d);
	if 0
  notes.sp(i*2) = subplot(length(cam),2,i*2);
  end
  
  figure('position',[10 10 1000 550]);
	imshow(undistorted,[])
  hold on;
  
  % Compute the distorted normalized point
  xx_u = inv(K)*[digitizedCoords(:,1)';digitizedCoords(:,2)';ones(1,noPnts)];
	x_u = xx_u(1, :);
	y_u = xx_u(2, :);
	r = sqrt(x_u.^2 + y_u.^2);
  k1 = 0; k2 =0; k3 = 0;
  if length(d) > 0
    k1 = d(1);
   end
  if length(d) >1
    k2 = d(2);
  end
  if length(d) ==5
    k3 = d(5);
  end
  radial = (1 + k1*r.^2 + k2*r.^4 + k3*r.^6);
  xx_d(1, :) =x_u./radial;
  xx_d(2, :) = y_u./radial;
  uc = K*[xx_d; ones(1,noPnts)];
  plot(uc(1,:),uc(2,:),'ro','linestyle','none');
  keyboard
	%Plot undistorted digitized points
	hold on;
  notes.uc(i).coords = undistortCoordinates(digitizedCoords,K,d);
	plot(notes.uc(i).coords(1,:),notes.uc(i).coords(2,:),'g.')
  notes.coeffs(i).coeff = getDLTcoeffs(calibrationFrame,notes.uc(i).coords');
  notes.bp(i).bp = backproject(notes.coeffs(i).coeff,[1.5,0,1.5]);
  plot(notes.bp(i).bp(1),notes.bp(i).bp(2),'k*')
  %Backproject all calibrationframe coordinates
  for cc = 1:size(calibrationFrame,1)
    bbc = backproject(notes.coeffs(i).coeff,calibrationFrame(cc,:));
   plot(bbc(1),bbc(2),'ro');
  end
  keyboard
  
end



