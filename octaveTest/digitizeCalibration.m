function [digitizedCoordinates]=digitizeCalibration(imageIn,objectCoordinates);
	figure('position',[10,10,800,600]);
	imshow(imageIn);
	hold on;
	digitizedCoordinates = zeros(size(objectCoordinates,1),2);
	for c = 1:size(objectCoordinates,1)
		%coordinateString = ['Digit X' num2str(objectCoordinates(c,1)) ' Y ' num2str(objectCoordinates(c,2)) ' Z ' num2str(objectCoordinates(c,3))];
		%title(coordinateString);
		%disp(coordinateString);
		[x,y,discard]= ginput(1);		
		digitizedCoordinates(c,:) = [x,y];
		%plot(x,y,'ro','linewidth',5)
		disp(['marker ' num2str(c) ' done'])	
	end
	close;
end
%imshow(cam(3).image),c = 17;,coordinateString = ['Digit X' num2str(calibrationFrame(c,1)) ' Y ' num2str(calibrationFrame(c,2)) ' Z ' num2str(calibrationFrame(c,3))];,title(coordinateString);,[x,y,discard]= ginput(1);,cam(3).digitizedCoordinates(c,:)=[x,y];
