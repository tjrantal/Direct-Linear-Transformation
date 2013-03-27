%Digitize point and calculate point location in global coordinates
coefficients = [];
coordinates = [];
figure
for i = 1:length(cam)
	coefficients(i,:) = cam(i).coeffs;
	imshow(cam(i).image);
	[x,y,discard]= ginput(1);		
	coordinates(i,:) = [x,y];
end
close;
globalCoordinates = getGlobalCoordinates(coefficients,coordinates)
