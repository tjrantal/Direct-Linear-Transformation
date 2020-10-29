function globalCoordinates = getGlobalCoordinates(coefficients,coordinates)
	L1 = zeros(2*size(coefficients,1),3);
	L2 = zeros(2*size(coordinates,1),1);
	for i =1:size(coordinates,1)
		L2(2*i-1)	= coefficients(i,4)- coordinates(i,1);
		L2(2*i)		= coefficients(i,8)- coordinates(i,2);
	end
	
	for i =1:size(coefficients,1)
		L1(2*i-1,1)	=coefficients(i,9)*coordinates(i,1)-coefficients(i,1);
		L1(2*i-1,2)	=coefficients(i,10)*coordinates(i,1)-coefficients(i,2);
		L1(2*i-1,3)	=coefficients(i,11)*coordinates(i,1)-coefficients(i,3);
		L1(2*i,1)	=coefficients(i,9)*coordinates(i,2)-coefficients(i,5);
		L1(2*i,2)	=coefficients(i,10)*coordinates(i,2)-coefficients(i,6);
		L1(2*i,3)	=coefficients(i,11)*coordinates(i,2)-coefficients(i,7);
	end
	globalCoordinates = L1\L2;
end
