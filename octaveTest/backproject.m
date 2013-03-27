function backprojectedPoint = backproject(coefficients, globalCoordinates)
	u = (coefficients(1)*globalCoordinates(1)+coefficients(2)*globalCoordinates(2)+coefficients(3)*globalCoordinates(3)+coefficients(4))/(coefficients(9)*globalCoordinates(1)+coefficients(10)*globalCoordinates(2)+coefficients(11)*globalCoordinates(3)+1);
	v = (coefficients(5)*globalCoordinates(1)+coefficients(6)*globalCoordinates(2)+coefficients(7)*globalCoordinates(3)+coefficients(8))/(coefficients(9)*globalCoordinates(1)+coefficients(10)*globalCoordinates(2)+coefficients(11)*globalCoordinates(3)+1);
	backprojectedPoint = [u,v];
end
