function point = getPrincipalPoint(coefficients)
	u = (coefficients(1)*coefficients(9)+coefficients(2)*coefficients(10)+coefficients(3)*coefficients(11))/(coefficients(9)^2+coefficients(10)^2+coefficients(11)^2);
	v = (coefficients(5)*coefficients(9)+coefficients(6)*coefficients(10)+coefficients(7)*coefficients(11))/(coefficients(9)^2+coefficients(10)^2+coefficients(11)^2);
	point =[u,v];
end
