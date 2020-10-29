function position = getCamPosition(coefficients)
	A = [coefficients(1), coefficients(2), coefficients(3); ...
	coefficients(5), coefficients(6), coefficients(7); ...
	coefficients(9), coefficients(10), coefficients(11)];
	B = [-coefficients(4); -coefficients(8); -1];
	position =A\B;
end
