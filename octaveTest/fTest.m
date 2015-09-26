	%Define the functions
	function y = fTest(x)
		y = zeros(2,1);
		y(1) = -2*x(1)^2 + 3*x(1)*x(2)   + 4*sin(x(2)) - 6;
	 	y(2) =  3*x(1)^2 - 2*x(1)*x(2)^2 + 3*cos(x(1)) + 4;
	endfunction
