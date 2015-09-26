	%Define the functions
%f(x,y) = [x^2*y;5*x+sin(y)]
	function y = fTest2(x)
		y = zeros(2,1);
		y(1) = x(1)^2*x(2);
	 	y(2) =  5*x(1) + sin(x(2));
	endfunction
