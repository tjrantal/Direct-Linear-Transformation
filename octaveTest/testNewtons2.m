%Testing Newton's method to solve a group of equations
%https://en.wikipedia.org/wiki/Newton%27s_method
%http://math.stackexchange.com/questions/466809/solving-a-set-of-equations-with-newton-raphson
%https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant
%Example taken from https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant

%solve equations
%-2x^2 + 3xy   + 4 sin(y)-6 = 0
% 3x^2 - 2xy^2 + 3 cos(x)+4 = 0

%f(x,y) = [-2x^2 + 3xy   + 4 sin(y)-6;3x^2 - 2xy^2 + 3 cos(x)+4]
%Therefore
%f1(x,y) = -2x^2 + 3xy   + 4 sin(y)-6
%f2(x,y) = 3x^2 - 2xy^2 + 3 cos(x)+4
%Derivatives
%df1(x,y)/dx = -4*xy(1) + 3*xy(2)
%df1(x,y)/dy = 3*xy(1)+4*cos(xy(2))
%df2(x,y)/dx = 6*xy(1) - 2*xy(2)^2 - 3*sin(xy(1))
%df2(x,y)/dy = -4*xy(1)*xy(2)
%Jacobian
%f'(x,y) = [-4*xy(1) + 3*xy(2), 3*xy(1)+4*cos(xy(2));6*xy(1) - 2*xy(2)^2 - 3*sin(xy(1)), -4*xy(1)*xy(2)]

%Newton's method
%x1 = x0 - f(x0) / f '(x0)
%x(k)=x(k−1)−J(x(k−1))^−1*F(x(k−1)). 


function testNewtons2
	%Define the function
	function fXY = f(xy)
		fXY = [-2*xy(1)^2+3*xy(1)*xy(2)+ 4*sin(xy(2))-6 ; ...
			 3*xy(1)^2-2*xy(1)*xy(2)^2+3*cos(xy(1))+4];
	endfunction
	%Define the jacobian of the function
	function jacXY = jac(xy)
		jacXY = [-4*xy(1) + 3*xy(2), 3*xy(1)+4*cos(xy(2)) ; ...
			6*xy(1) - 2*xy(2)^2 - 3*sin(xy(1)), -4*xy(1)*xy(2) ];
	endfunction

	%Start solving from x = 1, y = 2
	x = 1;
	y = 2;
	xy = [x;y]
	fXY = f(xy)
	diffXY = 100;
	iterations = 0;
	while diffXY > 0.001
		xy1 = xy-inv(jac(xy))*f(xy);
		fXY1 = f(xy1);
		diffXY = sqrt(sum((xy1-xy).^2));
		xy = xy1;
		iterations = iterations+1;
	end
	disp(['iterations ' num2str(iterations) ' x, y '])%
	xy
	fXY1
endfunction

