%Function to undistort digitized coordinates iteratively using Deasik's coefficients
%Only 5 distortion coefficients implemented at the moment
%Uses Gaus-Newton (or could be just Newton's method, not sure on the terminology)
%@argument coordinates the digitized coordinates Nx2
%@argument K camera matrix
%@argument d Daesik's distortion coefficients
%@return uc undistorted coordinates
%Written by Timo Rantalainen 2015 tjrantal at gmail dot com
function dc = undistortCoordinates(coordinates,K,d)

	%Define the functions for Newton's method
	%Define the function
	%@param uxy, undistorted coordinates
	%@param dxy, distorted coordinates
	function fXY = f(u,v,x,y,k1,k2,p1,p2,k3)
		fXY = [u*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+2*p1*u*v+p2*((u^2+v^2)+2*u^2)-x; ...
			 	 v*(1+k1*(u^2+v^2)+k2*(u^2+v^2)^2+k3*(u^2+v^2)^3)+p1*((u^2+v^2)+2*v^2)+2*p2*u*v-y];
	endfunction
	%Define the jacobian of the function
	function jacXY = jac(u,v,k1,k2,p1,p2,k3)
		jacXY = [ ...
		(u^2 + v^2)^3*k3 + (u^2 + v^2)^2*k2 + (u^2 + v^2)*k1 + 2*(3*(u^2 + v^2)^2*k3*u + 2*(u^2 + v^2)*k2*u + k1*u)*u + 6*p2*u + 2*p1*v + 1, ...
2*(3*(u^2 + v^2)^2*k3*v + 2*(u^2 + v^2)*k2*v + k1*v)*u + 2*p1*u + 2*p2*v; ...
2*p1*u + 2*(3*(u^2 + v^2)^2*k3*u + 2*(u^2 + v^2)*k2*u + k1*u)*v + 2*p2*v, ...
(u^2 + v^2)^3*k3 + (u^2 + v^2)^2*k2 + (u^2 + v^2)*k1 + 2*p2*u + 2*(3*(u^2 + v^2)^2*k3*v + 2*(u^2 + v^2)*k2*v + k1*v)*v + 6*p1*v + 1 ...
		];

	endfunction


  	% Get distortion coefficients from d (the if's are a remnant, not to implement the various lengths)
	wh = size(coordinates,1);
	k1 = d(1);
  if length(d) >=  2
    k2 = d(2);
  end
  if length(d) == 5
    p1 = d(3); 
    p2 = d(4);  
    k3 = d(5);
  end
  % Compute the distorted normalized point
	xx_u = inv(K)*[coordinates(:,1)';coordinates(:,2)';ones(1,wh)];
	
	%calculate undistorted coordinates one at a time
	for i = 1:size(xx_u,2)
		x_u = xx_u(1, i);
		y_u = xx_u(2, i); 
		x_d = x_u;
		y_d = y_u;
		diffXY = 100;
		iterations = 0;
		while diffXY > eps && iterations < 100
			%disp(['it ' num2str(iterations) ' x ' num2str(x_u) ' y ' num2str(y_u)]);
			xy1 = [x_u;y_u]-inv(jac(x_u,y_u,k1,k2,p1,p2,k3))*f(x_u,y_u,x_d,y_d,k1,k2,p1,p2,k3);
			diffXY = sqrt(sum((xy1-[x_u;y_u]).^2));
			x_u = xy1(1);
			y_u = xy1(2);
			iterations = iterations+1;
		end
		xx_d(1:2,i) = [x_u;y_u];
	end
	xx_d(3, :) = ones(1, wh);
	% Compute the distorted image points
	dc = K*xx_d;
endfunction
