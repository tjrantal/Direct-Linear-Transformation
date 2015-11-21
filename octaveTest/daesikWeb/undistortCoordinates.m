%Function to undistort digitized coordinates using Deasik's coefficients
%@argument coordinates the digitized coordinates Nx2
%@argument K camera matrix
%@argument d Daesik's distortion coefficients
%@return uc undistorted coordinates
%Written by Timo Rantalainen 2015 tjrantal at gmail dot com
function uc = undistortCoordinates(coordinates,K,d)
  	% Compute the undistorted normalized point
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
	xx_u = inv(K)*[coordinates(:,1)';coordinates(:,2)';ones(1,wh)];
	% Compute the distorted normalized point
	x_u = xx_u(1, :);
	y_u = xx_u(2, :);
	r = sqrt(x_u.^2 + y_u.^2);
  if length(d) ==5
    %radial = (1 + k1*r.^2 + k2*r.^4  + k3*r.^6);
	radial = 1./(1 + k1*r.^2 + k2*r.^4  + k3*r.^6);
    %xx_d(1, :) = radial.*x_u + 2*p1*x_u.*y_u + p2*(r.^2 + 2*x_u.^2);
    %xx_d(2, :) = radial.*y_u + p1*(r.^2 + 2*y_u.^2) + 2*p2*x_u.*y_u;
	  xx_d(1, :) = radial.*x_u - 2*p1*x_u.*y_u - p2*(r.^2 + 2*x_u.^2);
    xx_d(2, :) = radial.*y_u - p1*(r.^2 + 2*y_u.^2) - 2*p2*x_u.*y_u;
  end
  if length(d) ==2
    %radial = (1 + k1*r.^2 + k2*r.^4);
    %Use the inverse to visualise the result
    radial = 1./(1 + k1*r.^2 + k2*r.^4);
    xx_d(1, :) = radial.*x_u;
    xx_d(2, :) = radial.*y_u;
  end
  if length(d) ==1
    radial = 1./(1 + k1*r.^2);
    xx_d(1, :) = radial.*x_u;
    xx_d(2, :) = radial.*y_u; 
  end
	xx_d(3, :) = ones(1, wh);
	% Compute the distorted image point
	uc = K*xx_d;
endfunction