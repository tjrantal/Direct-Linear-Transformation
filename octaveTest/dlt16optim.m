%@x are the 16 DLT coefficients and x0 and y0
%Equations taken from http://kwon3d.com/theory/dlt/dlt.html
function y = dlt16optim(x)
    global calibObject digitizedCoords iterations
    %Create Matrix A (2*N x 16) to be multiplied with x' (DLT coefficients), and vector B (2*N x 1),
    A = zeros(2*size(calibObject,1),16);
    B = zeros(2*size(calibObject,1),1);
%     keyboard;
    for r = 1:size(calibObject,1)
        R = x(9)*calibObject(r,1)+x(10)*calibObject(r,2)+x(11)*calibObject(r,3);
        c = digitizedCoords(r,1)-x(17); %Coordinate with respect to the image centre
        n = digitizedCoords(r,2)-x(18); %Coordinate with respect to the image centre
        rad = sqrt(c*c+n*n);
        A(2*r-1,1)		=calibObject(r,1);
		A(2*r-1,2)		=calibObject(r,2);
		A(2*r-1,3)		=calibObject(r,3);
		A(2*r-1,4)		=1;
		A(2*r-1,5)		=0;
		A(2*r-1,6)		=0;
		A(2*r-1,7)		=0;
		A(2*r-1,8)		=0;
		A(2*r-1,9)		=-calibObject(r,1)*digitizedCoords(r,1);
		A(2*r-1,10)		=-calibObject(r,2)*digitizedCoords(r,1);
		A(2*r-1,11)     =-calibObject(r,3)*digitizedCoords(r,1);
        
        A(2*r-1,12)		=c*rad^2*R;
		A(2*r-1,13)		=c*rad^4*R;
		A(2*r-1,14)		=c*rad^6*R;
		A(2*r-1,15)		=rad^2+2*c^2*R;
		A(2*r-1,16)     =c*n*R;
        
        
		A(2*r,1)		=0;
		A(2*r,2)		=0;
		A(2*r,3)		=0;
		A(2*r,4)		=0;
		A(2*r,5)		=calibObject(r,1);
		A(2*r,6)		=calibObject(r,2);
		A(2*r,7)		=calibObject(r,3);
		A(2*r,8)		=1;
		A(2*r,9)		=-calibObject(r,1)*digitizedCoords(r,2);
		A(2*r,10)		=-calibObject(r,2)*digitizedCoords(r,2);
		A(2*r,11)		=-calibObject(r,3)*digitizedCoords(r,2);
        
        A(2*r-1,12)		=n*rad^2*R;
		A(2*r-1,13)		=n*rad^4*R;
		A(2*r-1,14)		=n*rad^6*R;
		A(2*r-1,15)		=n*c*R;
		A(2*r-1,16)     =rad^2+2*n^2*R;

        B(2*r-1) = digitizedCoords(r,1);
        B(2*r) = digitizedCoords(r,2);
    end
    
    y = A*x(1:16)'-B;    
    iterations = iterations+1;
%     keyboard;
end