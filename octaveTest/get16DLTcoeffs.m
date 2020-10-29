%Calculate camera calibration iteratively using levenbeg-marquardt
%optimisation. Matlab defaults to numerical estimation of the Jacobian so
%don't have to figure out the partial derivatives yourself (I think?).
function coefficients = get16DLTcoeffs(calibrationObjectGlobalCoordinates,digitizedCoordinates,dtl11coeffs,imSize)
    global calibObject digitizedCoords iterations;
    calibObject = calibrationObjectGlobalCoordinates;
    digitizedCoords = digitizedCoordinates;
    iterations = 0;
%     keyboard;
%     opts1=  optimset('display','off', 'FunctionTolerance', 1e-12);
	opts1=  optimoptions('lsqnonlin','display','off','FunctionTolerance', 1e-12,'StepTolerance',1e-12,'OptimalityTolerance', 1.0000e-12);
	
    coefficients = lsqnonlin(@dlt16optim,[dtl11coeffs' zeros(1,5) imSize(2) imSize(1)],[],[],opts1);
    disp(sprintf('Function calls %d',iterations));
end

