%Calculate camera calibration iteratively using levenbeg-marquardt
%optimisation. Matlab defaults to numerical estimation of the Jacobian so
%don't have to figure out the partial derivatives yourself (I think?).
function coefficients = get16DLTcoeffs(calibrationObjectGlobalCoordinates,digitizedCoordinates,dtl11coeffs,imSize)
    global calibObject digitizedCoords iterations eqsToUse imCentre
    calibObject = calibrationObjectGlobalCoordinates;
    digitizedCoords = digitizedCoordinates;
    iterations = 0;
%     keyboard;
%     opts1=  optimset('display','off', 'FunctionTolerance', 1e-12);
	opts1=  optimoptions('lsqnonlin','display','off','FunctionTolerance', 1e-14,'StepTolerance',1e-14,'OptimalityTolerance', 1.0000e-14);
	imCentre = [imSize(2) imSize(1)]./2;
%     coefficients = lsqnonlin(@dlt16optim,[dtl11coeffs' zeros(1,eqsToUse-11) imSize(2) imSize(1)],[],[],opts1);
    initVals = [dtl11coeffs' zeros(1,eqsToUse-11)];
%     coefficients = lsqnonlin(@dlt16optim,initVals,[],[],opts1);

    coefficients = lsqnonlin(@backProjectOptim,initVals,[],[],opts1);
    keyboard;
%     diffs = zeros(size(initVals,1),size(initVals,2));
    diffs = coefficients-initVals;
    compText = '';
    for i = 1:length(diffs)
       compText = [compText ' ' sprintf('%.03f',diffs(i))]; 
    end
    
    disp(sprintf('Function calls %d',iterations));
    disp(sprintf('diffs %s',compText));
    
end

