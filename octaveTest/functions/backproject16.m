function backprojectedPoint = backproject16(coefficients, globalCoordinates,coeffsIn)
    global eqsToUse coeffs16 backprojected
    coeffs16 = coeffsIn;
    backprojected = backproject(coefficients,globalCoordinates);
    opts1=  optimoptions('lsqnonlin','display','off', 'FunctionTolerance', 1e-12,'StepTolerance',1e-12,'OptimalityTolerance', 1.0000e-12);
    backprojectedPoint = lsqnonlin(@adjustError,backprojected,[],[],opts1);
end
