function backprojectedPoint = backproject16(coefficients, globalCoordinates,coeffs16)
    global eqsToUse
    function y = adjustError(x)
        deltas = getDeltas(x,coeffs16);
        y = backprojected-(x-deltas);
%         y = backprojected-(x+deltas);
    end
    backprojected = backproject(coefficients,globalCoordinates);
    opts1=  optimoptions('lsqnonlin','display','off', 'FunctionTolerance', 1e-12,'StepTolerance',1e-12,'OptimalityTolerance', 1.0000e-12);
    backprojectedPoint = lsqnonlin(@adjustError,backprojected,[],[],opts1);

end
