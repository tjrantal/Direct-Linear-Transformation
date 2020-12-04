%@x are the 16 DLT coefficients and x0 and y0
%Equations taken from http://kwon3d.com/theory/dlt/dlt.html
function y = backProjectOptim(x,eqsToUse)
    global calibObject digitizedCoords iterations eqsToUse imCentre
    
    %Backproject calibObject
    distorted = [];
    imCentre = getPrincipalPoint(x);
    
    opts1=  optimoptions('lsqnonlin','display','off', 'FunctionTolerance', 1e-12,'StepTolerance',1e-12,'OptimalityTolerance', 1.0000e-12);
    tic
    for r = 1:size(calibObject,1)
        backProjected = backproject16(x, calibObject(r,:),x); %These are error-less backprojections
       	distorted(r,:) = lsqnonlin(@adjustError,backProjected,[],[],opts1);
    end
    toc
%     keyboard;
    
    y = [distorted(:)-digitizedCoords(:)]';    
    iterations = iterations+1;
%     keyboard;
end