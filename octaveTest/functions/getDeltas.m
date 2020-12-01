function  deltas = getDeltas(digitizedCoordinates,coeffs)
    global eqsToUse
    deltas = zeros(size(digitizedCoordinates,1),size(digitizedCoordinates,2));
    for r = 1:size(digitizedCoordinates,1)
        c = digitizedCoordinates(r,1)-coeffs(end-1); %Coordinate with respect to the image centre
        n = digitizedCoordinates(r,2)-coeffs(end); %Coordinate with respect to the image centre
        rad = sqrt(c*c+n*n);

        %16 coeff version
%         deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*(r^2+2*c^2)+coeffs(16)*c*n;
%        deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*n*c+coeffs(16)*(r^2+2*n^2);
        
        if eqsToUse == 12
           deltas(r,1) = c*coeffs(12)*rad^2;
           deltas(r,2) = n*coeffs(12)*rad^2;
        end
        if eqsToUse == 13
           deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4);
           deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4);
        end
        if eqsToUse == 14
           deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6);
           deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6);
        end
        if eqsToUse == 15
           deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*(r^2+2*c^2);
           deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*n*c;
        end
        if eqsToUse == 16
           deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*(r^2+2*c^2)+coeffs(16)*c*n;
           deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*n*c+coeffs(16)*(r^2+2*n^2);
        end
        
       
    end
end