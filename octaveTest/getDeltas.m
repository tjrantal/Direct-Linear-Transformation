function  deltas = getDeltas(digitizedCoordinates,coeffs)
    deltas = zeros(size(digitizedCoordinates,1),size(digitizedCoordinates,2));
    for r = 1:size(digitizedCoordinates,1)
        c = digitizedCoordinates(r,1)-coeffs(17); %Coordinate with respect to the image centre
        n = digitizedCoordinates(r,2)-coeffs(18); %Coordinate with respect to the image centre
        rad = sqrt(c*c+n*n);
       deltas(r,1) = c*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*(r^2+2*c^2)+coeffs(16)*c*n;
       deltas(r,2) = n*(coeffs(12)*rad^2+coeffs(13)*rad^4+coeffs(14)*rad^6)+coeffs(15)*n*c+coeffs(16)*(r^2+2*n^2);
    end
end