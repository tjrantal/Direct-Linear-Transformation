function distorted = distort(image,coeffs)
    
    function y = mapDistortion(x)
        temp = getDist(coords);
        disp(sprintf('%.02f %.02f',temp(1),temp(2)));
        y = x-(coords-temp); 
%         keyboard;
    end

    function deltas = getDist(in)
        c = in(1)-imCentre(1);%coeffs(end-1); %Coordinate with respect to the image centre
        n = in(2)-imCentre(2);%coeffs(end); %Coordinate with respect to the image centre
        rad = sqrt(c*c+n*n);
        deltas(1) = c*(coeffs(1)*rad^2+coeffs(2)*rad^4+coeffs(3)*rad^6)+coeffs(4)*(r^2+2*c^2)+coeffs(5)*c*n;
        deltas(2) = n*(coeffs(1)*rad^2+coeffs(2)*rad^4+coeffs(3)*rad^6)+coeffs(4)*n*c+coeffs(5)*(r^2+2*n^2);
        keyboard;
    end

    opts1=  optimoptions('lsqnonlin','display','off', 'FunctionTolerance', 1e-12,'StepTolerance',1e-12,'OptimalityTolerance', 1.0000e-12);
    imCentre = [size(image,1)/2,size(image,2)/2];    

	xCoords = 1:20:size(image,2);
	yCoords = 1:20:size(image,1);
%     [XX,YY] = meshgrid(xCoords,yCoords);
	distorted = uint8(zeros(size(image,1),size(image,2),size(image,3)));
    for r = yCoords
        for col = xCoords
%             lookup = [c,r];
%             delta =  getDeltas(lookup,coeffs);
            coords = [r, col];
            lookup = lsqnonlin(@mapDistortion,[r, col],[],[],opts1);
            if coords(1) ~= lookup(1) && coords(2) ~= lookup(2)
               disp('optimisation worked'); 
            end
            try
%                 lookup = lookup-delta;
%                 lookup = lookup+delta;
                rowrange = floor(lookup(1)):floor(lookup(1)+1);
                colrange = floor(lookup(2)):floor(lookup(2)+1);
                [XX,YY] = meshgrid(rowrange,colrange);
%                 keyboard;
                distorted(r,col,:) =  [interp2(XX,YY,image(rowrange,colrange,1),lookup(2),lookup(1)) , ...
                                    interp2(XX,YY,image(rowrange,colrange,2),lookup(2),lookup(1)), ...
                                    interp2(XX,YY,image(rowrange,colrange,3),lookup(2),lookup(1))];
            catch
            end
%             disp(sprintf('r %d c %d',r,c));
        end
        disp(sprintf('r %d c %d',r,col));
%         keyboard;
    end
    figure
    imshow(image);
    temp = distorted(yCoords,xCoords,:);
    figure
    imshow(temp);
    keyboard;

end