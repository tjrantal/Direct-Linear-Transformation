function rectified = undistort(image,coeffs)
	xCoords = 1:5:size(image,2);
	yCoords = 1:5:size(image,1);
%     [XX,YY] = meshgrid(xCoords,yCoords);
	rectified = zeros(size(image,1),size(image,2),size(image,3));
    for r = yCoords
        for c = xCoords
            lookup = [c,r];
            delta =  getDeltas(lookup,coeffs);
            try
                lookup = lookup-delta;
                rowrange = floor(lookup(2)):floor(lookup(2)+1);
                colrange = floor(lookup(1)):floor(lookup(1)+1);
                [XX,YY] = meshgrid(rowrange,colrange);
%                 keyboard;
                rectified(r,c,:) =  [interp2(XX,YY,image(rowrange,colrange,1),lookup(2),lookup(1)) , ...
                                    interp2(XX,YY,image(rowrange,colrange,2),lookup(2),lookup(1)), ...
                                    interp2(XX,YY,image(rowrange,colrange,3),lookup(2),lookup(1))];
            catch
            end
%             disp(sprintf('r %d c %d',r,c));
        end
%         keyboard;
    end
    figure
    imshow(image);
    temp = rectified(yCoords,xCoords,:);
    figure
    imshow(temp);

end