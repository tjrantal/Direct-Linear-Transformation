	for j = 1:size(calibrationFrame,1)
		bp = backproject(cam(i).coeffs, calibrationFrame(j,:))
 		cam(i).backprojectedPoints(j,:) = bp;
 	end
