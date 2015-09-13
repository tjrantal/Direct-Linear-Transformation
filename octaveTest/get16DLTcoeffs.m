function coefficients = get16DLTcoeffs(calibrationObjectGlobalCoordinates,digitizedCoordinates,coefficients,cam)
	C = zeros(2*length(digitizedCoordinates),1); %Digitized calibrationObjectGlobalCoordinates coordinates
	monta = 0;
	for j =1:length(digitizedCoordinates)
		for i =1:2
			monta = monta+1;
			C(monta) = digitizedCoordinates(j,i);
		end
	end

	%ITERATE until a converged solution is found
	if 0
		testOutput = fopen(['testing16Coeffs' num2str(cam) '.tab'],'w');
		for p = 1:length(coefficients)
			if p == length(coefficients)
				fprintf(testOutput,"%f\n",coefficients(p));
			else
				fprintf(testOutput,"%f\t",coefficients(p));
			end
		end
	else 
		for p = 1:length(coefficients)
			if p == length(coefficients)
				printf("%.2f\n",coefficients(p));
			else
				printf("%.2f ",coefficients(p));
			end
		end
	end
	for iteration = 1:50	
		B = zeros(2*length(calibrationObjectGlobalCoordinates),16);%Matrix for solving DLT-parameters
		principalPoint = getPrincipalPoint(coefficients);
		disp(['cam ' num2str(cam) ' iter ' num2str(iteration) ' P1 ' num2str(principalPoint(1)) ' P2 ' num2str(principalPoint(2))]);  
		for i=1:length(calibrationObjectGlobalCoordinates)
			R = coefficients(9)*calibrationObjectGlobalCoordinates(i,1)+coefficients(10)*calibrationObjectGlobalCoordinates(i,2)+coefficients(11)*calibrationObjectGlobalCoordinates(i,3)+1;
			c = digitizedCoordinates(i,1)-principalPoint(1);
			n = digitizedCoordinates(i,2)-principalPoint(2);
			r = sqrt(c^2+n^2);
			B(2*i-1,1)		=calibrationObjectGlobalCoordinates(i,1);
			B(2*i-1,2)		=calibrationObjectGlobalCoordinates(i,2);
			B(2*i-1,3)		=calibrationObjectGlobalCoordinates(i,3);
			B(2*i-1,4)		=1;
			B(2*i-1,5)		=0;
			B(2*i-1,6)		=0;
			B(2*i-1,7)		=0;
			B(2*i-1,8)		=0;
			B(2*i-1,9)		=-calibrationObjectGlobalCoordinates(i,1)*digitizedCoordinates(i,1);
			B(2*i-1,10)		=-calibrationObjectGlobalCoordinates(i,2)*digitizedCoordinates(i,1);
			B(2*i-1,11)		=-calibrationObjectGlobalCoordinates(i,3)*digitizedCoordinates(i,1);
			B(2*i-1,12)		=c*r^2*R;
			B(2*i-1,13)		=c*r^4*R;
			B(2*i-1,14)		=c*r^6*R;
			B(2*i-1,15)		=(r^2+2*c^2)*R;
			B(2*i-1,16)		=2*c*n*R;
			B(2*i,1)		=0;
			B(2*i,2)		=0;
			B(2*i,3)		=0;
			B(2*i,4)		=0;
			B(2*i,5)		=calibrationObjectGlobalCoordinates(i,1);
			B(2*i,6)		=calibrationObjectGlobalCoordinates(i,2);
			B(2*i,7)		=calibrationObjectGlobalCoordinates(i,3);
			B(2*i,8)		=1;
			B(2*i,9)		=-calibrationObjectGlobalCoordinates(i,1)*digitizedCoordinates(i,2);
			B(2*i,10)		=-calibrationObjectGlobalCoordinates(i,2)*digitizedCoordinates(i,2);
			B(2*i,11)		=-calibrationObjectGlobalCoordinates(i,3)*digitizedCoordinates(i,2);
			B(2*i,12)		=n*r^2*R;
			B(2*i,13)		=n*r^4*R;
			B(2*i,14)		=n*r^6*R;
			B(2*i,15)		=2*n*c*R;
			B(2*i,16)		=(r^2+2*n^2)*R;
		end
		coefficients = B\C; %Solve the coefficients w/ least squares method
		if 0
			for p = 1:length(coefficients)
				if p == length(coefficients)
					fprintf(testOutput,"%f\n",coefficients(p));
				else
					fprintf(testOutput,"%f\t",coefficients(p));
				end
			end
		else 
			for p = 1:length(coefficients)
				if p == length(coefficients)
					printf("%.2f\n",coefficients(p));
				else
					printf("%.2f ",coefficients(p));
				end
			end
		end
		keyboard;
	end
	%Converged solution was found
	if 0
		fclose(testOutput);
	end
end

