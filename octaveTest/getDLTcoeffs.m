function coefficients = getDLTcoeffs(calibrationObjectGlobalCoordinates,digitizedCoordinates)
	B = zeros(2*length(calibrationObjectGlobalCoordinates),11);%Matrix for solving DLT-parameters
	C = zeros(2*length(digitizedCoordinates),1); %Digitized calibrationObjectGlobalCoordinates coordinates
	monta = 0;
	for j =1:length(digitizedCoordinates)
		for i =1:2
			monta = monta+1;
			C(monta) = digitizedCoordinates(j,i);
		end
	end
	for i=1:length(calibrationObjectGlobalCoordinates)
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
		B(2*i-1,11)     =-calibrationObjectGlobalCoordinates(i,3)*digitizedCoordinates(i,1);
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
	end
	coefficients = B\C; %Solve the coefficients w/ least squares method
	if 0
		%Compare to manually calculating, and pseudoinverse
		coeffManual = inv(B'*B)*(B'*C);
		coeffPInv = pinv(B)*C;
	
		%Calculate Moore-Pentrose pseudoinverse manually
		%http://people.revoledu.com/kardi/tutorial/LinearAlgebra/SVD.html
		Bb = B*B';
		Bc = B'*B;
		[U, discard] = eig(Bb);
		[V, discard] = eig(Bc);
		D = U'*B*V;
		tol = eps*max(size(D))*max(D(:));
		Dinv = zeros(size(D,1),size(D,2));
		Dinv(find(D > tol)) = 1./D(find(D > tol));
		Dinv = Dinv';
		coeffPInvManual = (V*Dinv*U')*C;
	end
end

