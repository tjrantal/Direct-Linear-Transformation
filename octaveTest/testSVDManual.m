	close all;
	clear all;
	clc;

	B = [1,0,0,0,2;0,0,3,0,0;0,0,0,0,0;0,4,0,0,0]';
	Bb = B*B';
	Bc = B'*B;
	[U, discard] = eig(Bb);
	[V, discard] = eig(Bc);
	D = U'*B*V;
	%B2 = U*D*V';
	tol = eps*max(size(D))*max(D(:));
	Dinv = zeros(size(D,1),size(D,2));
	Dinv(find(D > tol)) = 1./D(find(D > tol));
	Dinv = Dinv';
	PInvOct = pinv(B)	
	PInvManual = (V*Dinv*U')

	A = [1 2 0; 2 5 -1; 4 10 -1];
	svd(A)
	B=A;
	Bb = B*B';
	Bc = B'*B;
	[U, discard] = eig(Bb);
	[V, discard] = eig(Bc);
	D = U'*B*V
	[u,s,v] =svd(A)
