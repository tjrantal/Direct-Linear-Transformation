% RefinePMat - refines the 3x4 perspective projection matrix 
%              with LM(Levenberg-Marquardt) nonlinear least squares algorithm.
%
% Usage:
%           [P, rperr] = RefinePMat(x, X, P)
%           [P, rperr] = RefinePMat(x, X, P, rerr, iter)
%
% Input:
%           x: 2xn, nx2, 3xn or nx3 image points 
%           X: 3xn, nx3, 4xn or nx4 3D points
%           P: (initial) 3x4 perspective projection matrix
%           rerr: relative error between the last and preceding iteration.
%                 (default: 2^(-52) (It is close to 2.2204e-016.))
%           iter: the number of maximum iteration (default : 30)
%
% Output:
%           P: refined 3x4 perspective projection matrix
%           rperr: average of the re-projection error
%
% cf.:
%           x ~ P*X
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Feb. 2012 - Original version.
% Sep. 2012 - Input arguments are modified.


function [P, rperr] = RefinePMat(varargin)

if (nargin == 3)
    [x, X, P] = varargin{:};
    rerr = 2^(-52);
    iter = 30;
elseif (nargin == 5)
    [x, X, P, rerr, iter] = varargin{:};
end

%% The Number of Points
xSize = size(x);
XSize = size(X);


%% source points
if (xSize(1) == 2)
    x = [x; ones(1,xSize(2))];
elseif (xSize(2) == 2)
    x = [x'; ones(1,xSize(1))];
elseif (xSize(2) == 3)
    x = x';
end


%% destination points
if (XSize(1) == 3)
    X = [X; ones(1,XSize(2))];
elseif (XSize(2) == 3)
    X = [X'; ones(1,XSize(1))];
elseif (XSize(2) == 4)
    X = X';
end


if (length(x) ~= length(X))
    error('The size of x and/or X are wrong.');
end

noPnts = length(x);


%% Evaluate the Jacobian
noParam = 12; % 12 parameters (p11, p12, p13, p14, p21, ..., p33)
J     = zeros(2*noPnts, noParam);  
dist  = zeros(2*noPnts,1);
delta = zeros(noParam,1);
rperr = inf;    % initial error

for n = 1:iter
    % 3x4 perspective projection Matrix
    p = reshape(P',12,1);
    
    p11 = p(1) + delta(1);  p12 = p(2 ) + delta(2 );  p13 = p(3 ) + delta(3 ); p14 = p(4 ) + delta(4 );
    p21 = p(5) + delta(5);  p22 = p(6 ) + delta(6 );  p23 = p(7 ) + delta(7 ); p24 = p(8 ) + delta(8 );
    p31 = p(9) + delta(9);  p32 = p(10) + delta(10);  p33 = p(11) + delta(11); p34 = p(12) + delta(12);
    
    P_lm = [p11 p12 p13 p14;
            p21 p22 p23 p24;
            p31 p32 p33 p34];

        
    % Re-Projection    
    x_rp = P_lm*X;

    
    % Cost Function: Geometric error between the re-projected abd measured points
    dist_lm(1:2:2*noPnts,1) = x_rp(1,:)./x_rp(3,:) - x(1,:);
    dist_lm(2:2:2*noPnts,1) = x_rp(2,:)./x_rp(3,:) - x(2,:);
    
    
    % Re-Projection Error
    rperr_lm = sqrt(dot(dist_lm,dist_lm)/noPnts/2);
    
    
    if (rperr_lm <= rperr)
        if (((n > 1) && sqrt(dot(delta,delta)/dot(p,p)) < rerr))
            P     = P_lm;
            rperr = rperr_lm;
            break;
        end
        
        % Update
        P     = P_lm;
        dist  = dist_lm;
        rperr = rperr_lm;
       
        for (i=1:noPnts)
            %% Jabobian of the 3x4 perspective projection matrix
            df1_dp11 = X(1,i)/x_rp(3,i);
            df1_dp12 = X(2,i)/x_rp(3,i);
            df1_dp13 = X(3,i)/x_rp(3,i);
            df1_dp14 =      1/x_rp(3,i);
            
            df1_dp21 = 0;
            df1_dp22 = 0;
            df1_dp23 = 0;
            df1_dp24 = 0;
            
            df1_dp31 = -(X(1,i)*x_rp(1,i))/x_rp(3,i)^2;
            df1_dp32 = -(X(2,i)*x_rp(1,i))/x_rp(3,i)^2;
            df1_dp33 = -(X(3,i)*x_rp(1,i))/x_rp(3,i)^2;
            df1_dp34 = -(       x_rp(1,i))/x_rp(3,i)^2;

            df2_dp11 = 0;
            df2_dp12 = 0;
            df2_dp13 = 0;
            df2_dp14 = 0;
            
            df2_dp21 = X(1,i)/x_rp(3,i);
            df2_dp22 = X(2,i)/x_rp(3,i);
            df2_dp23 = X(3,i)/x_rp(3,i);
            df2_dp24 =      1/x_rp(3,i);
            
            df2_dp31 = -(X(1,i)*x_rp(2,i))/x_rp(3,i)^2;
            df2_dp32 = -(X(2,i)*x_rp(2,i))/x_rp(3,i)^2;
            df2_dp33 = -(X(3,i)*x_rp(2,i))/x_rp(3,i)^2;
            df2_dp34 = -(       x_rp(2,i))/x_rp(3,i)^2;
            
            
            J(2*i-1, :)  = [df1_dp11 df1_dp12 df1_dp13 df1_dp14 ...
                            df1_dp21 df1_dp22 df1_dp23 df1_dp24 ...
                            df1_dp31 df1_dp32 df1_dp33 df1_dp34];
                        
            J(2*i  , :)  = [df2_dp11 df2_dp12 df2_dp13 df2_dp14 ...
                            df2_dp21 df2_dp22 df2_dp23 df2_dp24 ...
                            df2_dp31 df2_dp32 df2_dp33 df2_dp34];
        end

        
        % Compute the approximated Hessian matrix
        He = J'*J;
        
        if (n == 1)
            lambda = 0.001*trace(He)/noParam;
        else
            lambda = lambda/10;
        end
    else
        lambda = lambda*10;
    end
    
   
    % Apply the damping factor to the Hessian matrix
    He_lm = He + (lambda * eye(noParam, noParam));
    
    
    % Prevent the matrix from being singular
    if (rcond(He_lm) < eps)
        lambda = lambda*10;
        He_lm = He + (lambda * eye(noParam, noParam));
    end

    
    % Compute the updated parameters
    delta = -inv(He_lm)*(J'*dist(:));
end