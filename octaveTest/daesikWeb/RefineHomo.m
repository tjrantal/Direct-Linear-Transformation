% RefineHomo - refines the 2D homography with LM(Levenberg-Marquardt) nonlinear least squares algorithm.
%
% Usage:
%           [H, rperr] = RefineHomo(x1, x2, H, rerr, iter)
%
% Input:
%           x1: 2xn, nx2, 3xn or nx3 source points 
%           x2: 2xn, nx2, 3xn or nx3 destination points
%           H: 3x3 (initial) 2D Homography Matrix
%           rerr: relative error between the solutions estimated
%                 in the last and preceding iteration. 
%                 (default: 2^(-52) (It is close to 2.2204e-016.))
%           iter: the number of maximum iteration
%
% Output:
%           H: 3x3 refined 2D Homography Matrix
%           rperr: average of the re-projection error
%
% cf.:
%           x2 ~ H*x1
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% May  2011 - Original version.
% Jan. 2012 - Bug (jacobian computation) is fixed.
%             Input and output argument type is modified.


function varargout = RefineHomo(varargin)

if (nargin == 3)
    [x1, x2, H] = varargin{:};
    rerr = 2^(-52);
    iter = 30;
end


%% The Number of Points
x1Size = size(x1);
x2Size = size(x2);


%% source points
if (x1Size(1) == 2)
    x1 = [x1; ones(1,x1Size(2))];
elseif (x1Size(2) == 2)
    x1 = [x1'; ones(1,x1Size(1))];
elseif (x1Size(2) == 3)
    x1 = x1';
end


%% destination points
if (x2Size(1) == 2)
    x2 = [x2; ones(1,x2Size(2))];
elseif (x2Size(2) == 2)
    x2 = [x2'; ones(1,x2Size(1))];
elseif (x2Size(2) == 3)
    x2 = x2';
end


if (length(x1) ~= length(x2))
    error('The size of x1 and/or x2 are wrong.');
end

noPnts = length(x1);


%% Evaluate the Jacobian
noParam = 9;    %% 9 parameters (h11, h21, h31, h12, ... , h33)
J     = zeros(noPnts, noParam);  
dist  = zeros(noPnts,1);
delta = zeros(noParam,1);
bperr = inf;    % initial error

for n = 1:iter
    % Homography Matrix
    h = reshape(H,9,1);
    
    h11 = h(1) + delta(1);  h12 = h(4) + delta(4);  h13 = h(7) + delta(7);
    h21 = h(2) + delta(2);  h22 = h(5) + delta(5);  h23 = h(8) + delta(8);
    h31 = h(3) + delta(3);  h32 = h(6) + delta(6);  h33 = h(9) + delta(9);
    
    H_lm = [h11 h12 h13;
            h21 h22 h23;
            h31 h32 h33];

        
    % Back-projection    
    x2_bp = H_lm*x1;

    
    % Cost Function: Geometric error between the measured and back-projected points
    dist_lm(1:2:2*noPnts,1) = x2(1,:) - x2_bp(1,:)./x2_bp(3,:);
    dist_lm(2:2:2*noPnts,1) = x2(2,:) - x2_bp(2,:)./x2_bp(3,:);
    
    
    % Back-projection Error
    bperr_lm = sqrt(dot(dist_lm,dist_lm))/noPnts;
    

    if (bperr_lm <= bperr)
        if (((n > 1) && sqrt(dot(delta,delta)/dot(h,h)) < rerr))
            H     = H_lm;
            bperr = bperr_lm;
            break;
        end
        
        % Update
        H     = H_lm;
        dist  = dist_lm;
        bperr = bperr_lm;
       
        for (i=1:noPnts)
            %% Jabobian of the homography matrix
            df1_dh11 = -x1(1,i)/x2_bp(3,i);
            df1_dh21 = 0;
            df1_dh31 = x2_bp(1,i)/x2_bp(3,i)^2*x1(1,i);
            df1_dh12 = -x1(2,i)/x2_bp(3,i);
            df1_dh22 = 0;
            df1_dh32 = x2_bp(1,i)/x2_bp(3,i)^2*x1(2,i);
            df1_dh13 = -1/x2_bp(3,i);
            df1_dh23 = 0;
            df1_dh33 = x2_bp(1,i)/x2_bp(3,i)^2;

            df2_dh11 = 0;
            df2_dh21 = -x1(1,i)/x2_bp(3,i);
            df2_dh31 = x2_bp(2,i)/x2_bp(3,i)^2*x1(1,i);
            df2_dh12 = 0;
            df2_dh22 = -x1(2,i)/x2_bp(3,i);
            df2_dh32 = x2_bp(2,i)/x2_bp(3,i)^2*x1(2,i);
            df2_dh13 = 0;
            df2_dh23 = -1/x2_bp(3,i);
            df2_dh33 = x2_bp(2,i)/x2_bp(3,i)^2;
            
            J(2*i-1, :)  = [df1_dh11 df1_dh21 df1_dh31 ...
                            df1_dh12 df1_dh22 df1_dh32 ...
                            df1_dh13 df1_dh23 df1_dh33];
                        
            J(2*i  , :)  = [df2_dh11 df2_dh21 df2_dh31 ...
                            df2_dh12 df2_dh22 df2_dh32 ...
                            df2_dh13 df2_dh23 df2_dh33];
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


%% Output
if (nargout == 1)
    varargout = {H};
elseif (nargout == 2)
    varargout = {H, bperr};
end