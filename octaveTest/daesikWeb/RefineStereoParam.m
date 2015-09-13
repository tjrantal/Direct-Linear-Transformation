% RefineStereoParam - refines the stereo camera parameters with LM 
%                     (Levenberg-Marquardt) nonlinear least squares 
%                     algorithm.
%
% Usage:
%        [K1, K2, R, t, rperr] = ...
%           RefineCamParam2(x1, x2, X, K1, R1, t1, K2, R2, t2)
%        [K1, K2, R, t, rperr] = ...
%           RefineCamParam2(x1, x2, X, K1, R1, t1, K2, R2, t2, rerr, iter)
%        [K1, d1, K2, d2, R, t, rperr] = ...
%           RefineCamParam2(x1, x2, X, K1, d1, R1, t1, K2, d2, R2, t2)
%        [K1, d1, K2, d2, R, t, rperr] = ...
%           RefineCamParam2(x1, x2, X, K1, d1, R1, t1, K2, d2, R2, t2, rerr, iter)
%
% Input:
%        x1: 2xn, nx2, 3xn or nx3 image points of the 1st camera
%        x2: 2xn, nx2, 3xn or nx3 image points of the 2nd camera
%        X : 2xn, nx2, 3xn or nx3 planar object points
%            (n is the total number of points)
%        K1: initial 3x3 intrinsic parameters matrix of the 1st camera
%        K2: initial 3x3 intrinsic parameters matrix of the 2nd camera
%            (If the skew is 0, the function does not refine the skew.)
%        R1: initial 3x3m or 3mx3 rotation matrices of the 1st camera
%        R2: initial 3x3m or 3mx3 rotation matrices of the 2nd camera
%        t1: initial 3xm or mx3 translation vectors of the 1st camera
%        t2: initial 3xm or mx3 translation vectors of the 2nd camera
%            (m is the number of planes. If m is 3, it is assumed that 
%             the translation vectors are stacked in the same way as R. For 
%             example, the R is 3x3m, t is regarded as 3xm and vice versa.)
%        d1: initial 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%            lens distortion vector of the 1st camera
%        d2: initial 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%            lens distortion vector of the 2nd camera
%           (1x1        : k1, 
%            2x1(or 1x2): k1, k2,
%            4x1(or 4x1): k1, k2, p1, p2, 
%            5x1(or 5x1): k1, k2, p1, p2, k3,
%            where k1, k2, k3 are the radial distortion parameters and
%            p1, p2 are the tangential distortion parameters.) 
%        rerr: relative error between the last and preceding iteration.
%              (default: 2^(-52) (It is close to 2.2204e-016.))
%        iter: the number of maximum iteration (default : 30)
%
% Output:
%         K1: refined 3x3 intrinsic parameters matrix of the 1st camera
%         K2: refined 3x3 intrinsic parameters matrix of the 2nd camera
%         R : refined 3x3m or 3mx3 rotation matrix between two cameras
%         t : refined 3xm or mx3 translation vectors between two cameras
%         d1: refined 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%            lens distortion vector of the 1st camera
%         d2: refined 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%            lens distortion vector of the 2nd camera
%         rperr: re-projection error
%
%
% Kim, Daesik
% Intelligent Systems Research Institute
% Sungkyunkwan Univ. (SKKU), South Korea
% E-mail  : daesik80@skku.edu
% Homepage: http://www.3DRobotVision.com
%           http://www.daesik80.com
%
% Jun. 2010 - Original version.
% Feb. 2012 - Bug (jacobian computation) is fixed.
%             Input and ouput arguments are modified.


function varargout = RefineStereoParam(varargin)

if (nargin == 9)
    [x1, x2, X, K1, R1, t1, K2, R2, t2] = varargin{:};
    d1 = [];
    d2 = [];
    rerr = 2^(-52);
    iter = 30;
elseif (nargin == 11) && (nargout == 5)
    d1 = [];
    d2 = [];
    [x1, x2, X, K1, R1, t1, K2, R2, t2, rerr, iter] = varargin{:};
elseif (nargin == 11) && (nargout == 7)
    [x1, x2, X, K1, d1, R1, t1, K2, d2, R2, t2] = varargin{:};
    rerr = 2^(-52);
    iter = 30;
elseif (nargin == 13)
    [x1, x2, X, K1, d1, R1, t1, K2, d2, R2, t2, rerr, iter] = varargin{:};
end


%% The Number of Planes
if ((length(R1) == length(R2)) && (length(t1) == length(t2)))
    m = length(R1)/3; % number of plane
else
    error('The size of the R1, R2, t1, and/or t2 are wrong.');
end


%% The Number of Points
XSize  = size(X);
x1Size = size(x1); x2Size = size(x2);
R1Size = size(R1); R2Size  = size(R2);
t1Size = size(t1); t2Size  = size(t2);
d1Size = size(d1); d2Size  = size(d2);

if (length(x1) == length(x2)) && (length(x1) == length(X))
    noPnts = length(x1)/m;  % number of points
else
    error('The size of x1, x2 and/or X are wrong.');
end


%% Image Points
if (x1Size(2) == 2)
    x1 = x1';
elseif (x1Size(1) == 3)
    x1 = x1(1:2,:);
elseif (x1Size(2) == 3)
    x1 = x1(:,1:2)';
end

if (x2Size(2) == 2)
    x2 = x2';
elseif (x2Size(1) == 3)
    x2 = x2(1:2,:);
elseif (x2Size(2) == 3)
    x2 = x2(:,1:2)';
end


%% Object Points
if (XSize(1) == 2)
    XXw = [X; zeros(1,m*noPnts)];
elseif (XSize(2) == 2)
    XXw = [X'; zeros(1,m*noPnts)];
elseif (XSize(1) == 3)
    XXw = [X(1:2,:); zeros(1,m*noPnts)];
elseif (XSize(2) == 3)
    XXw = [X(:,1:2)'; zeros(1,m*noPnts)];
end


%% Rotation Matrix
if (R1Size(2) == 3)
    R1 = R1';

    for (j=1:m)
        R1(:,(j-1)*3+1:j*3) = R1(:,(j-1)*3+1:j*3)';
    end
end

if (R2Size(2) == 3)
    R2 = R2';

    for (j=1:m)
        R2(:,(j-1)*3+1:j*3) = R2(:,(j-1)*3+1:j*3)';
    end
end


%% Translation Vector
if ((t1Size(1) ~= 3) && (t1Size(2) == 3))
    t1 = t1';
elseif ((t1Size(1) == 3) && (t1Size(2) == 3) && (R1Size(2) == 3))
    t1 = t1';
end

if ((t2Size(1) ~= 3) && (t2Size(2) == 3))
    t2 = t2';
elseif ((t2Size(1) == 3) && (t2Size(2) == 3) && (R2Size(2) == 3))
    t2 = t2';
end


%% Distortion Vector
if (d1Size(1) == 1)
    d1 = d1';
end

if (d2Size(1) == 1)
    d2 = d2';
end


%% Convert the 3x3 rotation matrix into the 3x1 vector of Rodigrues representation
for (j=1:m)
    % Rotation of the 1st camera
    w1(:,j) = Rodrigues(R1(:,(j-1)*3+1:j*3));
    
    % Rotation of the 2nd camera
    w2(:,j) = Rodrigues(R2(:,(j-1)*3+1:j*3));
    
    % Rotation between two cameras
    R_temp =  R2(:,3*(j-1)+1:3*j)*R1(:,3*(j-1)+1:3*j)';
    w_temp(:,j) = Rodrigues(R_temp);
    t_temp(:,j) = -R_temp*t1(:,j) + t2(:,j);
end

%% Rotation between two cameras (median value)
if (m == 1)
    w = w_temp;
else
    w = median(w_temp')';
end

%% Translation vector between two cameras (median value)
if (m == 1)
    t = t_temp;
else
    t = median(t_temp')';
end


%% LM(Levenberg-Marquardt) nonlinear least squares algorithm
% The number of parameters: 10 + 6*m + 6 + 10
%  10 (intrinsic of the 1st and 2nd camera) + 
%  6xm (rotation and translation of the 1st camera) + 
%  6 (rotation and translation between two cameras) +
%  10 (distortion of the 1st and 2nd camera)
noParam = 10 + 6*m + 6 + length(d1) + length(d2);
param   = zeros(noParam, 1);
K1_lm   = zeros(3,3);
K2_lm   = zeros(3,3);
w1_lm   = zeros(3,m);
R1_lm   = zeros(3,3*m);
t1_lm   = zeros(3,m);
R_lm    = zeros(3,3);
w_lm    = zeros(3,1);
t_lm    = zeros(3,1);
XXc1    = zeros(3,noPnts*m);
XXc2    = zeros(3,noPnts*m);
J       = zeros(4*m*noPnts, noParam);  
dist_lm = zeros(4*m*noPnts,       1);
delta   = zeros(noParam, 1);
rperr  = inf;    % initial error


for n = 1:iter
    %% Intrinsic Parameters Matrix of the 1st Camera
    K1_lm(1,1) = K1(1,1) + delta(1);  % fx
    K1_lm(2,2) = K1(2,2) + delta(2);  % fy  
    K1_lm(1,3) = K1(1,3) + delta(3);  % u0
    K1_lm(2,3) = K1(2,3) + delta(4);  % v0
    K1_lm(1,2) = K1(1,2) + delta(5);  % s
    K1_lm(3,3) = 1;
    
    K2_lm(1,1) = K2(1,1) + delta(6);  % fx
    K2_lm(2,2) = K2(2,2) + delta(7);  % fy  
    K2_lm(1,3) = K2(1,3) + delta(8);  % u0
    K2_lm(2,3) = K2(2,3) + delta(9);  % v0
    K2_lm(1,2) = K2(1,2) + delta(10); % s
    K2_lm(3,3) = 1;
    
    %% Convert the 3x1 vector of the Rodigrues representation 
    %% into the 3x3 rotation matrix
    w_lm(1) = w(1) + delta(10 + 6*m + 1);  
    w_lm(2) = w(2) + delta(10 + 6*m + 2);  
    w_lm(3) = w(3) + delta(10 + 6*m + 3);
    
    wx = w_lm(1); wy = w_lm(2); wz = w_lm(3);
    theta = sqrt(wx^2 + wy^2 + wz^2);
    
    if (theta < eps)
        R_lm = eye(3);
    else
        wh_sk = [0 -wz wy;wz 0 -wx;-wy wx 0]/theta;

        %% 3x3 Rotation Matrix of the 2nd Camera
        R_lm = eye(3) + sin(theta)*wh_sk + (1-cos(theta))*(wh_sk*wh_sk);
    end

    
    %% Translation Vector of the 2nd Camera
    t_lm(1) = t(1) + delta(10 + 6*m + 4);
    t_lm(2) = t(2) + delta(10 + 6*m + 5);
    t_lm(3) = t(3) + delta(10 + 6*m + 6);
    
        
    for (j=1:m)
        %% Convert the 3x1 vector of the Rodigrues representation 
        %% into the 3x3 rotation matrix
        w1_lm(1, j) = w1(1, j) + delta(10 + 3*(j-1)+1);  
        w1_lm(2, j) = w1(2, j) + delta(10 + 3*(j-1)+2);  
        w1_lm(3, j) = w1(3, j) + delta(10 + 3*(j-1)+3);
        
        wx1 = w1_lm(1, j); wy1 = w1_lm(2, j); wz1 = w1_lm(3, j);
        theta1 = sqrt(wx1^2 + wy1^2 + wz1^2);
        
        if (theta1 < eps)
            R1_lm(:,(j-1)*3+1:j*3) = eye(3);
        else
            wh_sk1 = [0 -wz1 wy1;wz1 0 -wx1;-wy1 wx1 0]/theta1;

            %% 3x3 Rotation Matrix of the 1st Camera
            R1_lm(:,(j-1)*3+1:j*3) = eye(3) + sin(theta1)*wh_sk1 + (1-cos(theta1))*(wh_sk1*wh_sk1);
        end
        
        
        %% Translation Vector of the 1st Camera
        t1_lm(1, j) = t1(1, j) + delta(10 + 3*m + 3*(j-1)+1);
        t1_lm(2, j) = t1(2, j) + delta(10 + 3*m + 3*(j-1)+2);
        t1_lm(3, j) = t1(3, j) + delta(10 + 3*m + 3*(j-1)+3);
        
        
        %% 3D points represented with respect to the camera coordinate frame
        XXc1(:,(j-1)*noPnts+1:j*noPnts) = [R1_lm(:,(j-1)*3+1:j*3), t1_lm(:, j)]*[XXw(:,(j-1)*noPnts+1:j*noPnts); ones(1,noPnts)];
        XXc2(:,(j-1)*noPnts+1:j*noPnts) = [R_lm, t_lm]*[XXc1(:,(j-1)*noPnts+1:j*noPnts); ones(1,noPnts)];
    end
    
    
    %% undistorted normalized points
    xu1 = XXc1(1,:)./XXc1(3,:);
    yu1 = XXc1(2,:)./XXc1(3,:);
    
    xu2 = XXc2(1,:)./XXc2(3,:);
    yu2 = XXc2(2,:)./XXc2(3,:);


    %% Lens Distortion Parameters of the 1st Camera
    if (length(d1) > 0)
        r = sqrt(xu1.^2 + yu1.^2);

        if     (length(d1) == 1) 
            k1 = d1(1) + delta(10+6*m+6+1); k2 = 0 ; 
            p1 = 0; p2 = 0; k3 = 0 ;
            d1_lm = [k1];
        elseif (length(d1) == 2) 
            k1 = d1(1) + delta(10+6*m+6+1); k2 = d1(2) + delta(10+6*m+6+2); 
            p1 = 0; p2 = 0; k3 = 0;
            d1_lm = [k1; k2];
        elseif (length(d1) == 4) 
            k1 = d1(1) + delta(10+6*m+6+1); k2 = d1(2) + delta(10+6*m+6+2); 
            p1 = d1(3) + delta(10+6*m+6+3); p2 = d1(4) + delta(10+6*m+6+4); 
            k3 = 0;
            d1_lm = [k1; k2; p1; p2];
        elseif (length(d1) == 5) 
            k1 = d1(1) + delta(10+6*m+6+1); k2 = d1(2) + delta(10+6*m+6+2); 
            p1 = d1(3) + delta(10+6*m+6+3); p2 = d1(4) + delta(10+6*m+6+4); 
            k3 = d1(5) + delta(10+6*m+6+5);
            d1_lm = [k1; k2; p1; p2; k3];
        end

        xd1 = xu1.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1*xu1.*yu1 + p2*(r.^2 + 2*xu1.^2);
        yd1 = yu1.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*yu1.^2) + 2*p2*xu1.*yu1;
    else
        xd1 = xu1;
        yd1 = yu1;
    end
    
    
    %% Lens Distortion Parameters of the 2nd Camera
    if (length(d2) > 0)
        r = sqrt(xu2.^2 + yu2.^2);

        if     (length(d2) == 1) 
            k1 = d2(1) + delta(10+6*m+6+length(d1)+1); k2 = 0 ; 
            p1 = 0; p2 = 0; k3 = 0 ;
            d2_lm = [k1];
        elseif (length(d2) == 2) 
            k1 = d2(1) + delta(10+6*m+6+length(d1)+1); k2 = d2(2) + delta(10+6*m+6+length(d1)+2); 
            p1 = 0; p2 = 0; k3 = 0;
            d2_lm = [k1; k2];
        elseif (length(d2) == 4) 
            k1 = d2(1) + delta(10+6*m+6+length(d1)+1); k2 = d2(2) + delta(10+6*m+6+length(d1)+2); 
            p1 = d2(3) + delta(10+6*m+6+length(d1)+3); p2 = d2(4) + delta(10+6*m+6+length(d1)+4); 
            k3 = 0;
            d2_lm = [k1; k2; p1; p2];
        elseif (length(d2) == 5) 
            k1 = d2(1) + delta(10+6*m+6+length(d1)+1); k2 = d2(2) + delta(10+6*m+6+length(d1)+2); 
            p1 = d2(3) + delta(10+6*m+6+length(d1)+3); p2 = d2(4) + delta(10+6*m+6+length(d1)+4); 
            k3 = d2(5) + delta(10+6*m+6+length(d1)+5);
            d2_lm = [k1; k2; p1; p2; k3];
        end


        xd2 = xu2.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1*xu2.*yu2 + p2*(r.^2 + 2*xu2.^2);
        yd2 = yu2.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*yu2.^2) + 2*p2*xu2.*yu2;
    else
        xd2 = xu2;
        yd2 = yu2;
    end

    %% Projected 2D Image Points of the 1st and 2nd Camera
    u1 = K1_lm(1,1)*xd1 + K1_lm(1,2)*yd1 + K1_lm(1,3).*ones(1,m*noPnts);
    v1 = K1_lm(2,2)*yd1 + K1_lm(2,3).*ones(1,m*noPnts);
    
    u2 = K2_lm(1,1)*xd2 + K2_lm(1,2)*yd2 + K2_lm(1,3).*ones(1,m*noPnts);
    v2 = K2_lm(2,2)*yd2 + K2_lm(2,3).*ones(1,m*noPnts);


    %% Distance between the re-projected points and the measured points
    dist_lm(1:4:4*m*noPnts,1) = u1 - x1(1,:);
    dist_lm(2:4:4*m*noPnts,1) = v1 - x1(2,:);
    dist_lm(3:4:4*m*noPnts,1) = u2 - x2(1,:);
    dist_lm(4:4:4*m*noPnts,1) = v2 - x2(2,:);
    
    
    % Re-projection Error
    rperr_lm = sqrt(dot(dist_lm,dist_lm)/(4*m*noPnts));

    
    if (rperr_lm <= rperr)
        param = [K1(1,1); K1(2,2); K1(1,3); K1(2,3); K1(1,2); 
                 K2(1,1); K2(2,2); K2(1,3); K2(2,3); K2(1,2); 
                 reshape(w1,3*m,1); reshape(t1,3*m,1);
                 w; t; d1; d2];
        
       if (((n > 1) && sqrt(dot(delta,delta)/dot(param,param)) < rerr))
            K1 = K1_lm; 
            K2 = K2_lm; 
            R1 = R1_lm; 
            t1 = t1_lm;
            R  = R_lm; 
            t  = t_lm;
            
            if (length(d1 > 0)) d1 = d1_lm; end
            if (length(d2 > 0)) d2 = d2_lm; end
            
            rperr = rperr_lm;
            break;
       end
        
        % Update
        K1 = K1_lm; 
        K2 = K2_lm;
        R1 = R1_lm; 
        w1 = w1_lm; 
        t1 = t1_lm; 
        R  = R_lm;
        w  = w_lm;
        t  = t_lm;
        
        if (length(d1) > 0) d1 = d1_lm; end
        if (length(d2) > 0) d2 = d2_lm; end
        
        dist  = dist_lm;
        rperr = rperr_lm;
       
        
        %% Compute the Jacobian
        for (i=1:m*noPnts)
             % Plane Count
            j = floor((i-1)/noPnts)+1;
            
            xxu1 = [xu1(i); yu1(i)];
            xxd1 = [xd1(i); yd1(i)];
            xxu2 = [xu2(i); yu2(i)];
            xxd2 = [xd2(i); yd2(i)];
            
            %% The derivative of a undistorted normalized point
            [dxxu1_dw1, dxxu1_dt1] = ...
                Compute_dxxu(XXw(:,i), R1(:,(j-1)*3+1:j*3), t1(:,j), w1(:,j));
            
            [dxxu2_dw1, dxxu2_dt1, dxxu2_dw, dxxu2_dt] = ...
                Compute_dxxu2(XXw(:,i), XXc1(:,i), R, t, w1(:,j), w);
            
            dxxu1_dw = zeros(2,3);  dxxu1_dt = zeros(2,3);
            
            
            %% The derivative of a distorted normalized point
            [dxxd1_dw1, dxxd1_dt1, dxxd1_dd1] = ...
                Compute_dxxd(xxu1, d1, dxxu1_dw1, dxxu1_dt1);
            
            [dxxd2_dw1, dxxd2_dt1, dxxd2_dw, dxxd2_dt, dxxd2_dd2] = ...
                Compute_dxxd2(xxu2, d2, dxxu2_dw1, dxxu2_dt1, dxxu2_dw, dxxu2_dt);
            
            dxxd1_dw  = zeros(2,3); 
            dxxd1_dt  = zeros(2,3); 
            dxxd1_dd2 = zeros(2,length(d2));
            dxxd2_dd1 = zeros(2,length(d1));
            
            
            %% The derivative of a distotred 2D pixel points
            [dxx1_dk1, dxx1_dw1, dxx1_dt1, dxx1_dd1] = ...
                Compute_dxx(xxd1, K1, d1, dxxd1_dw1, dxxd1_dt1, dxxd1_dd1);
                 
            [dxx2_dk2, dxx2_dw1, dxx2_dt1, dxx2_dw, dxx2_dt, dxx2_dd2] = ...
                Compute_dxx2(xxd2, K2, d2, dxxd2_dw1, dxxd2_dt1, dxxd2_dw, dxxd2_dt, dxxd2_dd2);
            
            dxx1_dk2 = zeros(2,5); 
            dxx1_dw  = zeros(2,3); 
            dxx1_dt  = zeros(2,3); 
            dxx1_dd2 = zeros(2,length(d2));
            
            dxx2_dk1 = zeros(2,5); 
            dxx2_dd1 = zeros(2,length(d1));
            
            dxx1_dw1_all = zeros(2,3*m);
            dxx1_dt1_all = zeros(2,3*m);
            dxx2_dw1_all = zeros(2,3*m);
            dxx2_dt1_all = zeros(2,3*m);
            
            dxx1_dw1_all(:,(j-1)*3+1:j*3) = dxx1_dw1;
            dxx1_dt1_all(:,(j-1)*3+1:j*3) = dxx1_dt1;
            dxx2_dw1_all(:,(j-1)*3+1:j*3) = dxx2_dw1;
            dxx2_dt1_all(:,(j-1)*3+1:j*3) = dxx2_dt1;
            
            
            %% Jacobian
            J(4*(i-1)+1:4*i, 1:noParam) = ...
                [dxx1_dk1 dxx1_dk2 dxx1_dw1_all dxx1_dt1_all ...
                 dxx1_dw dxx1_dt dxx1_dd1 dxx1_dd2;
                 dxx2_dk1 dxx2_dk2 dxx2_dw1_all dxx2_dt1_all ...
                 dxx2_dw dxx2_dt dxx2_dd1 dxx2_dd2];
        end
          
        % Compute the approximated Hessian matrix
        H = J'*J;
        
        if (n == 1)
            lambda = 0.001*trace(H)/noParam;
        else
             lambda = lambda/10;
        end
    else
        lambda = lambda*10;
    end
    
    
    % Apply the damping factor to the Hessian matrix
    H_lm = H + (lambda * eye(noParam, noParam));

    % Prevent the matrix from being singular
    if (rcond(H_lm) < eps)
        lambda = lambda*10;
        H_lm = H + (lambda * eye(noParam, noParam));
    end
    
    % Compute the updated parameters
    delta = -inv(H_lm)*(J'*dist(:));
end


%% Output
% Make the output same size as the input
if ((length(d1) > 0) && (d1Size(1) == 1))
    d1 = d1';
end

if ((length(d2) > 0) && (d2Size(1) == 1))
    d2 = d2';
end

if (nargout == 5)
    varargout = {K1, K2, R, t, rperr};
elseif (nargout == 7)
    varargout = {K1, d1, K2, d2, R, t, rperr};
end


%% Sub Functions
function [dxxu_dw, dxxu_dt] = Compute_dxxu(XXw, R, t, w)
    XXc = R*XXw + t;

    dxxu_dXXc = Compute_dxxu_dXXc(XXc);
    dXXc_dr   = Compute_dXXc_dr(XXw);
    dr_dw     = Compute_dr_dw(w);
    dXXc_dt   = eye(3,3);

    dxxu_dw = dxxu_dXXc*dXXc_dr*dr_dw;
    dxxu_dt = dxxu_dXXc*dXXc_dt;
    
    
function [dxxu2_dw1, dxxu2_dt1, dxxu2_dw, dxxu2_dt] = Compute_dxxu2(XXw, XXc1, R, t, w1, w)
    XXc2 = R*XXc1 + t;

    dxxu2_dXXc2 = Compute_dxxu_dXXc(XXc2);
    
    dXXc2_dr1  = Compute_dXXc2_dr1(XXw, R);
    dXXc2_dr   = Compute_dXXc_dr(XXc1);
    
    dr1_dw1    = Compute_dr_dw(w1);
    dr_dw      = Compute_dr_dw(w);
    
    dXXc2_dt1  = R;
    dXXc2_dt   = eye(3,3);

    dxxu2_dw1  = dxxu2_dXXc2*dXXc2_dr1*dr1_dw1;
    dxxu2_dt1  = dxxu2_dXXc2*dXXc2_dt1;
    dxxu2_dw   = dxxu2_dXXc2*dXXc2_dr*dr_dw;
    dxxu2_dt   = dxxu2_dXXc2*dXXc2_dt;
    
    
function [dxxu_dXXc] = Compute_dxxu_dXXc(XXc)
Xc = XXc(1);
Yc = XXc(2);
Zc = XXc(3);

dxxu_dXXc = zeros(2,3);
dxxu_dXXc(1,1) = 1/Zc;  
dxxu_dXXc(2,1) = 0; 
dxxu_dXXc(1,2) = 0;
dxxu_dXXc(2,2) = 1/Zc;
dxxu_dXXc(1,3) = -Xc/Zc^2;
dxxu_dXXc(2,3) = -Yc/Zc^2;
    
    
function [dXXc_dr] = Compute_dXXc_dr(XXw)
dXXc_dr = zeros(3,9);
dXXc_dr(1,1:3:end) = XXw;
dXXc_dr(2,2:3:end) = XXw;
dXXc_dr(3,3:3:end) = XXw;


function [dXXc2_dr1] = Compute_dXXc2_dr1(XXw, R)
dXXc2_dr1 = [XXw(1)*R XXw(2)*R XXw(3)*R];

    
function [dr_dw] = Compute_dr_dw(w)
wx = w(1); wy = w(2); wz = w(3);
theta = sqrt(wx^2 + wy^2 + wz^2);

if (theta < eps)
    dr_dw = [ 0  0  0;
              0  0  1;
              0 -1  0;
              0  0 -1;
              0  0  0;
              1  0  0;
              0  1  0;
             -1  0  0;
              0  0  0];
else
    wxh = wx/theta; wyh = wy/theta; wzh = wz/theta; 

    dsth_dw   = [wx*cos(theta)/theta wy*cos(theta)/theta wz*cos(theta)/theta];
    domcth_dw = [wx*sin(theta)/theta, wy*sin(theta)/theta, wz*sin(theta)/theta];

    dwxh_dw = [1/theta-wx^2/theta^3, -wx*wy/theta^3, -wx*wz/theta^3];
    dwyh_dw = [-wx*wy/theta^3, 1/theta-wy^2/theta^3, -wy*wz/theta^3];
    dwzh_dw = [-wx*wz/theta^3, -wy*wz/theta^3, 1/theta-wz^2/theta^3];

    dwxh2pwyh2_dw = [2*wx/theta^2-2*(wx^3+wy^2*wx)/theta^4, 2*wy/theta^2-2*(wy^3+wx^2*wy)/theta^4, -2*(wx^2*wz+wy^2*wz)/theta^4];
    dwxh2pwzh2_dw = [2*wx/theta^2-2*(wx^3+wz^2*wx)/theta^4, -2*(wx^2*wy+wz^2*wy)/theta^4, 2*wz/theta^2-2*(wz^3+wx^2*wz)/theta^4];
    dwyh2pwzh2_dw = [-2*(wy^2*wx+wz^2*wx)/theta^4, 2*wy/theta^2-2*(wy^3+wz^2*wy)/theta^4, 2*wz/theta^2-2*(wz^3+wy^2*wz)/theta^4];

    wxh2pwyh2 = wxh^2 + wyh^2;
    wxh2pwzh2 = wxh^2 + wzh^2;
    wyh2pwzh2 = wyh^2 + wzh^2;

    omcth = 1-cos(theta);
    sth   = sin(theta);

    dr_dw = [-omcth*dwyh2pwzh2_dw - wyh2pwzh2*domcth_dw;
              sth*dwzh_dw + wzh*dsth_dw + wyh*omcth*dwxh_dw + wxh*omcth*dwyh_dw + wxh*wyh*domcth_dw;
             -sth*dwyh_dw - wyh*dsth_dw + wzh*omcth*dwxh_dw + wxh*omcth*dwzh_dw + wxh*wzh*domcth_dw;
             -sth*dwzh_dw - wzh*dsth_dw + wyh*omcth*dwxh_dw + wxh*omcth*dwyh_dw + wxh*wyh*domcth_dw;
             -omcth*dwxh2pwzh2_dw - wxh2pwzh2*domcth_dw;
              sth*dwxh_dw + wxh*dsth_dw + wzh*omcth*dwyh_dw + wyh*omcth*dwzh_dw + wyh*wzh*domcth_dw;
              sth*dwyh_dw + wyh*dsth_dw + wzh*omcth*dwxh_dw + wxh*omcth*dwzh_dw + wxh*wzh*domcth_dw;
             -sth*dwxh_dw - wxh*dsth_dw + wzh*omcth*dwyh_dw + wyh*omcth*dwzh_dw + wyh*wzh*domcth_dw;
             -omcth*dwxh2pwyh2_dw - wxh2pwyh2*domcth_dw];
end

           
function [dxxd_dw, dxxd_dt, dxxd_dd] = Compute_dxxd(xxu, d, dxxu_dw, dxxu_dt)
xu = xxu(1);
yu = xxu(2);
r  = sqrt(xu^2 + yu^2);

dxu_dw = dxxu_dw(1,:);
dyu_dw = dxxu_dw(2,:);

dxu_dt = dxxu_dt(1,:);
dyu_dt = dxxu_dt(2,:);


% The derivative of a radial distortion function
%   c = 1 + k1*r^2 + k2*r^4 + k3*r^6;
if (length(d) >= 1)
    % First radial distortion coefficients, k1
    k1 = d(1);
    c = 1 + k1*r^2;

    dr2_dw = 2*xu*dxu_dw + 2*yu*dyu_dw;
    dr2_dt = 2*xu*dxu_dt + 2*yu*dyu_dt;
    dc_dk1 = r^2;

    dc_dw = k1*dr2_dw;
    dc_dt = k1*dr2_dt;
    dc_dk = dc_dk1;
end

 if (length(d) >= 2)
     % Second radial distortion coefficients, k2
     k2 = d(2);
     c = c + k2*r^4;

     dr4_dw = 2*r^2*dr2_dw;
     dr4_dt = 2*r^2*dr2_dt;
     dc_dk2 = r^4;

     dc_dw = dc_dw + k2*dr4_dw;
     dc_dt = dc_dt + k2*dr4_dt;
     dc_dk = [dc_dk dc_dk2];
 end

 if (length(d) == 5)
     % Third radial distortion coefficients, k3
     k3 = d(5);
     c = c + k3*r^6;

     dr6_dw = 3*r^4*dr2_dw;
     dr6_dt = 3*r^4*dr2_dt;
     dc_dk3 = r^6;

     dc_dw = dc_dw + k3*dr6_dw;
     dc_dt = dc_dt + k3*dr6_dt;
     dc_dk = [dc_dk dc_dk3];
 end

 if (length(d) > 0)
     % The derivative of a radially distorted normalized point
     dxr_dw = c*dxu_dw + xu*dc_dw;
     dyr_dw = c*dyu_dw + yu*dc_dw;

     dxr_dt = c*dxu_dt + xu*dc_dt;
     dyr_dt = c*dyu_dt + yu*dc_dt;

     dxr_dk = xu*dc_dk;
     dyr_dk = yu*dc_dk;
 end


 % The derivative of a tangentially distorted normalized point
 if (length(d) >= 4)
     p1 = d(3);
     p2 = d(4);
     
     dxt_dxxu = [2*p1*yu+6*p2*xu 2*p1*xu+2*p2*yu];
     dyt_dxxu = [2*p1*xu+2*p2*yu 6*p1*yu+2*p2*xu];

     dxxu_dw = [dxu_dw; dyu_dw];
     dxxu_dt = [dxu_dt; dyu_dt];

     dxt_dw = dxt_dxxu*dxxu_dw;
     dyt_dw = dyt_dxxu*dxxu_dw;

     dxt_dt = dxt_dxxu*dxxu_dt;
     dyt_dt = dyt_dxxu*dxxu_dt;

     dxt_dp = [    2*xu*yu r^2+2*xu^2 ];
     dyt_dp = [r^2+2*yu^2      2*xu*yu];
 end


% The derivative of a distorted normalized point
if (length(d) == 0)
    dxd_dw = dxu_dw;
    dyd_dw = dyu_dw;

    dxd_dt = dxu_dt;
    dyd_dt = dyu_dt;
    
    dxd_dd = [];
    dyd_dd = [];

elseif (length(d) <= 2)
    dxd_dw = dxr_dw;
    dyd_dw = dyr_dw;

    dxd_dt = dxr_dt;
    dyd_dt = dyr_dt;

    dxd_dd = dxr_dk;
    dyd_dd = dyr_dk;

elseif (length(d) <= 5)
    dxd_dw = dxr_dw + dxt_dw;
    dyd_dw = dyr_dw + dyt_dw;

    dxd_dt = dxr_dt + dxt_dt;
    dyd_dt = dyr_dt + dyt_dt;

    if (length(d) == 4)
        dxd_dd = [dxr_dk dxt_dp];
        dyd_dd = [dyr_dk dyt_dp];
    elseif (length(d) == 5)
        dxd_dd = [dxr_dk(1,1:2) dxt_dp dxr_dk(1,3)];
        dyd_dd = [dyr_dk(1,1:2) dyt_dp dyr_dk(1,3)];
    end
end

dxxd_dw = [dxd_dw; dyd_dw];
dxxd_dt = [dxd_dt; dyd_dt];
dxxd_dd = [dxd_dd; dyd_dd];


function [dxxd2_dw1, dxxd2_dt1, dxxd2_dw, dxxd2_dt, dxxd2_dd2] = Compute_dxxd2(xxu2, d2, dxxu2_dw1, dxxu2_dt1, dxxu2_dw, dxxu2_dt)
xu2 = xxu2(1);
yu2 = xxu2(2);
r  = sqrt(xu2^2 + yu2^2);

dxu2_dw1 = dxxu2_dw1(1,:);
dyu2_dw1 = dxxu2_dw1(2,:);

dxu2_dt1 = dxxu2_dt1(1,:);
dyu2_dt1 = dxxu2_dt1(2,:);

dxu2_dw = dxxu2_dw(1,:);
dyu2_dw = dxxu2_dw(2,:);

dxu2_dt = dxxu2_dt(1,:);
dyu2_dt = dxxu2_dt(2,:);


% The derivative of a radial distortion function
%   c = 1 + k1*r^2 + k2*r^4 + k3*r^6;
if (length(d2) >= 1)
    % First radial distortion coefficients, k1
    k1 = d2(1);
    c = 1 + k1*r^2;

    dr2_dw1 = 2*xu2*dxu2_dw1 + 2*yu2*dyu2_dw1;
    dr2_dt1 = 2*xu2*dxu2_dt1 + 2*yu2*dyu2_dt1;
    dr2_dw  = 2*xu2*dxu2_dw  + 2*yu2*dyu2_dw;
    dr2_dt  = 2*xu2*dxu2_dt  + 2*yu2*dyu2_dt;
    dc_dk1  = r^2;

    dc_dw1 = k1*dr2_dw1;
    dc_dt1 = k1*dr2_dt1;
    dc_dw  = k1*dr2_dw;
    dc_dt  = k1*dr2_dt;
    dc_dk  = dc_dk1;
end

 if (length(d2) >= 2)
     % Second radial distortion coefficients, k2
     k2 = d2(2);
     c = c + k2*r^4;

     dr4_dw1 = 2*r^2*dr2_dw1;
     dr4_dt1 = 2*r^2*dr2_dt1;
     dr4_dw  = 2*r^2*dr2_dw;
     dr4_dt  = 2*r^2*dr2_dt;
     dc_dk2  = r^4;

     dc_dw1 = dc_dw1 + k2*dr4_dw1;
     dc_dt1 = dc_dt1 + k2*dr4_dt1;
     dc_dw  = dc_dw  + k2*dr4_dw;
     dc_dt  = dc_dt  + k2*dr4_dt;
     dc_dk = [dc_dk dc_dk2];
 end

 if (length(d2) == 5)
     % Third radial distortion coefficients, k3
     k3 = d2(5);
     c = c + k3*r^6;

     dr6_dw1 = 3*r^4*dr2_dw1;
     dr6_dt1 = 3*r^4*dr2_dt1;
     dr6_dw  = 3*r^4*dr2_dw;
     dr6_dt  = 3*r^4*dr2_dt;
     dc_dk3 = r^6;

     dc_dw1 = dc_dw1 + k3*dr6_dw1;
     dc_dt1 = dc_dt1 + k3*dr6_dt1;
     dc_dw  = dc_dw  + k3*dr6_dw;
     dc_dt  = dc_dt  + k3*dr6_dt;
     dc_dk = [dc_dk dc_dk3];
 end

 if (length(d2) > 0)
     % The derivative of a radially distorted normalized point
     dxr_dw1 = c*dxu2_dw1 + xu2*dc_dw1;
     dyr_dw1 = c*dyu2_dw1 + yu2*dc_dw1;
     dxr_dt1 = c*dxu2_dt1 + xu2*dc_dt1;
     dyr_dt1 = c*dyu2_dt1 + yu2*dc_dt1;
     
     dxr_dw  = c*dxu2_dw + xu2*dc_dw;
     dyr_dw  = c*dyu2_dw + yu2*dc_dw;
     dxr_dt  = c*dxu2_dt + xu2*dc_dt;
     dyr_dt  = c*dyu2_dt + yu2*dc_dt;

     dxr_dk = xu2*dc_dk;
     dyr_dk = yu2*dc_dk;
 end


 % The derivative of a tangentially distorted normalized point
 if (length(d2) >= 4)
     p1 = d2(3);
     p2 = d2(4);
     
     dxt_dxxu2 = [2*p1*yu2+6*p2*xu2 2*p1*xu2+2*p2*yu2];
     dyt_dxxu2 = [2*p1*xu2+2*p2*yu2 6*p1*yu2+2*p2*xu2];

     dxxu2_dw1 = [dxu2_dw1; dyu2_dw1];
     dxxu2_dt1 = [dxu2_dt1; dyu2_dt1];
     dxxu2_dw  = [dxu2_dw;  dyu2_dw ];
     dxxu2_dt  = [dxu2_dt;  dyu2_dt ];

     dxt_dw1 = dxt_dxxu2*dxxu2_dw1;
     dyt_dw1 = dyt_dxxu2*dxxu2_dw1;
     dxt_dt1 = dxt_dxxu2*dxxu2_dt1;
     dyt_dt1 = dyt_dxxu2*dxxu2_dt1;
     
     dxt_dw = dxt_dxxu2*dxxu2_dw;
     dyt_dw = dyt_dxxu2*dxxu2_dw;
     dxt_dt = dxt_dxxu2*dxxu2_dt;
     dyt_dt = dyt_dxxu2*dxxu2_dt;

     dxt_dp = [    2*xu2*yu2 r^2+2*xu2^2 ];
     dyt_dp = [r^2+2*yu2^2      2*xu2*yu2];
 end


% The derivative of a distorted normalized point
if (length(d2) == 0)
    dxd2_dw1 = dxu2_dw1;
    dyd2_dw1 = dyu2_dw1;
    dxd2_dt1 = dxu2_dt1;
    dyd2_dt1 = dyu2_dt1;
    
    dxd2_dw  = dxu2_dw;
    dyd2_dw  = dyu2_dw;
    dxd2_dt  = dxu2_dt;
    dyd2_dt  = dyu2_dt;
    
    dxd2_dd2 = [];
    dyd2_dd2 = [];

elseif (length(d2) <= 2)
    dxd2_dw1 = dxr_dw1;
    dyd2_dw1 = dyr_dw1;
    dxd2_dt1 = dxr_dt1;
    dyd2_dt1 = dyr_dt1;
    
    dxd2_dw  = dxr_dw;
    dyd2_dw  = dyr_dw;
    dxd2_dt  = dxr_dt;
    dyd2_dt  = dyr_dt;

    dxd2_dd2 = dxr_dk;
    dyd2_dd2 = dyr_dk;

elseif (length(d2) <= 5)
    dxd2_dw1 = dxr_dw1 + dxt_dw1;
    dyd2_dw1 = dyr_dw1 + dyt_dw1;
    dxd2_dt1 = dxr_dt1 + dxt_dt1;
    dyd2_dt1 = dyr_dt1 + dyt_dt1;
    
    dxd2_dw  = dxr_dw + dxt_dw;
    dyd2_dw  = dyr_dw + dyt_dw;
    dxd2_dt  = dxr_dt + dxt_dt;
    dyd2_dt  = dyr_dt + dyt_dt;

    if (length(d2) == 4)
        dxd2_dd2 = [dxr_dk dxt_dp];
        dyd2_dd2 = [dyr_dk dyt_dp];
    elseif (length(d2) == 5)
        dxd2_dd2 = [dxr_dk(1,1:2) dxt_dp dxr_dk(1,3)];
        dyd2_dd2 = [dyr_dk(1,1:2) dyt_dp dyr_dk(1,3)];
    end
end

dxxd2_dw1 = [dxd2_dw1; dyd2_dw1];
dxxd2_dt1 = [dxd2_dt1; dyd2_dt1];
dxxd2_dw  = [dxd2_dw ; dyd2_dw ];
dxxd2_dt  = [dxd2_dt ; dyd2_dt ];
dxxd2_dd2 = [dxd2_dd2; dyd2_dd2];

    
function [dxx_dk, dxx_dw, dxx_dt, dxx_dd] = Compute_dxx(xxd, K, d, dxxd_dw, dxxd_dt, dxxd_dd)
xd = xxd(1);
yd = xxd(2);

fx = K(1,1);
fy = K(2,2);
u0 = K(1,3);
v0 = K(2,3);
s  = K(1,2);

dxd_dw = dxxd_dw(1,:);
dyd_dw = dxxd_dw(2,:);

dxd_dt = dxxd_dt(1,:);
dyd_dt = dxxd_dt(2,:);

du_dfx = xd; du_dfy =  0;
du_du0 =  1; du_dv0 =  0;
dv_dfx =  0; dv_dfy = yd;
dv_du0 =  0; dv_dv0 =  1; 

if (s == 0)
    du_ds = 0;
    dv_ds = 0;
else
    du_ds = yd;
    dv_ds =  0;
end

if (length(d) ~= 0)
    dxd_dd = dxxd_dd(1,:);
    dyd_dd = dxxd_dd(2,:);

    du_dd = fx*dxd_dd + s*dyd_dd;
    dv_dd = fy*dyd_dd;
end


du_dk = [du_dfx du_dfy du_du0 du_dv0 du_ds];
dv_dk = [dv_dfx dv_dfy dv_du0 dv_dv0 dv_ds];

du_dw = fx*dxd_dw + s*dyd_dw;
dv_dw = fy*dyd_dw;

du_dt = fx*dxd_dt + s*dyd_dt;
dv_dt = fy*dyd_dt;

dxx_dk = [du_dk; dv_dk];
dxx_dw = [du_dw; dv_dw];
dxx_dt = [du_dt; dv_dt];

if (length(d) == 0)
    dxx_dd = [];
else
    dxx_dd = [du_dd; dv_dd];
end


function [dxx2_dk2, dxx2_dw1, dxx2_dt1, dxx2_dw, dxx2_dt, dxx2_dd2] = ...
    Compute_dxx2(xxd2, K2, d2, dxxd2_dw1, dxxd2_dt1, dxxd2_dw, dxxd2_dt, dxxd2_dd2)
xd2 = xxd2(1);
yd2 = xxd2(2);

fx = K2(1,1);
fy = K2(2,2);
u0 = K2(1,3);
v0 = K2(2,3);
s  = K2(1,2);

dxd2_dw1 = dxxd2_dw1(1,:);
dyd2_dw1 = dxxd2_dw1(2,:);

dxd2_dt1 = dxxd2_dt1(1,:);
dyd2_dt1 = dxxd2_dt1(2,:);

dxd2_dw  = dxxd2_dw(1,:);
dyd2_dw  = dxxd2_dw(2,:);

dxd2_dt  = dxxd2_dt(1,:);
dyd2_dt  = dxxd2_dt(2,:);

du2_dfx = xd2; du2_dfy =   0;
du2_du0 =   1; du2_dv0 =   0;
dv2_dfx =   0; dv2_dfy = yd2;
dv2_du0 =   0; dv2_dv0 =   1; 

if (s == 0)
    du2_ds = 0;
    dv2_ds = 0;
else
    du2_ds = yd2;
    dv2_ds =   0;
end

if (length(d2) ~= 0)
    dxd2_dd2 = dxxd2_dd2(1,:);
    dyd2_dd2 = dxxd2_dd2(2,:);

    du2_dd2 = fx*dxd2_dd2 + s*dyd2_dd2;
    dv2_dd2 = fy*dyd2_dd2;
end


du2_dk2 = [du2_dfx du2_dfy du2_du0 du2_dv0 du2_ds];
dv2_dk2 = [dv2_dfx dv2_dfy dv2_du0 dv2_dv0 dv2_ds];

du2_dw1 = fx*dxd2_dw1 + s*dyd2_dw1;
dv2_dw1 = fy*dyd2_dw1;

du2_dt1 = fx*dxd2_dt1 + s*dyd2_dt1;
dv2_dt1 = fy*dyd2_dt1;

du2_dw  = fx*dxd2_dw + s*dyd2_dw;
dv2_dw  = fy*dyd2_dw;

du2_dt  = fx*dxd2_dt + s*dyd2_dt;
dv2_dt  = fy*dyd2_dt;

dxx2_dk2 = [du2_dk2; dv2_dk2];
dxx2_dw1 = [du2_dw1; dv2_dw1];
dxx2_dt1 = [du2_dt1; dv2_dt1];
dxx2_dw  = [du2_dw ; dv2_dw ];
dxx2_dt  = [du2_dt ; dv2_dt ];

if (length(d2) == 0)
    dxx2_dd2 = [];
else
    dxx2_dd2 = [du2_dd2; dv2_dd2];
end