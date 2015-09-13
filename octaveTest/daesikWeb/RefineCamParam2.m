% RefineCamParam2 - refines the camera parameters for planar object based 
%                   calibration with LM (Levenberg-Marquardt) nonlinear 
%                   least squares algorithm.
%
% Usage:
%        [K, R, t, rperr] = RefineCamParam2(x, X, K, R, t)
%        [K, R, t, rperr] = RefineCamParam2(x, X, K, R, t, rerr, iter)
%        [K, d, R, t, rperr] = RefineCamParam2(x, X, K, d, R, t)
%        [K, d, R, t, rperr] = RefineCamParam2(x, X, K, d, R, t, rerr, iter)
%           
%
% Input:
%        x: 2xn, nx2, 3xn or nx3 image points
%           (n is the total number of points)
%        X: 2xn, nx2, 3xn or nx3 planar object points
%           (n is the total number of points)
%        K: initial 3x3 camera intrinsic parameters matrix
%           (If the skew is 0, the function does not refine the skew.)
%        d: initial 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%           lens distortion vector
%           (1x1        : k1, 
%            2x1(or 1x2): k1, k2,
%            4x1(or 4x1): k1, k2, p1, p2, 
%            5x1(or 5x1): k1, k2, p1, p2, k3,
%            where k1, k2, k3 are the radial distortion parameters and
%            p1, p2 are the tangential distortion parameters.) 
%        R: initial 3x3m or 3mx3 rotation matrices 
%           (m is the number of planes)
%        t: initial 3xm or mx3 translation vectors 
%           (m is the number of planes. If m is 3, it is assumed that 
%            the translation vectors are stacked in the same way as R. For 
%            example, the R is 3x3m, t is regarded as 3xm and vice versa.)
%        rerr: relative error between the last and preceding iteration.
%              (default: 2^(-52) (It is close to 2.2204e-016.))
%        iter: the number of maximum iteration (default : 30)
%
% Output:
%         K: refined 3x3 camera intrinsic parameters matrix
%         d: refined 1x1, 2x1(or 1x2), 4x1(or 1x4), or 5x1(or 1x5) 
%            lens distortion vector
%         R: refined 3x3m or 3mx3 rotation matrix
%         t: refined 3xm or mx3 translation vectors
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
% Jun. 2011 - Lens distortion parameters are added.
%             Input and output arguments are modified.
% Feb. 2012 - Bug (jacobian computation) is fixed.
%             Input and ouput arguments are modified.


function varargout = RefineCamParam2(varargin)

if (nargin == 5)
    [x, X, K, R, t] = varargin{:};
    d = [];
    rerr = 2^(-52);
    iter = 30;
elseif (nargin == 6)
    [x, X, K, d, R, t] = varargin{:};
    rerr = 2^(-52);
    iter = 30;
elseif (nargin == 7)
    [x, X, K, R, t, rerr, iter] = varargin{:};
    d = [];
elseif (nargin == 8)
    [x, X, K, d, R, t, rerr, iter] = varargin{:};
end


%% The Number of Planes and Points
m = length(R)/3; % number of plane

xSize  = size(x);
XSize  = size(X);
RSize  = size(R);
tSize  = size(t);
dSize  = size(d);
noPnts = length(x)/m;


%% Image Points
if (xSize(2) == 2)
    x = x';
elseif (xSize(1) == 3)
    x = x(1:2,:);
elseif (xSize(2) == 3)
    x = x(:,1:2)';
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
if (RSize(2) == 3)
    R = R';

    for (j=1:m)
        R(:,(j-1)*3+1:j*3) = R(:,(j-1)*3+1:j*3)';
    end
end


%% Translation Vector
if ((tSize(1) ~= 3) && (tSize(2) == 3))
    t = t';
elseif ((tSize(1) == 3) && (tSize(2) == 3) && (RSize(2) == 3))
    t = t';
end


%% Distortion Vector
if (dSize(1) == 1)
    d = d';
end


%% Convert the 3x3 rotation matrix into the 3x1 vector of Rodigrues representation
for (j=1:m)
    theta  = acos((trace(R(:,(j-1)*3+1:j*3))-1)/2);
    if (theta < eps)
        w(:,j) = [0 0 0]';
    else
        w(:,j) = (theta/(2*sin(theta)))*[R(3,(j-1)*3+2)-R(2,(j-1)*3+3); R(1,(j-1)*3+3)-R(3,(j-1)*3+1); R(2,(j-1)*3+1)-R(1,(j-1)*3+2)];
    end
end


%% LM(Levenberg-Marquardt) nonlinear least squares algorithm
% The number of parameters
%  - fx, fy, u0, v0, s, m*(wx, wy, wz), m*(tx, ty, tz), (k1, k2, p1, p2, k3)
noParam = 5 + 6*m + length(d);
param   = zeros(noParam, 1);
K_lm    = zeros(3,3);
w_lm    = zeros(3,m);
R_lm    = zeros(3,3*m);
J       = zeros(2*m*noPnts, noParam);  
dist_lm = zeros(2*m*noPnts,       1);
delta   = zeros(noParam, 1);
rperr  = inf;    % initial error


for n = 1:iter
    %% Camera Intrinsic Parameters Matrix
    K_lm(1,1) = K(1,1) + delta(1);  % fx
    K_lm(2,2) = K(2,2) + delta(2);  % fy  
    K_lm(1,3) = K(1,3) + delta(3);  % u0
    K_lm(2,3) = K(2,3) + delta(4);  % v0
    K_lm(1,2) = K(1,2) + delta(5);  % s
	K_lm(3,3) = 1;
    
    
    for (j=1:m)
        %% Convert the 3x1 vector of the Rodigrues representation 
        %% into the 3x3 rotation matrix
        w_lm(1, j) = w(1, j) + delta(5 + 3*(j-1)+1);  
        w_lm(2, j) = w(2, j) + delta(5 + 3*(j-1)+2);  
        w_lm(3, j) = w(3, j) + delta(5 + 3*(j-1)+3);
        
        wx = w_lm(1, j); wy = w_lm(2, j); wz = w_lm(3, j);
        theta = sqrt(wx^2 + wy^2 + wz^2);
        
        if (theta < eps)
            R_lm(:,(j-1)*3+1:j*3) = eye(3);
        else
            wh_sk = [  0 -wz  wy;
                      wz   0 -wx;
                     -wy  wx   0]/theta;

            %% 3x3 Rotation Matrix
            R_lm(:,(j-1)*3+1:j*3) = eye(3) + sin(theta)*wh_sk + (1-cos(theta))*(wh_sk*wh_sk);
        end
        
        
        %% Translation Vector
        t_lm(1, j) = t(1, j) + delta(5 + 3*m + 3*(j-1)+1);
        t_lm(2, j) = t(2, j) + delta(5 + 3*m + 3*(j-1)+2);
        t_lm(3, j) = t(3, j) + delta(5 + 3*m + 3*(j-1)+3);
    
        
        %% 3D points represented with respect to the camera coordinate frame
        XXc(:,(j-1)*noPnts+1:j*noPnts) = [R_lm(:,(j-1)*3+1:j*3), t_lm(:, j)]*[XXw(:,(j-1)*noPnts+1:j*noPnts); ones(1,noPnts)];
    end
        
    
    %% undistorted normalized points
    xu = XXc(1,:)./XXc(3,:);
    yu = XXc(2,:)./XXc(3,:);


    if (length(d) > 0)
        r = sqrt(xu.^2 + yu.^2);


        %% Lens Distortion Parameters
        if     (length(d) == 1) 
            k1 = d(1) + delta(5+6*m+1); k2 = 0 ; 
            p1 = 0; p2 = 0; k3 = 0 ;
            d_lm = [k1];
        elseif (length(d) == 2) 
            k1 = d(1) + delta(5+6*m+1); k2 = d(2) + delta(5+6*m+2); 
            p1 = 0; p2 = 0; k3 = 0;
            d_lm = [k1; k2];
        elseif (length(d) == 4) 
            k1 = d(1) + delta(5+6*m+1); k2 = d(2) + delta(5+6*m+2); 
            p1 = d(3) + delta(5+6*m+3); p2 = d(4) + delta(5+6*m+4); 
            k3 = 0;
            d_lm = [k1; k2; p1; p2];
        elseif (length(d) == 5) 
            k1 = d(1) + delta(5+6*m+1); k2 = d(2) + delta(5+6*m+2); 
            p1 = d(3) + delta(5+6*m+3); p2 = d(4) + delta(5+6*m+4); 
            k3 = d(5) + delta(5+6*m+5);
            d_lm = [k1; k2; p1; p2; k3];
        end


        xd = xu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1*xu.*yu + p2*(r.^2 + 2*xu.^2);
        yd = yu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*yu.^2) + 2*p2*xu.*yu;
    else
        xd = xu;
        yd = yu;
    end

    u = K_lm(1,1)*xd + K_lm(1,2)*yd + K_lm(1,3).*ones(1,m*noPnts);
    v = K_lm(2,2)*yd + K_lm(2,3).*ones(1,m*noPnts);


    % Distance between the re-projected points and the measured points
    dist_lm(1:2:2*m*noPnts,1) = u - x(1,:);
    dist_lm(2:2:2*m*noPnts,1) = v - x(2,:);
    
    
    % Re-projection Error
    rperr_lm = sqrt(dot(dist_lm,dist_lm)/(2*m*noPnts));

    
    if (rperr_lm <= rperr)
        param = [K(1,1); K(2,2); K(1,3); K(2,3); K(1,2); 
                 reshape(w,3*m,1); reshape(t,3*m,1); d];
        
       if (((n > 1) && sqrt(dot(delta,delta)/dot(param,param)) < rerr))
            K = K_lm;
            R = R_lm;
            t = t_lm;
            
            if (length(d) > 0)
                d = d_lm;
            end
            
            rperr = rperr_lm;
            break;
       end
        
        % Update
        K = K_lm;
        w = w_lm;
        R = R_lm;
        t = t_lm;
        
        
        if (length(d) > 0)
            d = d_lm;
        end
        
        dist  = dist_lm;
        rperr = rperr_lm;
       
        
        % Compute the Jacobian
        for (i=1:m*noPnts)
             % Plane Count
            j = floor((i-1)/noPnts)+1;
            
            xxu = [xu(i); yu(i)];
            xxd = [xd(i); yd(i)];
            
            % The derivative of a undistorted normalized point
            [dxxu_dw, dxxu_dt] = Compute_dxxu(XXw(:,i), R(:,(j-1)*3+1:j*3), t(:,j), w(:,j));
            
            % The derivative of a distorted normalized point
            [dxxd_dw, dxxd_dt, dxxd_dd] = Compute_dxxd(xxu, d, dxxu_dw, dxxu_dt);
            
            % The derivative of a distotred 2D pixel points
            [dxx_dk, dxx_dw, dxx_dt, dxx_dd] = Compute_dxx(xxd, K, d, dxxd_dw, dxxd_dt, dxxd_dd);
            
            dxx_dw_all = zeros(2,3*m);
            dxx_dt_all = zeros(2,3*m);
            dxx_dw_all(:,(j-1)*3+1:j*3) = dxx_dw;
            dxx_dt_all(:,(j-1)*3+1:j*3) = dxx_dt;
            
            % Jacobian
            J(2*i-1:2*i, 1:noParam) = [dxx_dk dxx_dw_all dxx_dt_all dxx_dd];
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
if (RSize(2) == 3)
    R = R';
    
    for (j=1:m)
        R((j-1)*3+1:j*3,:) = R((j-1)*3+1:j*3,:)';
    end
end

% Translation Vector
if ((tSize(1) ~= 3) && (tSize(2) == 3))
    t = t';
elseif ((tSize(1) == 3) && (tSize(2) == 3) && (RSize(2) == 3))
    t = t';
end

if ((length(d) > 0) && (dSize(1) == 1))
    d = d';
end

if (nargout == 4)
    varargout = {K, R, t, rperr};
elseif (nargout == 5)
    varargout = {K, d, R, t, rperr};
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