% RefineCamParamPlane - refines the camera parameters for planar object
%   based calibration with LM (Levenberg-Marquardt) nonlinear least squares
%   algorithm.
%
% Usage
%   [K, d, R, t, rperr] = RefineCamParamPlane(x, X, K, d, R, t, K_flag)
%   [K, d, R, t, rperr] = RefineCamParamPlane(x, X, K, d, R, t, K_flag, rerr, maxiter)
%
% Input
%   x: 2xn, nx2, 3xn or nx3 image points
%      (n is the total number of points)
%      If some of the point is not detected, it must be set to [-1, -1].
%   X: 2xn, nx2, 3xn or nx3 planar object points
%      (n is the total number of points)
%   K: initial 3x3 camera intrinsic parameters matrix
%   d: initial 1x1, 2x1(or 1x2), 3x1(or 1x3), 4x1(or 1x4), or 5x1(or 1x5)
%      lens distortion vector. If it is empty(i.e., []), it is ignored.
%      (1x1        : k1, 
%       2x1(or 1x2): k1, k2,
%       3x1(or 1x3): k1, k2, k3
%       4x1(or 4x1): k1, k2, p1, p2, 
%       5x1(or 5x1): k1, k2, k3, p1, p2,
%       where k1, k2, k3 are the radial distortion parameters and
%             p1, p2 are the tangential distortion parameters.) 
%   R: initial 3x3m or 3mx3 rotation matrices (m is the number of planes)
%   t: initial 3xm or mx3 translation vectors (m is the number of planes.
%      If m is 3, it is assumed that the translation vectors are stacked
%      in the same way as R. For example, the R is 3x3m, t is regarded as
%      3xm and vice versa.)
%   K_flag: flags to refine the intrinsic parameters
%       [0 0 0 0 0]: no refinement
%       [1 0 0 0 0]: fx = fy
%       [1 1 0 0 0]: fx, fy
%       [1 0 1 1 0]: fx = fy, x0, y0
%       [1 1 1 1 0]: fx, fy, x0, y0
%       [1 1 1 1 1]: fx, fy, x0, y0, skew
%   d_flag: flags to refine the distortion parameters
%       0: no refinement (The initial values are remained.)
%       1: the distortion parameters are refined.
%   rerr: relative error between the last and preceding iteration.
%         (default: 10^(-6))
%   iter: the number of maximum iteration (default : 200)
%
% Output
%   K: refined 3x3 camera intrinsic parameters matrix
%   d: refined 1x1, 2x1, 3x1, 4x1, or 5x1 lens distortion vector
%   R: refined 3x3m or 3mx3 rotation matrix
%   t: refined 3xm or mx3 translation vectors
%   rperr: re-projection error
%   noiter: number of iterations
%
%
% Kim, Daesik, Ph. D.
% E-mail  : daesik80[at]gmail.com
%
% Jun. 2010 - Original version.
% Jun. 2011 - Lens distortion parameters are added.
%             Input and output arguments are modified.
% Feb. 2012 - Bug (jacobian computation) is fixed.
%             Input and ouput arguments are modified.
% Jul. 2019 - Relative error computation method is changed.
% Nov. 2019 - Default values are changed.
%             Flag option for intrinsic parameters refinement is added.
%             Radial distortion only (1x3 or 3x1) option is supported.
%             'lambda' computation method is changed.
% Apr. 2023 - Flag option for distortion parameters refinement is added.


function varargout = RefineCamParamPlane(varargin)

if (nargin == 8)
    [x, X, K, d, R, t, K_flag, d_flag] = varargin{:};
    rerr = 10^(-6);
    maxiter = 200;
elseif (nargin == 10)
    [x, X, K, d, R, t, K_flag, d_flag, rerr, maxiter] = varargin{:};
end


%% The number of planes (views) and points
m = length(R)/3; % number of planes (views)

x_size  = size(x);
X_size  = size(X);
R_size  = size(R);
t_size  = size(t);
d_size  = size(d);
num_p_all = length(x); % number of all points
num_p = num_p_all/m; % number of points on each plane


%% Image points
if (x_size(2) == 2)
    x = x';
elseif (x_size(1) == 3)
    x = x(1:2,:);
elseif (x_size(2) == 3)
    x = x(:,1:2)';
end


%% Object Points
if (X_size(1) == 2)
    XYZw = [X; zeros(1,m*num_p)];
elseif (X_size(2) == 2)
    XYZw = [X'; zeros(1,m*num_p)];
elseif (X_size(1) == 3)
    XYZw = [X(1:2,:); zeros(1,m*num_p)];
elseif (X_size(2) == 3)
    XYZw = [X(:,1:2)'; zeros(1,m*num_p)];
end


%% Rotation Matrix
if (R_size(2) == 3)
    R = R';

    for j=1:m
        R(:,(j-1)*3+1:j*3) = R(:,(j-1)*3+1:j*3)';
    end
end


%% Translation Vector
if ((t_size(1) ~= 3) && (t_size(2) == 3))
    t = t';
elseif ((t_size(1) == 3) && (t_size(2) == 3) && (R_size(2) == 3))
    t = t';
end


%% Distortion Vector
if (d_size(1) == 1)
    d = d';
end


%% Convert the 3x3 rotation matrix into the 3x1 vector of Rodigrues representation
w = zeros(3, m);
for j=1:m
    w(:,j) = Rot2Rod(R(:,(j-1)*3+1:j*3));
end


%% Check the 2D points if it is detected (valid)
is_valid = zeros(1,num_p_all);
valid_idx_temp = zeros(1,num_p_all);
num_valid_p = 0;
for i=1:num_p_all
    if ((x(1,i) ~= -1) && (x(2,i) ~= -1))
        num_valid_p = num_valid_p + 1;
        valid_idx_temp(1,num_valid_p) = i;
        is_valid(i) = 1;
    else
        is_valid(i) = 0;
    end
end
valid_idx = valid_idx_temp(1,1:num_valid_p);


%% LM(Levenberg-Marquardt) nonlinear least squares algorithm
% The number of parameters
%  (fx, fy, u0, v0, s), m*(wx, wy, wz), m*(tx, ty, tz), (k1, k2, (k3), p1, p2)
num_Kparam = sum(K_flag);
if (d_flag == 1)
    num_Param = num_Kparam + 6*m + length(d);
else
    num_Param = num_Kparam + 6*m;
end
    
w_lm    = zeros(3, m);
R_lm    = zeros(3, 3*m);
t_lm    = zeros(3, m);
J       = zeros(2*m*num_p, num_Param);  
dist_lm = zeros(2*m*num_p,       1);
delta   = zeros(num_Param, 1);
rperr  = inf;    % initial error


for n = 1:maxiter
    noIter = n;
    
    %% 
    % Camera Intrinsic Parameters Matrix
    K_lm = Update_K(K, K_flag, delta(1:num_Kparam));
   
    
    for j=1:m
        %% Convert the 3x1 vector of the Rodigrues representation 
        %%  into the 3x3 rotation matrix
        w_lm(1, j) = w(1, j) + delta(num_Kparam + 3*(j-1)+1);  
        w_lm(2, j) = w(2, j) + delta(num_Kparam + 3*(j-1)+2);  
        w_lm(3, j) = w(3, j) + delta(num_Kparam + 3*(j-1)+3);
        
        R_lm(:,(j-1)*3+1:j*3) = Rod2Rot(w_lm(:, j)); 
       
        
        %% Translation Vector
        t_lm(1, j) = t(1, j) + delta(num_Kparam + 3*m + 3*(j-1)+1);
        t_lm(2, j) = t(2, j) + delta(num_Kparam + 3*m + 3*(j-1)+2);
        t_lm(3, j) = t(3, j) + delta(num_Kparam + 3*m + 3*(j-1)+3);
    
        
        %% 3D points represented with respect to the camera coordinate frame
        XYZc(:,(j-1)*num_p+1:j*num_p) = [R_lm(:,(j-1)*3+1:j*3), t_lm(:, j)]*[XYZw(:,(j-1)*num_p+1:j*num_p); ones(1,num_p)];
    end
        
    
    %% undistorted normalized points
    xu = XYZc(1,:)./XYZc(3,:);
    yu = XYZc(2,:)./XYZc(3,:);

    if (~isempty(d))
        r = sqrt(xu.^2 + yu.^2);
        
        d_delta = delta(num_Param - length(d) + 1:num_Param);
        [k1, k2, k3, p1, p2, d_lm] = Update_d(d, d_flag, d_delta);

        xd = xu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + 2*p1*xu.*yu + p2*(r.^2 + 2*xu.^2);
        yd = yu.*(1 + k1*r.^2 + k2*r.^4 + k3*r.^6) + p1*(r.^2 + 2*yu.^2) + 2*p2*xu.*yu;
    else
        xd = xu;
        yd = yu;
    end

    u = K_lm(1,1)*xd + K_lm(1,2)*yd + K_lm(1,3).*ones(1,m*num_p);
    v = K_lm(2,2)*yd + K_lm(2,3).*ones(1,m*num_p);


    %% Distance between the re-projected points and the measured points
    dist_lm(1:2:2*num_valid_p,1) = u(1,valid_idx) - x(1,valid_idx);
    dist_lm(2:2:2*num_valid_p,1) = v(1,valid_idx) - x(2,valid_idx);
    
    
    %% Re-projection Error
    rperr_lm = sqrt(dot(dist_lm,dist_lm)/(2*num_valid_p));

    
    if (rperr_lm <= rperr)
        delta_rperr = rperr - rperr_lm;
       
       
        % Update
        K = K_lm;
        w = w_lm;
        R = R_lm;
        t = t_lm;
        
        
        if (~isempty(d))
            d = d_lm;
        end
        
        dist  = dist_lm;
        rperr = rperr_lm;
       
        % Display the Re-projection Error
        disp(['Iteration: ', num2str(n), ', Re-projection Error: ', num2str(rperr)]);
        
        % If the relative error is small, break.
        if (delta_rperr < rerr)
            break;
        end
        
        
        %% Compute the Jacobian
        k = 0;
        for i=1:num_p_all
            if (is_valid(i) == 1)
                k = k + 1;
                
                % Plane Count
                j = floor((i-1)/num_p) + 1;

                xyu = [xu(i); yu(i)];
                xyd = [xd(i); yd(i)];

                % The derivative of a undistorted normalized point
                [dxyu_dw, dxyu_dt] = Compute_dxxu(XYZw(:,i), R(:,(j-1)*3+1:j*3), t(:,j), w(:,j));

                % The derivative of a distorted normalized point
                [dxyd_dw, dxyd_dt, dxyd_dd] = Compute_dxyd(xyu, d, dxyu_dw, dxyu_dt);

                % The derivative of a distotred 2D pixel points
                [dxy_dk, dxy_dw, dxy_dt, dxy_dd] = ...
                    Compute_dxy(xyd, K, d, dxyd_dw, dxyd_dt, dxyd_dd, K_flag, d_flag);

                dxy_dw_all = zeros(2,3*m);
                dxy_dt_all = zeros(2,3*m);
                dxy_dw_all(:,(j-1)*3+1:j*3) = dxy_dw;
                dxy_dt_all(:,(j-1)*3+1:j*3) = dxy_dt;

                % Jacobian
                J(2*k-1:2*k, 1:num_Param) = [dxy_dk dxy_dw_all dxy_dt_all dxy_dd];
            end
        end
        
        
        % Compute the approximated Hessian matrix
        H = J'*J;
        
        if (n == 1)
            % Some of the diagonal elements of Hessian matrix 'H'
            % can be very large, and some can be very small. In this case,
            % the average of the diagonal elements of 'H' produce still
            % large lambda. If the large lambda is added to the diagonal
            % elements of 'H', it acts like a gradient descent and
            % the inverse of 'H' has very small diagonal elements.
            % Eventually it affects little to delta computation from the
            % first iteration, so it can cause very slow convergence
            % or no good results.
            % In this implementaion, the minimum value of the
            % diagonal elements of 'H' is used for lambda, so it acts
            % like a Gauss-Newton from the first iteration.
            
            % The following lambda value is used 
            % instead of 'lambda = 0.001*trace(H)/num_Param'.
            lambda = 0.001*min(diag(H));
          
            if (lambda < eps)
                lambda = eps;
            end
        else
             lambda = lambda/10;
             
             if (lambda < eps)
                lambda = eps;
            end
        end
    else
        lambda = lambda*10;
        
        % Display the Re-projection Error
        disp(['Iteration: ', num2str(n), ', Re-projection Error: ', num2str(rperr)]);
    end
    
    
    % Apply the damping factor to the Hessian matrix
    H_lm = H + (lambda * eye(num_Param, num_Param));

    % Prevent the matrix from being singular
    if (rcond(H_lm) < eps)
        lambda = lambda*10;
        H_lm = H + (lambda * eye(num_Param, num_Param));
    end
    
    % Compute the updated parameters
    delta = -inv(H_lm)*(J'*dist(:));
end


%% Output
varargout = {K, d, R, t, rperr, noIter};


%% Sub Functions
function [K_lm] = Update_K(K, kc, delta)
    if isequal(kc, [0 0 0 0 0])
        K_lm(1,1) = K(1,1);  % fx
        K_lm(2,2) = K(2,2);  % fy
        K_lm(1,3) = K(1,3);  % x0
        K_lm(2,3) = K(2,3);  % y0
        K_lm(1,2) = K(1,2);  % skew
        K_lm(3,3) = 1;
    elseif isequal(kc, [1 0 0 0 0])
        K_lm(1,1) = K(1,1) + delta(1);  % fx
        K_lm(2,2) = K(2,2) + delta(1);  % fy
        K_lm(1,3) = K(1,3);             % x0
        K_lm(2,3) = K(2,3);             % y0
        K_lm(1,2) = 0;                  % skew
        K_lm(3,3) = 1;
    elseif isequal(kc, [1 1 0 0 0])
        K_lm(1,1) = K(1,1) + delta(1);  % fx
        K_lm(2,2) = K(2,2) + delta(2);  % fy
        K_lm(1,3) = K(1,3);             % x0
        K_lm(2,3) = K(2,3);             % y0
        K_lm(1,2) = 0;                  % skew
        K_lm(3,3) = 1;
    elseif isequal(kc, [1 0 1 1 0])
        K_lm(1,1) = K(1,1) + delta(1);  % fx
        K_lm(2,2) = K(2,2) + delta(1);  % fy
        K_lm(1,3) = K(1,3) + delta(2);  % x0
        K_lm(2,3) = K(2,3) + delta(3);  % y0
        K_lm(1,2) = 0;                  % skew
        K_lm(3,3) = 1;
    elseif isequal(kc, [1 1 1 1 0])
        K_lm(1,1) = K(1,1) + delta(1);  % fx
        K_lm(2,2) = K(2,2) + delta(2);  % fy
        K_lm(1,3) = K(1,3) + delta(3);  % x0
        K_lm(2,3) = K(2,3) + delta(4);  % y0
        K_lm(1,2) = 0;                  % skew
        K_lm(3,3) = 1;
    elseif isequal(kc, [1 1 1 1 1])
        K_lm(1,1) = K(1,1) + delta(1);  % fx
        K_lm(2,2) = K(2,2) + delta(2);  % fy
        K_lm(1,3) = K(1,3) + delta(3);  % x0
        K_lm(2,3) = K(2,3) + delta(4);  % y0
        K_lm(1,2) = K(1,2) + delta(5);  % skew
        K_lm(3,3) = 1;
    else
        error('The constraint of the camera is wrong.');
    end
    
function [k1, k2, k3, p1, p2, d_lm] = Update_d(d, d_flag, d_delta)
    % Lens Distortion Parameters
    if (~isempty(d))
        if (d_flag == 1)
            if (length(d) == 1) 
                k1 = d(1) + d_delta(1);
                k2 = 0; k3 = 0; p1 = 0; p2 = 0;
                d_lm = k1;
            elseif (length(d) == 2) 
                k1 = d(1) + d_delta(1);
                k2 = d(2) + d_delta(2); 
                k3 = 0; p1 = 0; p2 = 0;
                d_lm = [k1; k2];
            elseif (length(d) == 3) 
                k1 = d(1) + d_delta(1);
                k2 = d(2) + d_delta(2); 
                k3 = d(3) + d_delta(3);
                p1 = 0; p2 = 0;
                d_lm = [k1; k2; k3];
            elseif (length(d) == 4) 
                k1 = d(1) + d_delta(1);
                k2 = d(2) + d_delta(2); 
                k3 = 0;
                p1 = d(3) + d_delta(3);
                p2 = d(4) + d_delta(4);
                d_lm = [k1; k2; p1; p2];
            elseif (length(d) == 5) 
                k1 = d(1) + d_delta(1);
                k2 = d(2) + d_delta(2); 
                k3 = d(3) + d_delta(3);
                p1 = d(4) + d_delta(4);
                p2 = d(5) + d_delta(5);
                d_lm = [k1; k2; k3; p1; p2];
            end
        else
            if (length(d) == 1) 
                k1 = d(1);
                k2 = 0; k3 = 0; p1 = 0; p2 = 0;
            elseif (length(d) == 2) 
                k1 = d(1);
                k2 = d(2);
                k3 = 0; p1 = 0; p2 = 0;
            elseif (length(d) == 3)
                k1 = d(1);
                k2 = d(2);
                k3 = d(3);
                p1 = 0; p2 = 0;
            elseif (length(d) == 4)
                k1 = d(1);
                k2 = d(2);
                k3 = 0;
                p1 = d(3);
                p2 = d(4);
            elseif (length(d) == 5) 
                k1 = d(1);
                k2 = d(2);
                k3 = d(3);
                p1 = d(4);
                p2 = d(5);
            end
            
            d_lm = d;
        end
    else
        k1 = 0;
        k2 = 0;
        k3 = 0;
        p1 = 0;
        p2 = 0;
        
        d_lm = [];
    end
        
function [dxyu_dw, dxyu_dt] = Compute_dxxu(XYZw, R, t, w)
    XYZc = R*XYZw + t;

    dxyu_dXYZc = Compute_dxyu_dXYZc(XYZc);
    dXYZc_dr   = Compute_dXYZc_dr(XYZw);
    dr_dw     = Compute_dr_dw(w);
    dXYZc_dt   = eye(3,3);

    dxyu_dw = dxyu_dXYZc*dXYZc_dr*dr_dw;
    dxyu_dt = dxyu_dXYZc*dXYZc_dt;
    

function [dxyu_dXYZc] = Compute_dxyu_dXYZc(XYZc)
Xc = XYZc(1);
Yc = XYZc(2);
Zc = XYZc(3);

dxyu_dXYZc = zeros(2,3);
dxyu_dXYZc(1,1) = 1/Zc;  
dxyu_dXYZc(2,1) = 0; 
dxyu_dXYZc(1,2) = 0;
dxyu_dXYZc(2,2) = 1/Zc;
dxyu_dXYZc(1,3) = -Xc/Zc^2;
dxyu_dXYZc(2,3) = -Yc/Zc^2;
    
    
function [dXYZc_dr] = Compute_dXYZc_dr(XYZw)
dXYZc_dr = zeros(3,9);
dXYZc_dr(1,1:3:end) = XYZw;
dXYZc_dr(2,2:3:end) = XYZw;
dXYZc_dr(3,3:3:end) = XYZw;

    
function [dr_dw] = Compute_dr_dw(w)
% r = [r11 r21 r31 r12 r22 r32 r13 r23 r33]

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

    dsth_dw   = [wx*cos(theta)/theta wy*cos(theta)/theta wz*cos(theta)/theta];   % d(sin(th))/d(w)
    domcth_dw = [wx*sin(theta)/theta, wy*sin(theta)/theta, wz*sin(theta)/theta]; % d(1-cos(th))/d(w)

    dwxh_dw = [1/theta-wx^2/theta^3, -wx*wy/theta^3, -wx*wz/theta^3];
    dwyh_dw = [-wx*wy/theta^3, 1/theta-wy^2/theta^3, -wy*wz/theta^3];
    dwzh_dw = [-wx*wz/theta^3, -wy*wz/theta^3, 1/theta-wz^2/theta^3];

    dwxh2pwyh2_dw = [2*wx/theta^2-2*(wx^3+wy^2*wx)/theta^4, 2*wy/theta^2-2*(wy^3+wx^2*wy)/theta^4, -2*(wx^2*wz+wy^2*wz)/theta^4]; %d(wxh^2 + wyh^2)/d(w)
    dwxh2pwzh2_dw = [2*wx/theta^2-2*(wx^3+wz^2*wx)/theta^4, -2*(wx^2*wy+wz^2*wy)/theta^4, 2*wz/theta^2-2*(wz^3+wx^2*wz)/theta^4]; %d(wxh^2 + wzh^2)/d(w)
    dwyh2pwzh2_dw = [-2*(wy^2*wx+wz^2*wx)/theta^4, 2*wy/theta^2-2*(wy^3+wz^2*wy)/theta^4, 2*wz/theta^2-2*(wz^3+wy^2*wz)/theta^4]; %d(wyh^2 + wzh^2)/d(w)

    wxh2pwyh2 = wxh^2 + wyh^2;
    wxh2pwzh2 = wxh^2 + wzh^2;
    wyh2pwzh2 = wyh^2 + wzh^2;

    omcth = 1-cos(theta);
    sth   = sin(theta);

    dr_dw = [-omcth*dwyh2pwzh2_dw - wyh2pwzh2*domcth_dw;                                             % dr11_dw
              sth*dwzh_dw + wzh*dsth_dw + wyh*omcth*dwxh_dw + wxh*omcth*dwyh_dw + wxh*wyh*domcth_dw; % dr21_dw
             -sth*dwyh_dw - wyh*dsth_dw + wzh*omcth*dwxh_dw + wxh*omcth*dwzh_dw + wxh*wzh*domcth_dw; % dr31_dw
             -sth*dwzh_dw - wzh*dsth_dw + wyh*omcth*dwxh_dw + wxh*omcth*dwyh_dw + wxh*wyh*domcth_dw; % dr12_dw
             -omcth*dwxh2pwzh2_dw - wxh2pwzh2*domcth_dw;                                             % dr22_dw
              sth*dwxh_dw + wxh*dsth_dw + wzh*omcth*dwyh_dw + wyh*omcth*dwzh_dw + wyh*wzh*domcth_dw; % dr32_dw
              sth*dwyh_dw + wyh*dsth_dw + wzh*omcth*dwxh_dw + wxh*omcth*dwzh_dw + wxh*wzh*domcth_dw; % dr13_dw
             -sth*dwxh_dw - wxh*dsth_dw + wzh*omcth*dwyh_dw + wyh*omcth*dwzh_dw + wyh*wzh*domcth_dw; % dr23_dw
             -omcth*dwxh2pwyh2_dw - wxh2pwyh2*domcth_dw];                                            % dr33_dw
end

           
function [dxyd_dw, dxyd_dt, dxyd_dd] = Compute_dxyd(xyu, d, dxyu_dw, dxyu_dt)
xu = xyu(1);
yu = xyu(2);
r  = sqrt(xu^2 + yu^2);

dxu_dw = dxyu_dw(1,:);
dyu_dw = dxyu_dw(2,:);

dxu_dt = dxyu_dt(1,:);
dyu_dt = dxyu_dt(2,:);


% The derivative of a radial distortion function
% rc = 1 + k1*r^2 + k2*r^4 + k3*r^6;
if (length(d) >= 1)
    % First radial distortion coefficients, k1
    k1 = d(1);
    rc = 1 + k1*r^2;

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
     rc = rc + k2*r^4;

     dr4_dw = 2*r^2*dr2_dw;
     dr4_dt = 2*r^2*dr2_dt;
     dc_dk2 = r^4;

     dc_dw = dc_dw + k2*dr4_dw;
     dc_dt = dc_dt + k2*dr4_dt;
     dc_dk = [dc_dk dc_dk2];
 end

 if ((length(d) == 3) || (length(d) == 5))
     % Third radial distortion coefficients, k3
     k3 = d(3);
     rc = rc + k3*r^6;

     dr6_dw = 3*r^4*dr2_dw;
     dr6_dt = 3*r^4*dr2_dt;
     dc_dk3 = r^6;

     dc_dw = dc_dw + k3*dr6_dw;
     dc_dt = dc_dt + k3*dr6_dt;
     dc_dk = [dc_dk dc_dk3];
 end

 if (~isempty(d))
     % The derivative of a radially distorted normalized point
     dxr_dw = rc*dxu_dw + xu*dc_dw;
     dyr_dw = rc*dyu_dw + yu*dc_dw;

     dxr_dt = rc*dxu_dt + xu*dc_dt;
     dyr_dt = rc*dyu_dt + yu*dc_dt;

     dxr_dk = xu*dc_dk;
     dyr_dk = yu*dc_dk;
 end


 % The derivative of a tangentially distorted normalized point
 if (length(d) >= 4)
     if (length(d) == 4)
         p1 = d(3);
         p2 = d(4);
     else
         p1 = d(4);
         p2 = d(5);
     end
     
     dxt_dxxu = [2*p1*yu+6*p2*xu 2*p1*xu+2*p2*yu];
     dyt_dxxu = [2*p1*xu+2*p2*yu 6*p1*yu+2*p2*xu];

     dxyu_dw = [dxu_dw; dyu_dw];
     dxyu_dt = [dxu_dt; dyu_dt];

     dxt_dw = dxt_dxxu*dxyu_dw;
     dyt_dw = dyt_dxxu*dxyu_dw;

     dxt_dt = dxt_dxxu*dxyu_dt;
     dyt_dt = dyt_dxxu*dxyu_dt;

     dxt_dp = [    2*xu*yu r^2+2*xu^2 ];
     dyt_dp = [r^2+2*yu^2      2*xu*yu];
 end


% The derivative of a distorted normalized point
if (isempty(d))
    dxd_dw = dxu_dw;
    dyd_dw = dyu_dw;

    dxd_dt = dxu_dt;
    dyd_dt = dyu_dt;
    
    dxd_dd = [];
    dyd_dd = [];

elseif (length(d) <= 3)
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

    dxd_dd = [dxr_dk dxt_dp];
    dyd_dd = [dyr_dk dyt_dp];
end

dxyd_dw = [dxd_dw; dyd_dw];
dxyd_dt = [dxd_dt; dyd_dt];
dxyd_dd = [dxd_dd; dyd_dd];
    
    
function [dxy_dk, dxy_dw, dxy_dt, dxy_dd] = Compute_dxy(xyd, K, d, dxyd_dw, dxyd_dt, dxyd_dd, K_flag, d_flag)
    % x = fx*xd + s*yd + x0
    % y = fy*yd + y0
    
    xd = xyd(1);
    yd = xyd(2);

    fx = K(1,1);
    fy = K(2,2);
    s  = K(1,2);

    dxd_dw = dxyd_dw(1,:);
    dyd_dw = dxyd_dw(2,:);

    dxd_dt = dxyd_dt(1,:);
    dyd_dt = dxyd_dt(2,:);

    dx_dfx = xd; dx_dfy =  0;
    dx_dx0 =  1; dx_dy0 =  0;
    dx_ds  = yd;
    dy_dfx =  0; dy_dfy = yd;
    dy_dx0 =  0; dy_dy0 =  1; 
    dy_ds  =  0;

    if (~isempty(d)) && (d_flag == 1)
        dxd_dd = dxyd_dd(1,:);
        dyd_dd = dxyd_dd(2,:);

        dx_dd = fx*dxd_dd + s*dyd_dd;
        dy_dd = fy*dyd_dd;
    end

    if isequal(K_flag, [0 0 0 0 0])
        dx_dk = [];
        dy_dk = [];
    elseif (isequal(K_flag, [1 0 0 0 0]) || isequal(K_flag, [0 1 0 0 0]))
        dx_dk = dx_dfx;
        dy_dk = dy_dfy;
    elseif isequal(K_flag, [1 1 0 0 0])
        dx_dk = [dx_dfx dx_dfy];
        dy_dk = [dy_dfx dy_dfy];
    elseif isequal(K_flag, [1 0 1 1 0])
        dx_dk = [dx_dfx dx_dx0 dx_dy0];
        dy_dk = [dy_dfy dy_dx0 dy_dy0];
    elseif isequal(K_flag, [1 1 1 1 0])
        dx_dk = [dx_dfx dx_dfy dx_dx0 dx_dy0];
        dy_dk = [dy_dfx dy_dfy dy_dx0 dy_dy0];
    elseif isequal(K_flag, [1 1 1 1 1])
        dx_dk = [dx_dfx dx_dfy dx_dx0 dx_dy0 dx_ds];
        dy_dk = [dy_dfx dy_dfy dy_dx0 dy_dy0 dy_ds];
    end

    dx_dw = fx*dxd_dw + s*dyd_dw;
    dy_dw = fy*dyd_dw;

    dx_dt = fx*dxd_dt + s*dyd_dt;
    dy_dt = fy*dyd_dt;

    dxy_dk = [dx_dk; dy_dk];
    dxy_dw = [dx_dw; dy_dw];
    dxy_dt = [dx_dt; dy_dt];

    if (isempty(d)) || (d_flag == 0)
        dxy_dd = [];
    else
        dxy_dd = [dx_dd; dy_dd];
    end