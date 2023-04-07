function [ AntData ] = RotateAntenna( AntData, rot_angle )
%% Setup
% Get the angles.
phi = AntData.azimuth_grid;
theta = AntData.elevation_grid';
no_az = AntData.no_az;
no_el = AntData.no_el;

% Change the unit from degree to rad.
rot_angle = rot_angle/180*pi;

%% Rotation matrix
% Calculate the rotation matrix.
rx = rot_angle(1);
ry = rot_angle(2);
rz = rot_angle(3);
rot = zeros(3,3);
rot(1,1) = cos(rz) * cos(ry);
rot(2,1) = sin(rz) * cos(ry);
rot(3,1) = -1 * sin(ry);
rot(1,2) = cos(rz) * sin(ry) * sin(rx) - sin(rz) * cos(rx);
rot(2,2) = sin(rz) * sin(ry) * sin(rx) + cos(rz) * cos(rx);
rot(3,2) = cos(ry) * sin(rx);
rot(1,3) = cos(rz) * sin(ry) * cos(rx) + sin(rz) * sin(rx);
rot(2,3) = sin(rz) * sin(ry) * cos(rx) - cos(rz) * sin(rx);
rot(3,3) = cos(ry) * cos(rx);

%% Theta and Phi of signal direction in the global and local coordinate system
% Change Theta from [0, pi] to [pi/2, -pi/2].
theta_mat = repmat(theta, [1, length(phi)]);
theta_mat = -1*theta_mat + pi/2;
% Change Phi from [0, 2*pi] to [pi, -pi].
phiTMP = phi;
phiTMP_idx = phiTMP >= pi; %%%
phiTMP(phiTMP_idx) = phiTMP(phiTMP_idx) - 2*pi;
phi_mat = repmat(phiTMP, [length(theta), 1]);
% Transform from spherical coordinates into Cartesian coordinates.
[ Fx, Fy, Fz ] = sph2cart(phi_mat, theta_mat, 1e10);

% Rotation of signal direction
Fx_r = rot(1,1)*Fx + rot(2,1)*Fy +rot(3,1)*Fz; % Transpose
Fy_r = rot(1,2)*Fx + rot(2,2)*Fy +rot(3,2)*Fz;
Fz_r = rot(1,3)*Fx + rot(2,3)*Fy +rot(3,3)*Fz;

% Transforme from Cartesian coordinates into spherical coordinates.
[phi_mat_new, theta_mat_new] = cart2sph( Fx_r, Fy_r, Fz_r ); % [pi, -pi] and [pi/2, -pi/2].
% The pattern is not defined when the theta > 90 degree and < -90 degree.
err_limit = 1e-5;
theta_out_idx = theta_mat_new < -pi/2 + err_limit;
theta_mat_new(theta_out_idx)  = -pi/2 + err_limit;
theta_out_idx = theta_mat_new > pi/2 - err_limit;
theta_mat_new(theta_out_idx)  = pi/2 - err_limit;
% Change Theta from [pi/2, -pi/2] to [0, pi].
theta_mat_new = -1*theta_mat_new + pi/2;
% Change Phi from [pi, -pi] to [0, 2*pi].
phi_mat_new_pi_idx = phi_mat_new < 0;
phi_mat_new(phi_mat_new_pi_idx) = phi_mat_new(phi_mat_new_pi_idx) + 2*pi;

%% Transformation matrix
% Global coordinate system in theta direction.
T_theta2x_g = cos(theta) * cos(phi);
T_theta2y_g = cos(theta) * sin(phi);
T_theta2z_g = -1 * sin(theta) * ones(1,no_az);

% Global coordinate system in phi direction.
T_phi2x_g = ones(no_el,1) * -1 * sin(phi);
T_phi2y_g = ones(no_el,1) * cos(phi);
T_phi2z_g = zeros( no_el , no_az );

% Local coordinate system in theta direction.
T_theta2x_l = cos(theta_mat_new) .* cos(phi_mat_new);
T_theta2y_l = cos(theta_mat_new) .* sin(phi_mat_new);
T_theta2z_l = -1 * sin(theta_mat_new);

% Local coordinate system in phi direction.
T_phi2x_l = -1 * sin(phi_mat_new);
T_phi2y_l = cos(phi_mat_new);
T_phi2z_l = zeros( no_el , no_az );

%% Polarization rotation matrix
T_theta2x_tmp = rot(1,1)*T_theta2x_l + rot(1,2)*T_theta2y_l + rot(1,3)*T_theta2z_l;
T_theta2y_tmp = rot(2,1)*T_theta2x_l + rot(2,2)*T_theta2y_l + rot(2,3)*T_theta2z_l;
T_theta2z_tmp = rot(3,1)*T_theta2x_l + rot(3,2)*T_theta2y_l + rot(3,3)*T_theta2z_l;
T_phi2x_tmp = rot(1,1)*T_phi2x_l + rot(1,2)*T_phi2y_l + rot(1,3)*T_phi2z_l;
T_phi2y_tmp = rot(2,1)*T_phi2x_l + rot(2,2)*T_phi2y_l + rot(2,3)*T_phi2z_l;
T_phi2z_tmp = rot(3,1)*T_phi2x_l + rot(3,2)*T_phi2y_l + rot(3,3)*T_phi2z_l;
Mo_11 = T_theta2x_tmp.*T_theta2x_g + T_theta2y_tmp.*T_theta2y_g + T_theta2z_tmp.*T_theta2z_g;
Mo_12 = T_phi2x_tmp.*T_theta2x_g + T_phi2y_tmp.*T_theta2y_g + T_phi2z_tmp.*T_theta2z_g;
Mo_21 = T_theta2x_tmp.*T_phi2x_g + T_theta2y_tmp.*T_phi2y_g + T_theta2z_tmp.*T_phi2z_g;
Mo_22 = T_phi2x_tmp.*T_phi2x_g + T_phi2y_tmp.*T_phi2y_g + T_phi2z_tmp.*T_phi2z_g;

%% Reading the original pattern at the rotated signal direction
[ V_org, H_org ] = LinearInterpolate( AntData , phi_mat_new , theta_mat_new );

%% Transformation of the polarization
V_new = Mo_11.*V_org + Mo_12.*H_org;
H_new = Mo_21.*V_org + Mo_22.*H_org;
        
%% Output
AntData.V = V_new;
AntData.H = H_new; 
AntData.T = sqrt(V_new.^2 + H_new.^2);
end

