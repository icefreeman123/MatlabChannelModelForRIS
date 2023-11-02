function [ AntData ] = LoadParameter( file, TxLoc, RxLoc, RISLoc )
%% Setup
TxPower_dBm = 0; % unit: dBm.
AntData.TxPower_dBm = TxPower_dBm;
AntData.TxPower = 1e-3 * 10^(TxPower_dBm/10); % unit: W.
AntData.LossType = 'withloss'; % 'withloss' / 'noloss'
AntData.freq = 3.5e9; % Center Frequency (unit: Hz)

elevation_range = [0, pi];
azimuth_range = [0, 2*pi];
AntData.elevation_range = elevation_range;
AntData.azimuth_range = azimuth_range;

if exist('RISLoc','var')
        AntNum = 3;
else
        AntNum = 2;
end
%% Load radiation pattern
startIDX = 2;
path_now = pwd;
path_data = [path_now,'\', file];
cd([path_data]);
listDateTMP = dir;
listDate = cell((length(listDateTMP)-2), 2);
for i=1:size(listDate,1)
        listDate(i,1) = {listDateTMP(i+2).name}; 
end
for i=1:size(listDate,1)
        filename = listDate{i,1};
        ANTcsvTMP = xlsread([filename]);
        ANTcsv = ANTcsvTMP(startIDX:end , startIDX:end);
        listDate(i,2) = {ANTcsv}; clear ANTcsv
end
cd([path_now]);

%% Digital Phase Shifter (DPS)
phase = [-33.64, -61.63, -102.0, -142.1, ...
            -178.0, -225.7, -311.1, -365.5, ...
            -395.4, -423.8, -462.1, -508.4, ...
            -545.1, -596.2, -681.1, -733.9];
loss = [3.601, 3.516, 5.175, 7.309, ...
            6.442, 7.578, 8.294, 6.973, ...
            3.289, 3.336, 5.031, 7.238, ...
            6.305, 7.194, 7.500, 6.228];
AntData.DPS.phase = phase;
AntData.DPS.loss = loss;

%%
typeAll = {'tx','rx','ris'};
AntNeedAll = [1, 3, 2];

for i = 1:AntNum
        type = typeAll{1,i};
        AntNeed = AntNeedAll(1,i);
        
        CustAntV = listDate(AntNeed + (size(listDate,1)/2),:);
        CustAntH = listDate(AntNeed,:);
        Vtmp = CustAntV{1,2};
        Htmp = CustAntH{1,2};
        Vtmp = sqrt(10.^(Vtmp/10));
        Htmp = sqrt(10.^(Htmp/10));
        no_elNeed = size(Vtmp, 1);
        no_azNeed = size(Vtmp, 2);        
        elevation_gridNeed = linspace(elevation_range(1), elevation_range(2), no_elNeed);
        azimuth_gridNeed = linspace(azimuth_range(1), azimuth_range(2), no_azNeed);

        switch type
                case 'tx'
                        AntData.tx.name = 'patch';
                        AntData.tx.V = {Vtmp};
                        AntData.tx.H = {Htmp};
                        AntData.tx.T = {sqrt(Vtmp.^2 + Htmp.^2)};
                        AntData.tx.position = {TxLoc}; 
                        AntData.tx.no_ant = 1;
                        AntData.tx.no_el = no_elNeed;
                        AntData.tx.no_az = no_azNeed;
                        AntData.tx.elevation_grid = elevation_gridNeed;
                        AntData.tx.azimuth_grid = azimuth_gridNeed;    
                        AntData.tx.normal = [0;0;1]; % Normal vector of Tx array
                        AntData.tx.orient = [0;1;0]; % Orientation of Tx
                case 'rx'
                        AntData.rx.name = 'patch';
                        AntData.rx.V = {Vtmp};
                        AntData.rx.H = {Htmp};
                        AntData.rx.T = {sqrt(Vtmp.^2 + Htmp.^2)};
                        AntData.rx.position = {RxLoc};
                        AntData.rx.no_ant = 1;
                        AntData.rx.no_el = no_elNeed;
                        AntData.rx.no_az = no_azNeed;
                        AntData.rx.elevation_grid = elevation_gridNeed;
                        AntData.rx.azimuth_grid = azimuth_gridNeed;    
                        AntData.rx.normal = [0;0;1]; % Normal vector of Rx array
                        AntData.rx.orient = [0;1;0]; % Orientation of Rx
                case 'ris'
                        AntData.ris.name = 'patch';
                        AntData.ris.V = {Vtmp};
                        AntData.ris.H = {Htmp};
                        AntData.ris.T = {sqrt(Vtmp.^2 + Htmp.^2)};
                        AntData.ris.position = {RISLoc};
                        AntData.ris.no_ant = 1;
                        AntData.ris.no_el = no_elNeed;
                        AntData.ris.no_az = no_azNeed;
                        AntData.ris.elevation_grid = elevation_gridNeed;
                        AntData.ris.azimuth_grid = azimuth_gridNeed;                        
                        AntData.ris.area_a = 0.094; % unit: m.
                        AntData.ris.area_b = 0.073;
                        AntData.ris.normal = [0;0;1]; % Normal vector of RIS array
                        AntData.ris.orient = [0;1;0]; % Orientation of RIS
        end
end
%%
TxRot = {[-10, 0, 0], [0, 10, 0]}; % Rotate in order. Unit: degree
RxRot = {[-10, 0, 0], [0, 10, 0]}; % Rotate in order. Unit: degree
RISRot = {[-10, 0, 0], [0, 10, 0]}; % Rotate in order. Unit: degree
TxRotIDX_num = length(TxRot);
for i = 1:TxRotIDX_num
        TxRotNeed = TxRot{i};
        AntData.tx = TuneRotateAntenna( AntData.tx, TxRotNeed );
end
RxRotIDX_num = length(RxRot);
for i = 1:RxRotIDX_num
        RxRotNeed = RxRot{i};
        AntData.rx = TuneRotateAntenna( AntData.rx, RxRotNeed );
end
RISRotIDX_num = length(RISRot);
for i = 1:RISRotIDX_num
        RISRotNeed = RISRot{i};
        AntData.ris = TuneRotateAntenna( AntData.ris, RISRotNeed );
end

end
function [ AntData ] = TuneRotateAntenna( AntData, rot_angle )
%% Setup
no_ant = AntData.no_ant;

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
Mo_21 = -1*Mo_12;
Mo_22 = Mo_11; 

%% Reading the original pattern at the rotated signal direction
[ V_org, H_org, T_org ] = LinearInterpolateGrid( AntData , phi_mat_new , theta_mat_new );

%% Transformation of the polarization        
Vtmp = V_org(:);
Htmp = H_org(:);
Vcell = cell(no_ant,1);
Hcell = cell(no_ant,1);
Tcell = cell(no_ant,1);
for i = 1:no_ant
        Vtmp_i = Vtmp{i};
        Htmp_i = Htmp{i};
        V_new = Mo_11.*Vtmp_i + Mo_12.*Htmp_i;
        H_new = Mo_21.*Vtmp_i + Mo_22.*Htmp_i;
        T_new = sqrt(V_new.^2 + H_new.^2);
        Vcell(i) = {V_new};        
        Hcell(i) = {H_new};        
        Tcell(i) = {T_new};
end
Vcell = reshape(Vcell, size(V_org));
Hcell = reshape(Hcell, size(H_org));
Tcell = reshape(Tcell, size(T_org));

%% Output
AntData.V = Vcell;
AntData.H = Hcell; 
AntData.T = Tcell;
end


