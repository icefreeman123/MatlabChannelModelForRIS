% =========================================================================
% Article Title:
%     A Novel Channel Model for Reconfigurable Intelligent Surfaces 
%     with Consideration of Polarization and Switch Impairments
% -------------------------------------------------------------------------
% Function:
%     CalculateCC.m / CalculateReflect.m / CalculateScatter.m / 
%     ChangeGrid.m / LinearInterpolate.m / 
%     LoadParameter.m / RotateAntenna.m
% -------------------------------------------------------------------------
% Data:
%     Pattern file (AntH2.xlsx / AntH3.xlsx / AntH4.xlsx / 
%                   AntV2.xlsx / AntV3.xlsx / AntV4.xlsx)
% -------------------------------------------------------------------------
% Copyright (c) 2023 De-Ming Chian and Chao-Kai Wen
% E-mail: icefreeman123@gmail.com and chaokai.wen@mail.nsysu.edu.tw
% Reference: QuaDriGa_2016.09.05_v1.4.8-571
% =========================================================================
close all; clear all; clc
%% Set up location & Load antenna gain
TxLoc = [ 0 ; -0.1 ; 0.1 ]; % unit: m
RxLoc = [ 0 ; 0.1 ; 0.1 ]; 
RISLoc = [ 0 ; 0 ; 0 ]; 
AntData = LoadParameter( 'Pattern', TxLoc, RxLoc, RISLoc );

%% Digital Phase Shifter (DPS)
DPS_phase = [-33.64, -61.63, -102.0, -142.1, ...
            -178.0, -225.7, -311.1, -365.5, ...
            -395.4, -423.8, -462.1, -508.4, ...
            -545.1, -596.2, -681.1, -733.9];
DPS_loss = [3.601, 3.516, 5.175, 7.309, ...
            6.442, 7.578, 8.294, 6.973, ...
            3.289, 3.336, 5.031, 7.238, ...
            6.305, 7.194, 7.500, 6.228];
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360) + 360;
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360) + 360;
DPS_phase = DPS_phase/180*pi;
DPS_need = sqrt(10.^(-0.1*DPS_loss)) .* exp(1j*DPS_phase + 1j*(2*pi)*5/12); 

%% Grid
ele = AntData.elevation_range;
azi = AntData.azimuth_range;
% Change grid with interpolate (Rough)
elevation_grid = ele(1):(10/180*pi):ele(2);
azimuth_grid = azi(1):(10/180*pi):azi(2);
AntData.tx = ChangeGrid( AntData.tx , azimuth_grid , elevation_grid );
AntData.rx = ChangeGrid( AntData.rx , azimuth_grid , elevation_grid );
AntData.ris = ChangeGrid( AntData.ris , azimuth_grid , elevation_grid );
% Change grid with interpolate (Need)
elevation_grid = ele(1):(1/180*pi):ele(2);
azimuth_grid = azi(1):(1/180*pi):azi(2);
AntData.tx = ChangeGrid( AntData.tx , azimuth_grid , elevation_grid );
AntData.rx = ChangeGrid( AntData.rx , azimuth_grid , elevation_grid );
AntData.ris = ChangeGrid( AntData.ris , azimuth_grid , elevation_grid );

%% Rotate pattern
% Rotate in order. unit: degree
TxRot = {[0, 90, 0], [0, 0, 90]}; 
RxRot = {[0, 90, 0], [0, 0, 90]};
RISRot_Rot0 = {[0, 90, 0], [180, 0, 0], [0, 0, 90]}; % Case 1 (without rotation)
RISRot_Rot90 = {[0, 90, 0], [180, 0, 0], [0, 0, 90], [0, 0, 90]}; % Case 2 (rotate 90 deg.)
TxRotIDX_num = length(TxRot);
for i = 1:TxRotIDX_num
        TxRotNeed = TxRot{i};
        AntData.tx = RotateAntenna( AntData.tx, TxRotNeed );
end
RxRotIDX_num = length(RxRot);
for i = 1:RxRotIDX_num
        RxRotNeed = RxRot{i};
        AntData.rx = RotateAntenna( AntData.rx, RxRotNeed );
end
AntData_Rot0 = AntData;
AntData_Rot90 = AntData;
RISRotIDX_num = length(RISRot_Rot0);
for i = 1:RISRotIDX_num
        RISRotNeed = RISRot_Rot0{i};
        AntData_Rot0.ris = RotateAntenna( AntData_Rot0.ris, RISRotNeed );
end
RISRotIDX_num = length(RISRot_Rot90);
for i = 1:RISRotIDX_num
        RISRotNeed = RISRot_Rot90{i};
        AntData_Rot90.ris = RotateAntenna( AntData_Rot90.ris, RISRotNeed );
end

%% Move (Y axis)
Dstack_Open_Rot0 = zeros(1, length(DPS_need));
Dstack_Open_Rot90 = zeros(1, length(DPS_need));
for j = 1:length(DPS_need)
        % Antenna position
        AntData_Rot0.tx.position = [ 0 ; -0.1 ; 0.25 ]; 
        AntData_Rot0.rx.position = [ 0 ; 0.1 ; 0.25 ];
        AntData_Rot90.tx.position = [ 0 ; -0.1 ; 0.25 ]; 
        AntData_Rot90.rx.position = [ 0 ; 0.1 ; 0.25 ];

        % Calculate LOS channel coefficient
        [ Coeff_Tx2Rx_Rot0, ~ ] = CalculateCC( AntData_Rot0, AntData_Rot0.tx, AntData_Rot0.rx );
        [ Coeff_Tx2Rx_Rot90, ~ ] = CalculateCC( AntData_Rot90, AntData_Rot90.tx, AntData_Rot90.rx );

        % Calculate RIS channel coefficient
        DPS = DPS_need(j);
        [ Coeff_Tx2RIS_Rot0, ~ ] = CalculateCC( AntData_Rot0, AntData_Rot0.tx, AntData_Rot0.ris );
        [ Coeff_RIS2Rx_Rot0, ~ ] = CalculateCC( AntData_Rot0, AntData_Rot0.ris, AntData_Rot0.rx );
        Coeff_Tx2RIS2Rx_Rot0 = Coeff_Tx2RIS_Rot0 * Coeff_RIS2Rx_Rot0 * DPS;
        [ Coeff_Tx2RIS_Rot90, ~ ] = CalculateCC( AntData_Rot90, AntData_Rot90.tx, AntData_Rot90.ris );
        [ Coeff_RIS2Rx_Rot90, ~ ] = CalculateCC( AntData_Rot90, AntData_Rot90.ris, AntData_Rot90.rx );
        Coeff_Tx2RIS2Rx_Rot90 = Coeff_Tx2RIS_Rot90 * Coeff_RIS2Rx_Rot90 * DPS;

        % Calculate Scatter coefficient on RIS
        [ Coeff_scatter_Rot0, ~ ] = CalculateScatter( AntData_Rot0 );
        [ Coeff_scatter_Rot90, ~ ] = CalculateScatter( AntData_Rot90 );

        % Calculate Reflection coefficient on RIS
        [ Coeff_reflect_Rot0, ~ ] = CalculateReflect( AntData_Rot0);
        [ Coeff_reflect_Rot90, ~ ] = CalculateReflect( AntData_Rot90);

        % Combine coefficient (Channel coefficient)
        Coeff_Open_Rot0 = Coeff_Tx2Rx_Rot0 + Coeff_reflect_Rot0 + Coeff_scatter_Rot0 + Coeff_Tx2RIS2Rx_Rot0;
        Coeff_Open_Rot90 = Coeff_Tx2Rx_Rot90 + Coeff_reflect_Rot90 + Coeff_scatter_Rot90 + Coeff_Tx2RIS2Rx_Rot90;

        % Stack data
        Dstack_Open_Rot0(1,j) = Coeff_Open_Rot0;
        Dstack_Open_Rot90(1,j) = Coeff_Open_Rot90;
end

%% Results
% Power Ratio
Power_stack = AntData.TxPower * (abs(Dstack_Open_Rot0).^2);
Power_stack_dBm = 10*log10(Power_stack / (1e-3)); % unit: dBm
Power_ratio_dB = Power_stack_dBm - AntData.TxPower_dBm; % Corresponding to the meas. of VNA
Power_ratio_dB_Open_Rot0 = Power_ratio_dB;
Power_stack = AntData.TxPower * (abs(Dstack_Open_Rot90).^2);
Power_stack_dBm = 10*log10(Power_stack / (1e-3)); % unit: dBm
Power_ratio_dB = Power_stack_dBm - AntData.TxPower_dBm; % Corresponding to the meas. of VNA
Power_ratio_dB_Open_Rot90 = Power_ratio_dB;

% Measured results
VNA_dB_Open_Rot0 = [-31.46 ; -30.53 ; -29.30 ; -28.68 ; ...
                -29.16 ; -30.31 ; -30.80 ; -31.96 ; ...
                -31.26 ; -29.57 ; -28.73 ; -28.85 ; ...
                -29.47 ; -29.39 ; -31.12 ; -31.77]; % Case 1 by VNA
VNA_dB_Open_Rot90 = [-30.33 ; -30.12 ; -30.21 ; -29.94 ; ...
                -30.07 ; -30.24 ; -30.08 ; -30.14 ; ...
                -30.23 ; -30.10 ; -30.49 ; -29.77 ; ...
                -30.19 ; -30.94 ; -30.47 ; -30.20]; % Case 2 by VNA

% Plot results
DPSstate_idx = 1:length(DPS_need);
figure;
a = plot(DPSstate_idx, Power_ratio_dB_Open_Rot0,'r-','LineWidth',1.5); hold on;
b = plot(DPSstate_idx, VNA_dB_Open_Rot0,'r--s','LineWidth',1.5,'MarkerSize',8); 
c = plot(DPSstate_idx, Power_ratio_dB_Open_Rot90,'b-','LineWidth',1.5);
d = plot(DPSstate_idx, VNA_dB_Open_Rot90,'b--x','LineWidth',1.5,'MarkerSize',8); 
hold off; grid on; axis square;
title('RIS element for the different polarization direction');
xlabel('DPS state');
ylabel('Loss (dB)');
legend([a,b,c,d], 'Rot. 0 deg / Sim.','Rot. 0 deg / Meas.', 'Rot. 90 deg / Sim.','Rot. 90 deg / Meas.'); 
axis([DPSstate_idx(1) DPSstate_idx(end) -40 -20]);
set(gca,'XTick', DPSstate_idx(1,1):3:DPSstate_idx(1,end));
set(gca,'YTick', -40:5:-20);
set(gca,'YTicklabel', {'40','35','30','25','20'});

