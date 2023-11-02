% =========================================================================
% Article Title:
%     A Novel Channel Model for Reconfigurable Intelligent Surfaces 
%     with Consideration of Polarization and Switch Impairments
%     (RIS array in Fig.12(b))(The case of receiver with rotation)
% -------------------------------------------------------------------------
% Function:
%     CalculateCC.m / CalculateReflect.m / CalculateScatter.m / 
%     ChangeGrid.m / LinearInterpolate.m / LinearInterpolateGrid.m / 
%     LoadParameter.m / ArrayGenerate.m / 
%     RodriguesRotVec.m / RotateAntenna.m / 
%     ControlRIS_Perfectbeam.m / ControlRIS_DPSbeam.m / 
%     ControlRIS_BGA.m / ControlRIS_BGApolar.m /
%     GetCC_LOS.m / GetCC_RIS.m 
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
%% Setup
TxLoc = [ 0 ; 0 ; 0.8 ]; % Array center. Unit: m.
RxLoc = [ 0 ; 0.2 ; 1.6 ]; % Array center. Unit: m.
RISLoc = [ 0 ; 0 ; 0 ]; % Array center. Unit: m.
AntData = LoadParameter( 'Pattern', TxLoc, RxLoc, RISLoc ); % Single antenna data

TxRot = {[0, 90, 0], [0, 0, 90]}; % Rotate in order. Unit: degree
RxRot = {[0, 90, 0], [0, 0, 90]}; % Rotate in order. Unit: degree
RISRot = {[0, 90, 0], [180, 0, 0], [0, 0, -90]}; % Rotate in order. Unit: degree

AntData.tx.size = [1, 1]; % Array size
AntData.tx = ArrayGenerate( AntData, AntData.tx ); % Duplicate array data
AntData.rx.size = [1, 1]; % Array size
AntData.rx = ArrayGenerate( AntData, AntData.rx ); % Duplicate array data
AntData.ris.size = [4, 4]; % Array size
AntData.ris = ArrayGenerate( AntData, AntData.ris ); % Duplicate array data

%%% Controlling algorithm:
AntData.CtlMethod = 'Perfect_beam'; % Perfect_beam / DPS_beam / BGA / BGA_polar
% BGA_polar needs lots of run time to rotate antenna polarization. (about 1.6 hr)
%% Start
NumPt_Rx = 21;
Track_Rx = repmat(RxLoc,1,NumPt_Rx) + [0;0.1;0]*linspace(0,5,NumPt_Rx);
TrackTheta_Rx = linspace(0,2*pi,NumPt_Rx);
ShowIdxNeed = 2; % Y axis

%%% The case of rotared receiver 
CoeffR_LOS_stack = zeros(1,NumPt_Rx);
CoeffR_RIS_stack = zeros(1,NumPt_Rx);
for Rx_idx = 1:NumPt_Rx
        AntData.rx.position = { Track_Rx(:, Rx_idx) };
        disp([ 'Rx Location: ', num2str(Track_Rx(ShowIdxNeed, Rx_idx)) ]);
        RxRot_new = {[0, 0, TrackTheta_Rx(Rx_idx)/pi*180]};
        RxRotNew = [RxRot, RxRot_new];
        disp([ 'Rx Rotation: ', num2str(TrackTheta_Rx(Rx_idx)/pi*180) ]);
        %
        Coeff_LOS = GetCC_LOS( AntData, TxRot, RxRotNew );
        Coeff_RIS = GetCC_RIS( Coeff_LOS, AntData, RISLoc, TxRot, RxRotNew, RISRot ); % Consider the controlling algorithm
        %
        CoeffR_LOS_stack(1,Rx_idx) = Coeff_LOS;
        CoeffR_RIS_stack(1,Rx_idx) = Coeff_RIS;
end

%% VNA
CoeffR_stack = CoeffR_LOS_stack + CoeffR_RIS_stack;
PowerR_dB = 10*log10(AntData.TxPower * (abs(CoeffR_stack).^2) / (1e-3)) - AntData.TxPower_dBm;

%% Plot
Track = Track_Rx(ShowIdxNeed,:);
figure;
PR = plot(Track, PowerR_dB, '--or','LineWidth',1.5);
grid on; axis square;
xlabel('Location of Rx on the y axis (m)');
ylabel('Loss (dB)');
legend([PR], 'Total loss with rotation');
axis([Track(1,1) Track(1,end) -60 -10]);
set(gca,'XTick', Track(1,1):0.1:Track(1,end));
set(gca,'YTick', -60:10:-10);
set(gca,'YTicklabel', {'60','50','40','30','20','10'});
