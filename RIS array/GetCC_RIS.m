function [ Coeff ] = GetCC_RIS( Coeff_LOS, AntData, RISLoc, TxRot, RxRot, RISRot )
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
RISRotIDX_num = length(RISRot);
for i = 1:RISRotIDX_num
        RISRotNeed = RISRot{i};
        AntData.ris = RotateAntenna( AntData.ris, RISRotNeed );
end

%% Calculate channel coefficient
% Calculate Scatter coefficient on RIS
ScatterEff = 11.2;
Coeff_scatter = CalculateScatter( AntData, RISLoc );
Coeff_scatter = ScatterEff * Coeff_scatter;

% Calculate Reflection coefficient on RIS
RefEff = 8;
Coeff_reflect_RIS = CalculateReflect( AntData, RISLoc );
Coeff_reflect_RIS = RefEff * Coeff_reflect_RIS;

% Calculate RIS channel coefficient
Coeff_Tx2RIS = CalculateCC( AntData, AntData.tx, AntData.ris );
Coeff_RIS2Rx = CalculateCC( AntData, AntData.ris, AntData.rx );

% Random initial phase of RIS
ris_initial = exp(1j*(2*pi)*repmat(5/12,AntData.ris.no_ant, 1));
ris_initial_diag = diag(ris_initial(:));
AntData.ris.initial_diag = ris_initial_diag;

% Control DPS
switch AntData.CtlMethod
        case 'Perfect_beam' % Without random initial phase of RIS
                Coeff = ControlRIS_Perfectbeam( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
        case 'DPS_beam' % Without random initial phase of RIS
                Coeff = ControlRIS_DPSbeam( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
        case 'BGA'
                Coeff = ControlRIS_BGA( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
        case 'BGA_polar'
                Coeff = ControlRIS_BGApolar( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
end
disp([ 'Controlling algorithm: ', AntData.CtlMethod ]);
end

