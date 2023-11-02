function [ Coeff ] = GetCC_LOS( AntData, TxRot, RxRot )
%% Grid
ele = AntData.elevation_range;
azi = AntData.azimuth_range;
% Change grid with interpolate (Rough)
elevation_grid = ele(1):(10/180*pi):ele(2);
azimuth_grid = azi(1):(10/180*pi):azi(2);
AntData.tx = ChangeGrid( AntData.tx , azimuth_grid , elevation_grid );
AntData.rx = ChangeGrid( AntData.rx , azimuth_grid , elevation_grid );

% Change grid with interpolate (Need)
elevation_grid = ele(1):(1/180*pi):ele(2);
azimuth_grid = azi(1):(1/180*pi):azi(2);
AntData.tx = ChangeGrid( AntData.tx , azimuth_grid , elevation_grid );
AntData.rx = ChangeGrid( AntData.rx , azimuth_grid , elevation_grid );

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

%% Calculate channel coefficient
Coeff = CalculateCC( AntData, AntData.tx, AntData.rx );

end

