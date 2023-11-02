function AntData = ChangeGrid( AntData , Azi_grid , Ele_grid )
%% Set grid
        Ele = repmat( Ele_grid' , 1 , numel(Azi_grid) );
        Azi = repmat( Azi_grid , numel(Ele_grid) , 1 );        
        EleIDX = Ele_grid <= max(AntData.elevation_grid) & Ele_grid >= min(AntData.elevation_grid);
        AziIDX = Azi_grid <= max(AntData.azimuth_grid) & Azi_grid >= min(AntData.azimuth_grid);
        [V, H, T] = LinearInterpolateGrid( AntData , Azi(EleIDX,AziIDX) , Ele(EleIDX,AziIDX) );
%% Output
        AntData.elevation_grid = Ele_grid(EleIDX);
        AntData.azimuth_grid = Azi_grid(AziIDX);
        AntData.V = V;
        AntData.H = H;
        AntData.T = T;
        AntData.no_el = length(AntData.elevation_grid);
        AntData.no_az = length(AntData.azimuth_grid);
end
