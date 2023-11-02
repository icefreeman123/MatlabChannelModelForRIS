function [ V, H, T ] = LinearInterpolateGrid( AntData, azimuth_meshgrid_new, elevation_meshgrid_new )
%% Linear Interpolation
        azimuth_grid  = AntData.azimuth_grid;
        elevation_grid = AntData.elevation_grid;
        azimuth_meshgrid = repmat(azimuth_grid, [length(elevation_grid),1]);
        elevation_meshgrid = repmat(elevation_grid.', [1,length(azimuth_grid)]);
        V = AntData.V;
        H = AntData.H;
        Vtmp = V(:);
        Htmp = H(:);
        no_ant = AntData.no_ant;
        Vcell = cell(no_ant,1);
        Hcell = cell(no_ant,1);
        Tcell = cell(no_ant,1);
        for i = 1:no_ant
                Vtmp_i = Vtmp{i};
                V_new = interp2(azimuth_meshgrid, elevation_meshgrid, Vtmp_i, ...
                                azimuth_meshgrid_new, elevation_meshgrid_new, 'linear');
                V_new(:,end) = V_new(:,1);
                Vcell(i) = {V_new};
                Htmp_i = Htmp{i};
                H_new = interp2(azimuth_meshgrid, elevation_meshgrid, Htmp_i, ...
                                azimuth_meshgrid_new, elevation_meshgrid_new, 'linear');
                H_new(:,end) = H_new(:,1);
                Hcell(i) = {H_new};
                T_new = sqrt(V_new.^2 + H_new.^2);
                Tcell(i) = {T_new};
        end
        Vcell = reshape(Vcell, size(V));
        Hcell = reshape(Hcell, size(H));
        Tcell = reshape(Tcell, size(H));
%% Output
        V = Vcell;
        H = Hcell;
        T = Tcell;
end

