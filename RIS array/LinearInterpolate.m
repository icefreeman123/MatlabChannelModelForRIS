function [ V, H ] = LinearInterpolate( V, H, AntData, azimuth_meshgrid_new, elevation_meshgrid_new )
%% Linear Interpolation
        azimuth_grid  = AntData.azimuth_grid;
        elevation_grid = AntData.elevation_grid;
        azimuth_meshgrid = repmat(azimuth_grid, [length(elevation_grid),1]);
        elevation_meshgrid = repmat(elevation_grid.', [1,length(azimuth_grid)]);
        V_new = interp2(azimuth_meshgrid, elevation_meshgrid, V, azimuth_meshgrid_new, elevation_meshgrid_new, 'linear');
        V_new(:,end) = V_new(:,1);
        H_new = interp2(azimuth_meshgrid, elevation_meshgrid, H, azimuth_meshgrid_new, elevation_meshgrid_new, 'linear');
        H_new(:,end) = H_new(:,1);
%% Output
        V = V_new;
        H = H_new;
end

