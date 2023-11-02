function [ Coeff, Results ] = CalculateCC( AntData, tx, rx )
%% Setup
        tx_position = tx.position;
        rx_position = rx.position;
        freq = AntData.freq;
        LossType = AntData.LossType;
        speed = physconst('lightspeed'); % speed of light(unit: m/s)
        lambda = speed/freq; % wavelength(unit: Hz)
        k = (2*pi)/lambda; % Wave number
%% Amplitude
        [AOD, AOA] = CalculateAngle( tx_position, rx_position ); % Angles between Tx and Rx
        aod_r = AOD(2).'*pi/180;
        aoa_r = AOA(2).'*pi/180;
        eod_r = AOD(1).'*pi/180;
        eoa_r = AOA(1).'*pi/180;
        [ Vt, Ht ] = LinearInterpolate( tx , aod_r , eod_r );
        [ Vr, Hr ] = LinearInterpolate( rx , aoa_r , eoa_r );
        V_final = Vr *Vt;
        H_final = Hr *Ht;
        Amp = V_final - H_final;
%% Phase        
        r_s = rx_position - tx_position;
        norm_r_s = sqrt( sum( r_s.^2 ,1 ) );
        psi_lms = k * mod(norm_r_s, lambda); 
%% Path Loss
        switch LossType
                case 'withloss'
                        loss_dB = pathloss( tx_position, rx_position, lambda );
                case 'noloss' % otherwise
                        loss_dB = 0;
        end
        loss = sqrt( 10.^( 0.1 * -1 * loss_dB ) );     
%% Calculate Coeff.
        Coeff = loss * Amp * exp(-1j*psi_lms);
%% Results
        Results.V_Amp = V_final;
        Results.H_Amp = H_final;
        Results.Amp = Amp;
        Results.Phase = psi_lms;
        Results.PathLoss = loss;
end
%% Pathloss
function [ loss ] = pathloss( tx_position, rx_position, lambda )
        txpos = tx_position;
        rxpos = rx_position;
        d_3d = sqrt(sum( (rxpos-txpos).^2 ));
        loss = 20 * log10(4*pi*d_3d / lambda);
end
%% Angle
function [ AOD, AOA ] = CalculateAngle( T_position, R_position )
        % Vector in Cart. coordinate
        X_vec = R_position(1) - T_position(1);
        Y_vec = R_position(2) - T_position(2);
        Z_vec = R_position(3) - T_position(3);        

        angles = zeros( 4, 1 ); 
        % AOD-Horizontal polarization
        angles(1,1) = atan2( Y_vec, X_vec ); % [¡Vpi, pi] 
        angles(1, isnan(angles(1,1)) ) = 0;
        % AOD-Vertical polarization
        angles(3,1) = atan2( Z_vec, sqrt(X_vec.^2 + Y_vec.^2) ); % [¡Vpi, pi] 
        angles(3, isnan(angles(3,1)) ) = 0;        
        % AOA-Horizontal polarization
        angles(2,1) = pi + angles(1,1);
        % AOA-Vertical polarization
        angles(4,1) = -1 * angles(3,1);
        
        % Change to (0,360) deg.
        angles = angles * 180/pi;
        angles(1,1) = mod( angles(1,1) , 360);
        angles(2,1) = mod( angles(2,1) , 360);
        
        % AOD (1. V / 2. H)
        AOD(1,1) = angles(3,1);
        AOD(2,1) = angles(1,1);
        % AOA (1. V / 2. H)
        AOA(1,1) = angles(4,1);
        AOA(2,1) = angles(2,1);
        
        % Vertical polarization changes to (0,180) deg.
        AOD(1,1) = -1*AOD(1,1) + 90;
        AOA(1,1) = -1*AOA(1,1) + 90;
        % Horizontal polarization changes to (0,360) deg.
        AODtmp = AOD(2,1);
        AODtmp_idx = AODtmp < 0;
        AODtmp(AODtmp_idx) = AODtmp(AODtmp_idx) + 360;
        AOD(2,1) = AODtmp;
        AOAtmp = AOA(2,1);
        AOAtmp_idx = AOAtmp < 0;
        AOAtmp(AOAtmp_idx) = AOAtmp(AOAtmp_idx) + 360;
        AOA(2,1) = AOAtmp;
end

