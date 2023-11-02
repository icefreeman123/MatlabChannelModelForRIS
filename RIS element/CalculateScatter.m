function [ Coeff, Results ] = CalculateScatter( AntData )
%% Setup
        tx = AntData.tx;
        ris = AntData.ris;
        rx = AntData.rx;
        tx_position = AntData.tx.position;
        ris_position = AntData.ris.position;
        rx_position = AntData.rx.position;

        freq = AntData.freq;
        LossType = AntData.LossType;
        speed = physconst('lightspeed'); % speed of light(unit: m/s)
        lambda = speed/freq; % wavelength(unit: Hz)
        k = (2*pi)/lambda; % Wave number

        a = AntData.ris.area_a;
        b = AntData.ris.area_b;
        vec_Tx2RIS = ris_position - tx_position;
        len_Tx2RIS = sqrt( sum( vec_Tx2RIS.^2 ,1 ) );
        vec_RIS2Rx = rx_position - ris_position;
        len_RIS2Rx = sqrt( sum( vec_RIS2Rx.^2 ,1 ) );
        len_total = len_Tx2RIS + len_RIS2Rx;
        ScatterEff = 1; % 11.2 for 4x4 array
        
%% Amplitude
        % Angle of Tx2RIS
        [AOD_Tx2RIS, AOA_Tx2RIS] = CalculateAngle( tx_position, ris_position ); % Angles between Tx and Rx
        aod_r_Tx2RIS = AOD_Tx2RIS(2).'*pi/180;
        aoa_r_Tx2RIS = AOA_Tx2RIS(2).'*pi/180; % phi_i
        eod_r_Tx2RIS = AOD_Tx2RIS(1).'*pi/180;
        eoa_r_Tx2RIS = AOA_Tx2RIS(1).'*pi/180; % theta_i
        % Sigmal of Tx2RIS
        [ Vt, Ht ] = LinearInterpolate( tx , aod_r_Tx2RIS , eod_r_Tx2RIS );
        V_Tx2RIS = 1 *Vt;
        H_Tx2RIS = 1 *Ht;

        % Angle of RIS2Rx
        [AOD_RIS2Rx, AOA_RIS2Rx] = CalculateAngle( ris_position, rx_position ); % Angles between RIS and Rx
        aod_r_RIS2Rx = AOD_RIS2Rx(2).'*pi/180; % phi_s
        aoa_r_RIS2Rx = AOA_RIS2Rx(2).'*pi/180;
        eod_r_RIS2Rx = AOD_RIS2Rx(1).'*pi/180; % theta_s
        eoa_r_RIS2Rx = AOA_RIS2Rx(1).'*pi/180;
        % Sigmal of RIS2Rx
        [ Vr, Hr ] = LinearInterpolate( rx , aoa_r_RIS2Rx , eoa_r_RIS2Rx );
        V_RIS2Rx = Vr *1;
        H_RIS2Rx = Hr *1;

        % Scatter Matrix (From C. A. Balanis, Advanced Engineering Electromagnetics. 
        % 2nd ed. Hoboken, NJ, USA: Wiley, 2012. pp.591-599.)
        X = (a*k/2) * sin(eod_r_RIS2Rx) * cos(aod_r_RIS2Rx);
        Y = (b*k/2) * (sin(eod_r_RIS2Rx) * sin(aod_r_RIS2Rx) - sin(eoa_r_Tx2RIS));
        if isnan(sin(X)/X), Xsinc = 1; else Xsinc = sin(X)/X; end
        if isnan(sin(Y)/Y), Ysinc = 1; else Ysinc = sin(Y)/Y; end
        C = (2*a*b*k) * Xsinc * Ysinc;
        S_V2V = C * cos(eoa_r_Tx2RIS) * cos(eod_r_RIS2Rx) * cos(aod_r_RIS2Rx);
        S_V2H = C * cos(eoa_r_Tx2RIS) * sin(aod_r_RIS2Rx);
        S_H2V = -1j * C * cos(eod_r_RIS2Rx) * sin(aod_r_RIS2Rx);
        S_H2H = -1j * C * cos(aod_r_RIS2Rx);
        S_matrix = [S_V2V, S_H2V; S_V2H, S_H2H];
        
        % Scatter coefficient
        Amp = V_RIS2Rx * (S_V2V*V_Tx2RIS - S_H2V*H_Tx2RIS) - H_RIS2Rx * (S_V2H*V_Tx2RIS - S_H2H*H_Tx2RIS);

%% Phase        
        psi_lms = k * mod(len_total, lambda); 

%% Path Loss
        switch LossType
                case 'withloss'
                        loss_dB_Tx2RIS = pathloss( tx_position, ris_position, lambda );
                        loss_dB_RIS2Rx = pathloss( ris_position, rx_position, lambda );    
                case 'noloss' % otherwise
                        loss_dB_Tx2RIS = 0;
                        loss_dB_RIS2Rx = 0;
        end
        loss = sqrt( 10.^( 0.1 * -1 * loss_dB_Tx2RIS ) ) * sqrt( 10.^( 0.1 * -1 * loss_dB_RIS2Rx ) );

%% Calculate Coeff.
        Coeff = loss * Amp * exp(-1j*psi_lms);
        Coeff = ScatterEff * Coeff;
        
%% Results
        Results.S_matrix = S_matrix;
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
