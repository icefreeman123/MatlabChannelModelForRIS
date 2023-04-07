function [ Coeff, Results ] = CalculateReflect( AntData )
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
        RefEff = 8;
        
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

        % Vector
        NorVec = AntData.ris.normal; % Normal vector of RIS
        Vec_RIS2Tx = (tx_position-ris_position)/norm(tx_position-ris_position);
        Vec_RIS2Rx = (rx_position-ris_position)/norm(rx_position-ris_position);
        Vec_Combine = (Vec_RIS2Tx+Vec_RIS2Rx)/norm(Vec_RIS2Tx+Vec_RIS2Rx);

        %%% Reflection Matrix (Good conductor)
        R_V = -1;
        R_H = 1;
        
        % Reflection coefficient
        Amp = V_RIS2Rx * (R_V*V_Tx2RIS - 0) - H_RIS2Rx * (0 - R_H*H_Tx2RIS);
%         if isequal(sum(abs(NorVec-Vec_Combine)), 0) % Incident Angle = Reflection Angle
%                 Amp = V_RIS2Rx * (R_V*V_Tx2RIS - 0) - H_RIS2Rx * (0 - R_H*H_Tx2RIS);
%         else % Incident Angle ~= Reflection Angle
%                 Amp = 0; 
%                 % It is not correct in the practical environment because 
%                 % the reflected signal exists except the reflection angle.
%         end

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
        Coeff = RefEff * Coeff;

%% Results
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

