function [ Coeff_mat ] = CalculateCC( AntData, tx, rx )
%% Setup
        freq = AntData.freq;
        LossType = AntData.LossType;
        speed = physconst('lightspeed'); % speed of light(unit: m/s)
        lambda = speed/freq; % wavelength(unit: Hz)
        k = (2*pi)/lambda; % Wave number     

        tx_position = tx.position; % Each antenna location.
        rx_position = rx.position; % Each antenna location.
        tx_V = tx.V;
        tx_H = tx.H;
        rx_V = rx.V;
        rx_H = rx.H;
        no_ant_tx = tx.no_ant;
        no_ant_rx = rx.no_ant;
        
%% Amplitude
        tx_position_tmp = tx_position(:);
        rx_position_tmp = rx_position(:);
        tx_Vtmp = tx_V(:);
        tx_Htmp = tx_H(:);
        rx_Vtmp = rx_V(:);
        rx_Htmp = rx_H(:);        
        V_final_mat = zeros(no_ant_rx, no_ant_tx);
        H_final_mat = zeros(no_ant_rx, no_ant_tx);
        Amp_mat = zeros(no_ant_rx, no_ant_tx);
        for i = 1:no_ant_tx
                tx_position_tmp2 = tx_position_tmp{i};
                tx_Vtmp2 = tx_Vtmp{i};
                tx_Htmp2 = tx_Htmp{i};
                for j = 1:no_ant_rx
                        rx_position_tmp2 = rx_position_tmp{j};
                        rx_Vtmp2 = rx_Vtmp{j};
                        rx_Htmp2 = rx_Htmp{j};
                        %                        
                        [AOD, AOA] = CalculateAngle( tx_position_tmp2, rx_position_tmp2 ); % Angles between Tx and Rx
                        aod_r = AOD(2).'*pi/180;
                        aoa_r = AOA(2).'*pi/180;
                        eod_r = AOD(1).'*pi/180;
                        eoa_r = AOA(1).'*pi/180;
                        %                        
                        [ Vt, Ht ] = LinearInterpolate( tx_Vtmp2, tx_Htmp2, tx , aod_r , eod_r );
                        [ Vr, Hr ] = LinearInterpolate( rx_Vtmp2, rx_Htmp2, rx , aoa_r , eoa_r );
                        V_final = Vr *Vt;
                        H_final = Hr *Ht;
                        Amp = V_final - H_final;
                        %                        
                        V_final_mat(j, i) = V_final;
                        H_final_mat(j, i) = H_final;
                        Amp_mat(j, i) = Amp;
                end
        end

%% Phase
        Phase_mat = zeros(no_ant_rx, no_ant_tx);
        for i = 1:no_ant_tx
                tx_position_tmp2 = tx_position_tmp{i};
                for j = 1:no_ant_rx
                        rx_position_tmp2 = rx_position_tmp{j};
                        r_s = rx_position_tmp2 - tx_position_tmp2;
                        norm_r_s = sqrt( sum( r_s.^2 ,1 ) );
                        Phase = k * mod(norm_r_s, lambda); 
                        Phase_mat(j, i) = Phase;
                end
        end       
        
%% Path Loss
        switch LossType
                case 'withloss'
                        loss_dB_mat = zeros(no_ant_rx, no_ant_tx);
                        for i = 1:no_ant_tx
                                tx_position_tmp2 = tx_position_tmp{i};
                                for j = 1:no_ant_rx
                                        rx_position_tmp2 = rx_position_tmp{j};
                                        loss_dB = pathloss( tx_position_tmp2, rx_position_tmp2, lambda );
                                        loss_dB_mat(j, i) = loss_dB;
                                end
                        end                        
                case 'noloss' % otherwise
                        loss_dB_mat = zeros(no_ant_rx, no_ant_tx);
        end
        loss_mat = sqrt( 10.^( 0.1 * -1 * loss_dB_mat ) );     
%% Calculate Coeff.
        Coeff_mat = loss_mat .* Amp_mat .* exp(-1j*Phase_mat);
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
        angles(1,1) = atan2( Y_vec, X_vec ); % [–pi, pi] 
        angles(1, isnan(angles(1,1)) ) = 0;
        % AOD-Vertical polarization
        angles(3,1) = atan2( Z_vec, sqrt(X_vec.^2 + Y_vec.^2) ); % [–pi, pi] 
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

