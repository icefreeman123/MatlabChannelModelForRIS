function [ Coeff_mat ] = CalculateScatter( AntData, RISLoc )
%% Setup
        freq = AntData.freq;
        LossType = AntData.LossType;
        speed = physconst('lightspeed'); % speed of light(unit: m/s)
        lambda = speed/freq; % wavelength(unit: Hz)
        k = (2*pi)/lambda; % Wave number
        
        tx = AntData.tx;
        ris = AntData.ris;
        rx = AntData.rx;
        tx_position = tx.position; % Each antenna location.
        ris_position = RISLoc; % Array center.
        rx_position = rx.position; % Each antenna location.
        tx_V = tx.V;
        tx_H = tx.H;
        rx_V = rx.V;
        rx_H = rx.H;
        no_ant_tx = tx.no_ant;
        no_ant_rx = rx.no_ant;
        
        size_ris = ris.size;
        a = AntData.ris.area_a;
        b = AntData.ris.area_b;
        a_all = a * size_ris(1,1);
        b_all = b * size_ris(1,2);
        
%% Amplitude
        tx_position_tmp = tx_position(:);
        rx_position_tmp = rx_position(:);
        tx_Vtmp = tx_V(:);
        tx_Htmp = tx_H(:);
        rx_Vtmp = rx_V(:);
        rx_Htmp = rx_H(:);
        V_Tx2RIS_mat = zeros(1, no_ant_tx);
        H_Tx2RIS_mat = zeros(1, no_ant_tx);
        V_RIS2Rx_mat = zeros(1, no_ant_rx);
        H_RIS2Rx_mat = zeros(1, no_ant_rx);
        aoa_r_Tx2RIS_mat = zeros(1, no_ant_tx);
        eoa_r_Tx2RIS_mat = zeros(1, no_ant_tx);
        aod_r_RIS2Rx_mat = zeros(1, no_ant_rx);
        eod_r_RIS2Rx_mat = zeros(1, no_ant_rx);
        S_matrix_cell = cell(no_ant_rx, no_ant_tx);
        Amp_mat = zeros(no_ant_rx, no_ant_tx);
        % Tx2RIS
        for i = 1:no_ant_tx
                tx_position_tmp2 = tx_position_tmp{i};
                tx_Vtmp2 = tx_Vtmp{i};
                tx_Htmp2 = tx_Htmp{i};                
                % Angle of Tx2RIS
                [AOD_Tx2RIS, AOA_Tx2RIS] = CalculateAngle( tx_position_tmp2, ris_position ); % Angles between Tx and Rx
                aod_r_Tx2RIS = AOD_Tx2RIS(2).'*pi/180;
                aoa_r_Tx2RIS = AOA_Tx2RIS(2).'*pi/180; % phi_i
                eod_r_Tx2RIS = AOD_Tx2RIS(1).'*pi/180;
                eoa_r_Tx2RIS = AOA_Tx2RIS(1).'*pi/180; % theta_i
                aoa_r_Tx2RIS_mat(1, i) = aoa_r_Tx2RIS;
                eoa_r_Tx2RIS_mat(1, i) = eoa_r_Tx2RIS;
                % Signal of Tx2RIS
                [ Vt, Ht ] = LinearInterpolate( tx_Vtmp2, tx_Htmp2, tx , aod_r_Tx2RIS , eod_r_Tx2RIS );                        
                V_Tx2RIS = 1 *Vt;
                H_Tx2RIS = 1 *Ht;
                V_Tx2RIS_mat(1, i) = V_Tx2RIS;
                H_Tx2RIS_mat(1, i) = H_Tx2RIS;
        end
        % RIS2Rx
        for j = 1:no_ant_rx
                rx_position_tmp2 = rx_position_tmp{j};
                rx_Vtmp2 = rx_Vtmp{j};
                rx_Htmp2 = rx_Htmp{j};                        
                % Angle of RIS2Rx
                [AOD_RIS2Rx, AOA_RIS2Rx] = CalculateAngle( ris_position, rx_position_tmp2 ); % Angles between RIS and Rx
                aod_r_RIS2Rx = AOD_RIS2Rx(2).'*pi/180; % phi_s
                aoa_r_RIS2Rx = AOA_RIS2Rx(2).'*pi/180;
                eod_r_RIS2Rx = AOD_RIS2Rx(1).'*pi/180; % theta_s
                eoa_r_RIS2Rx = AOA_RIS2Rx(1).'*pi/180;
                aod_r_RIS2Rx_mat(1, j) = aod_r_RIS2Rx;
                eod_r_RIS2Rx_mat(1, j) = eod_r_RIS2Rx;
                % Signal of RIS2Rx
                [ Vr, Hr ] = LinearInterpolate( rx_Vtmp2, rx_Htmp2, rx , aoa_r_RIS2Rx , eoa_r_RIS2Rx );
                V_RIS2Rx = Vr *1;
                H_RIS2Rx = Hr *1;
                V_RIS2Rx_mat(1, j) = V_RIS2Rx;
                H_RIS2Rx_mat(1, j) = H_RIS2Rx;
        end
        
        % RIS scatter matrix corresponding to the angle of Tx and Rx
        % (Scatter Matrix is from C. A. Balanis, Advanced Engineering Electromagnetics. 
        % 2nd ed. Hoboken, NJ, USA: Wiley, 2012, pp.591-599.)
        for i = 1:no_ant_tx
                V_Tx2RIS = V_Tx2RIS_mat(1, i);
                H_Tx2RIS = H_Tx2RIS_mat(1, i);
                aoa_r_Tx2RIS = aoa_r_Tx2RIS_mat(1, i);
                eoa_r_Tx2RIS = eoa_r_Tx2RIS_mat(1, i);
                for j = 1:no_ant_rx
                        V_RIS2Rx = V_RIS2Rx_mat(1, j);
                        H_RIS2Rx = H_RIS2Rx_mat(1, j);
                        aod_r_RIS2Rx = aod_r_RIS2Rx_mat(1, j);
                        eod_r_RIS2Rx = eod_r_RIS2Rx_mat(1, j);
                        % Calculate scatter matrix
                        X = (a_all*k/2) * sin(eod_r_RIS2Rx) * cos(aod_r_RIS2Rx);
                        Y = (b_all*k/2) * (sin(eod_r_RIS2Rx) * sin(aod_r_RIS2Rx) - sin(eoa_r_Tx2RIS));
                        if isnan(sin(X)/X), Xsinc = 1; else Xsinc = sin(X)/X; end
                        if isnan(sin(Y)/Y), Ysinc = 1; else Ysinc = sin(Y)/Y; end
                        C = (2*a_all*b_all*k) * Xsinc * Ysinc;
                        S_V2V = C * cos(eoa_r_Tx2RIS) * cos(eod_r_RIS2Rx) * cos(aod_r_RIS2Rx);
                        S_V2H = C * cos(eoa_r_Tx2RIS) * sin(aod_r_RIS2Rx);
                        S_H2V = -1j * C * cos(eod_r_RIS2Rx) * sin(aod_r_RIS2Rx);
                        S_H2H = -1j * C * cos(aod_r_RIS2Rx);
                        S_matrix = [S_V2V, S_H2V; S_V2H, S_H2H];
                        S_matrix_cell(j, i) = {S_matrix};
                        % Scatter coefficient
                        Amp = V_RIS2Rx * (S_V2V*V_Tx2RIS - S_H2V*H_Tx2RIS) - H_RIS2Rx * (S_V2H*V_Tx2RIS - S_H2H*H_Tx2RIS);
                        Amp_mat(j, i) = Amp;
                end
        end

%% Phase (Path length is "i-th Tx -> RIS center -> j-th Rx")
        tx_position_tmp = tx_position(:);
        rx_position_tmp = rx_position(:);
        Phase_mat = zeros(no_ant_rx, no_ant_tx);
        for i = 1:no_ant_tx
                tx_position_tmp2 = tx_position_tmp{i}; % i-th Tx location
                vec_Tx2RIS = ris_position - tx_position_tmp2;
                len_Tx2RIS = sqrt( sum( vec_Tx2RIS.^2 ,1 ) );
                for j = 1:no_ant_rx
                        rx_position_tmp2 = rx_position_tmp{j}; % j-th Rx location
                        vec_RIS2Rx = rx_position_tmp2 - ris_position;
                        len_RIS2Rx = sqrt( sum( vec_RIS2Rx.^2 ,1 ) );
                        len_total = len_Tx2RIS + len_RIS2Rx;
                        Phase = k * mod(len_total, lambda); 
                        Phase_mat(j, i) = Phase;
                end
        end

%% Path Loss
        switch LossType
                case 'withloss'
                        loss_dB_Tx2RIS_mat = zeros(no_ant_rx, no_ant_tx);
                        loss_dB_RIS2Rx_mat = zeros(no_ant_rx, no_ant_tx);
                        for i = 1:no_ant_tx
                                tx_position_tmp2 = tx_position_tmp{i};
                                for j = 1:no_ant_rx
                                        rx_position_tmp2 = rx_position_tmp{j};
                                        loss_dB_Tx2RIS = pathloss( tx_position_tmp2, ris_position, lambda );
                                        loss_dB_RIS2Rx = pathloss( ris_position, rx_position_tmp2, lambda );
                                        loss_dB_Tx2RIS_mat(j, i) = loss_dB_Tx2RIS;
                                        loss_dB_RIS2Rx_mat(j, i) = loss_dB_RIS2Rx;
                                end
                        end                        
                case 'noloss' % otherwise
                        loss_dB_Tx2RIS_mat = zeros(no_ant_rx, no_ant_tx);
                        loss_dB_RIS2Rx_mat = zeros(no_ant_rx, no_ant_tx);
        end
        loss_mat = (sqrt(10.^(0.1*-1*loss_dB_Tx2RIS_mat))) .* (sqrt(10.^(0.1*-1*loss_dB_RIS2Rx_mat)));

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
