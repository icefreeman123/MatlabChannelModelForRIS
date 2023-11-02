function [ Coeff ] = ControlDPS_BGApolar( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
%% BGA
%%% Input
ris_initial_diag = AntData.ris.initial_diag;
ris_num = AntData.ris.no_ant; % Number of RIS antennas
ris_size = AntData.ris.size;
Tr = 100;
Tg = 2;
state_need = 1:16; % DPS state needed to be controlled
Bbit = state_need;
state_num = length(Bbit);
bit_num = log10(state_num) / log10(2);
codebook_num = (2^bit_num)^ris_num;
polar_need = [0,45,90];
polar_num = length(polar_need); % Number of RIS polarization
%%% Initialize
Br = ones(1, ris_num);
Br_pol = zeros(ris_num, 3);
Btmp = ones(1, ris_num);
Btmp_pol = zeros(ris_num, 3);
Coeff_all_stack = [];
Pol_axis = 3;
%% DPS Setup
DPS_phase = AntData.DPS.phase;
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360)+360;
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360)+360;
DPS_loss = AntData.DPS.loss;
DPS_phase_need = DPS_phase(1,state_need);
DPS_loss_need = DPS_loss(1,state_need);
DPS_data_need = [DPS_phase_need ; DPS_loss_need].';
DPS_data_need = sortrows(DPS_data_need, -1).';
DPS_need = sqrt(10.^(-0.1*DPS_data_need(2,:))) .* exp(1j*(DPS_data_need(1,:)/180*pi));
%% RMS
% (1)
DPS = DPS_need(Btmp);
DPS_diag = diag(DPS(:));
Coeff_Tx2RIS2Rx = zeros(size(Coeff_RIS2Rx,1), size(Coeff_Tx2RIS,2));
for i = 1:size(Coeff_RIS2Rx,1)
        for j = 1:size(Coeff_Tx2RIS,2)
                Coeff_Tx2RIS2Rx(i, j) = Coeff_RIS2Rx(i,:) * DPS_diag * ris_initial_diag * Coeff_Tx2RIS(:,j);
        end
end
Coeff = Coeff_Tx2RIS2Rx + ...
        repmat(Coeff_scatter, size(Coeff_Tx2RIS2Rx)) + ...
        repmat(Coeff_reflect_RIS, size(Coeff_Tx2RIS2Rx)); %%% Combine all coefficients
Coeff_need = Coeff;
P = abs(Coeff).^2;
for tr = 1:Tr
        % (2)
        Btmp = [];
        for i = 1:ris_num
                rand( 'state', tr*i ); 
                Btmp_value = randi(state_num);
                Btmp = [Btmp, Btmp_value];
        end
        % (3)
        DPS = DPS_need( flip(Btmp) );
        DPS_diag = diag(DPS(:));
        Coeff_Tx2RIS2Rx = zeros(size(Coeff_RIS2Rx,1), size(Coeff_Tx2RIS,2));
        for i = 1:size(Coeff_RIS2Rx,1)
                for j = 1:size(Coeff_Tx2RIS,2)
                        Coeff_Tx2RIS2Rx(i, j) = Coeff_RIS2Rx(i,:) * DPS_diag * ris_initial_diag * Coeff_Tx2RIS(:,j);
                end
        end
        Coeff = Coeff_Tx2RIS2Rx + ...
                repmat(Coeff_scatter, size(Coeff_Tx2RIS2Rx)) + ...
                repmat(Coeff_reflect_RIS, size(Coeff_Tx2RIS2Rx)); %%% Combine all coefficients
        Ptmp = abs(Coeff).^2;
        % Stack data
        if (Ptmp > P)
                Coeff_need = Coeff;
        end
        % (4)
        if (Ptmp > P)
                P = Ptmp;
                Br = Btmp;
                Br_pol = Btmp_pol;
        end
end
% (5)
Bg = Br;
Bg_pol = Br_pol;
for i = 1:ris_num             
        Bg_pol_tmp = Bg_pol(i,:);
        AntData.ris = RotateAntenna( AntData.ris, Bg_pol_tmp, i );
end
%% GS
for tg = 1:Tg
        for nris = 1:ris_num
                Btmp = Bg;
                Btmp_pol = Bg_pol;
                for mbit = 1:state_num
                        for mpol = 1:polar_num
                                % (6)
                                Btmp(nris) = Bbit(mbit);
                                Btmp_pol(nris, Pol_axis) = polar_need(mpol);
                                % (7)
                                %%% Polarization
                                Btmp_pol_tmp = Btmp_pol(nris, :);
                                AntData.ris = RotateAntenna( AntData.ris, Btmp_pol_tmp, nris );
                                Coeff_Tx2RIS = CalculateCC( AntData, AntData.tx, AntData.ris );
                                Coeff_RIS2Rx = CalculateCC( AntData, AntData.ris, AntData.rx );        
                                %%% State
                                DPS = DPS_need( flip(Btmp) );
                                DPS_diag = diag(DPS(:));
                                Coeff_Tx2RIS2Rx = zeros(size(Coeff_RIS2Rx,1), size(Coeff_Tx2RIS,2));
                                for i = 1:size(Coeff_RIS2Rx,1)
                                        for j = 1:size(Coeff_Tx2RIS,2)
                                                Coeff_Tx2RIS2Rx(i, j) = Coeff_RIS2Rx(i,:) * DPS_diag * ris_initial_diag * Coeff_Tx2RIS(:,j);
                                        end
                                end
                                Coeff = Coeff_Tx2RIS2Rx + ...
                                        repmat(Coeff_scatter, size(Coeff_Tx2RIS2Rx)) + ...
                                        repmat(Coeff_reflect_RIS, size(Coeff_Tx2RIS2Rx)); %%% Combine all coefficients
                                Ptmp = abs(Coeff).^2;
                                % Stack data
                                if (Ptmp > P)
                                        Coeff_need = Coeff;
                                end
                                % (8)
                                if (Ptmp > P)
                                        P = Ptmp;
                                        Bg = Btmp;
                                        Bg_pol = Btmp_pol;
                                end
                                % Recover polarization
                                AntData.ris = RotateAntenna( AntData.ris, -1*Btmp_pol_tmp, nris );    
                        end
                end
        end
end
% (9)
B = Bg;
B_pol = Bg_pol(:,Pol_axis); 
B_reshape = reshape(flip(B), ris_size);
B_pol_reshape = reshape(flip(B_pol), ris_size);
%% Results
Coeff = Coeff_need;
end

