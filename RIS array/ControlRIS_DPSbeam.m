function [ Coeff ] = ControlRIS_DPSbeam( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
%% DPS Setup
DPS_phase = AntData.DPS.phase;
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360) + 360;
DPS_phase(DPS_phase<-360) = DPS_phase(DPS_phase<-360) + 360;
DPS_phase = DPS_phase/180*pi;
DPS_loss = AntData.DPS.loss;
DPS_need = sqrt(10.^(-0.1*DPS_loss)) .* exp(1j*DPS_phase);

%% Control DPS (Perfect beamforming with all phase information)
ris_size = AntData.ris.size;
Phase_LOS = angle(Coeff_LOS);
Phase_Tx2RIS = angle(Coeff_Tx2RIS);
Phase_RIS2Rx = angle(Coeff_RIS2Rx);
Phase_Tx2RIS2Rx = reshape(Phase_Tx2RIS, ris_size) + reshape(Phase_RIS2Rx, ris_size);
Phase_delta = repmat(Phase_LOS, ris_size) - Phase_Tx2RIS2Rx;

%% Corresponding to the closest DPS's phase
Phase_delta_mat = repmat(Phase_delta, [1,1,length(DPS_phase)]);
DPS_phase_mat = repmat(reshape(DPS_phase,[1,1,length(DPS_phase)]), [ris_size(1),ris_size(2),1]);
delta_mat = abs(Phase_delta_mat - DPS_phase_mat);
[~, DPSstate] = min(delta_mat, [], 3);
DPS = DPS_need(DPSstate);

%% Calculate coefficients
DPS_diag = diag(DPS(:));
Coeff_Tx2RIS2Rx = zeros(size(Coeff_RIS2Rx,1), size(Coeff_Tx2RIS,2));
for i = 1:size(Coeff_RIS2Rx,1)
        for j = 1:size(Coeff_Tx2RIS,2)
                Coeff_Tx2RIS2Rx(i, j) = Coeff_RIS2Rx(i,:) * DPS_diag * Coeff_Tx2RIS(:,j);
        end
end
%%% Combine all coefficients
Coeff = Coeff_Tx2RIS2Rx + ...
        repmat(Coeff_scatter, size(Coeff_Tx2RIS2Rx)) + ...
        repmat(Coeff_reflect_RIS, size(Coeff_Tx2RIS2Rx));
end

