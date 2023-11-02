function [ Coeff ] = ControlRIS_Perfectbeam( AntData, Coeff_scatter, Coeff_reflect_RIS, Coeff_RIS2Rx, Coeff_Tx2RIS, Coeff_LOS );
%% Control DPS (Perfect beamforming with all phase information)
ris_size = AntData.ris.size;
InjEff = 1 *ones(ris_size);
Phase_LOS = angle(Coeff_LOS);
Phase_Tx2RIS = angle(Coeff_Tx2RIS);
Phase_RIS2Rx = angle(Coeff_RIS2Rx);
Phase_Tx2RIS2Rx = reshape(Phase_Tx2RIS, ris_size) + reshape(Phase_RIS2Rx, ris_size);
Phase_delta = repmat(Phase_LOS, ris_size) - Phase_Tx2RIS2Rx;
DPS = InjEff .* exp(1j*Phase_delta);

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

