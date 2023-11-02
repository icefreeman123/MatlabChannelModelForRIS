function [ AntDataFold ] = ArrayGenerate( AntData, AntDataFold )
%% Setup
V = AntDataFold.V;
H = AntDataFold.H;
T = AntDataFold.T;
Array_a = AntDataFold.size(1); % Array size (Vertical direction)
Array_b = AntDataFold.size(2); % Array size (Horizontal direction)
no_ant = Array_a * Array_b;

ArrayCenterLoc = AntDataFold.position{1}; % Array center. Unit: m.
ArrayCenterNorVect = AntDataFold.normal; % Normal vector of array
freq = AntData.freq;
speed = physconst('lightspeed'); % speed of light(unit: m/s)
lambda = speed/freq; % wavelength(unit: Hz)
AntDist = lambda/2; % Distance between two antennas in the array.

%% Duplicate pattern
V_cell = cell(Array_a, Array_b);
[V_cell{:}] = deal(cell2mat(V));
H_cell = cell(Array_a, Array_b);
[H_cell{:}] = deal(cell2mat(H));
T_cell = cell(Array_a, Array_b);
[T_cell{:}] = deal(cell2mat(T));

%% Duplicate location 
% Original antenna location (default)
AllPt = zeros(no_ant, 3);
GridPt_a = AntDist * (0:(Array_a-1));
GridPt_b = AntDist * (0:(Array_b-1));
AllPt(:,3) = repmat(GridPt_a.', [Array_b,1]);
AllPt(:,2) = reshape(repmat(GridPt_b, [Array_a,1]), [],1);
AllPt_cell = mat2cell(AllPt, ones(no_ant,1), 3);
AllPt_cell = reshape(AllPt_cell, Array_a, Array_b); % default
AllPtCenter = [0, (AntDist*(Array_b-1)/2), (AntDist*(Array_a-1)/2)]; % default
AllPtAndCenter = [AllPt ; AllPtCenter].';

% Calculate rotation vector (based on Rodrigues' rotation formula)
NorVec_org = [1;0;0]; % default
AllPtAndCenter_new = RodriguesRotVec( NorVec_org, ArrayCenterNorVect, AllPtAndCenter );

% Move center
MoveDist = ArrayCenterLoc - AllPtAndCenter_new(:,end);
AllPtAndCenter_new2 = AllPtAndCenter_new + repmat(MoveDist, 1,size(AllPtAndCenter_new,2));
%%%%% Check: ArrayCenterLoc = AllPtAndCenter_new2(:,end);
AllPt_new = AllPtAndCenter_new2(:, 1:end-1);
AllPt_new_cell = mat2cell(AllPt_new, 3, ones(no_ant,1));
AllPt_new_cell = reshape(AllPt_new_cell, Array_a, Array_b);

%% Output
AntDataFold.V = V_cell;
AntDataFold.H = H_cell;
AntDataFold.T = T_cell;
AntDataFold.position = AllPt_new_cell;
AntDataFold.no_ant = no_ant;

end

