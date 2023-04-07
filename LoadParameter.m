function [ AntData ] = LoadParameter( file, TxLoc, RxLoc, RISLoc )
%% Setup
TxPower_dBm = 0; % unit: dBm.
AntData.TxPower_dBm = TxPower_dBm;
AntData.TxPower = 1e-3 * 10^(TxPower_dBm/10); % unit: W.
AntData.LossType = 'withloss'; % 'withloss' / 'noloss'
AntData.freq = 3.5e9; % Center Frequency (unit: Hz)

elevation_range = [0, pi];
azimuth_range = [0, 2*pi];
AntData.elevation_range = elevation_range;
AntData.azimuth_range = azimuth_range;

if exist('RISLoc','var')
        AntNum = 3;
else
        AntNum = 2;
end
%% Load radiation pattern
startIDX = 2;
path_now = pwd;
path_data = [path_now,'\', file];
cd([path_data]);
listDateTMP = dir;
listDate = cell((length(listDateTMP)-2), 2);
for i=1:size(listDate,1)
        listDate(i,1) = {listDateTMP(i+2).name}; 
end
for i=1:size(listDate,1)
        filename = listDate{i,1};
        ANTcsvTMP = xlsread([filename]);
        ANTcsv = ANTcsvTMP(startIDX:end , startIDX:end);
        listDate(i,2) = {ANTcsv}; clear ANTcsv
end
cd([path_now]);

%%
typeAll = {'tx','rx','ris'};
AntNeedAll = [1, 3, 2];

for i = 1:AntNum
        type = typeAll{1,i};
        AntNeed = AntNeedAll(1,i);
        
        CustAntV = listDate(AntNeed + (size(listDate,1)/2),:);
        CustAntH = listDate(AntNeed,:);
        Vtmp = CustAntV{1,2};
        Htmp = CustAntH{1,2};
        Vtmp = sqrt(10.^(Vtmp/10));
        Htmp = sqrt(10.^(Htmp/10));
        no_elNeed = size(Vtmp, 1);
        no_azNeed = size(Vtmp, 2);        
        elevation_gridNeed = linspace(elevation_range(1), elevation_range(2), no_elNeed);
        azimuth_gridNeed = linspace(azimuth_range(1), azimuth_range(2), no_azNeed);

        switch type
                case 'tx'
                        AntData.tx.name = 'patch';
                        AntData.tx.V = Vtmp;
                        AntData.tx.H = Htmp;
                        AntData.tx.T = sqrt(Vtmp.^2 + Htmp.^2);                        
                        AntData.tx.no_el = no_elNeed;
                        AntData.tx.no_az = no_azNeed;
                        AntData.tx.elevation_grid = elevation_gridNeed;
                        AntData.tx.azimuth_grid = azimuth_gridNeed;
                        AntData.tx.position = TxLoc; 
                case 'rx'
                        AntData.rx.name = 'patch';
                        AntData.rx.V = Vtmp;
                        AntData.rx.H = Htmp;
                        AntData.rx.T = sqrt(Vtmp.^2 + Htmp.^2);                        
                        AntData.rx.no_el = no_elNeed;
                        AntData.rx.no_az = no_azNeed;
                        AntData.rx.elevation_grid = elevation_gridNeed;
                        AntData.rx.azimuth_grid = azimuth_gridNeed;
                        AntData.rx.position = RxLoc;
                case 'ris'
                        AntData.ris.name = 'patch';
                        AntData.ris.V = Vtmp;
                        AntData.ris.H = Htmp;
                        AntData.ris.T = sqrt(Vtmp.^2 + Htmp.^2);
                        AntData.ris.no_el = no_elNeed;
                        AntData.ris.no_az = no_azNeed;
                        AntData.ris.elevation_grid = elevation_gridNeed;
                        AntData.ris.azimuth_grid = azimuth_gridNeed;
                        AntData.ris.position = RISLoc;
                        AntData.ris.area_a = 0.094; % unit: m.
                        AntData.ris.area_b = 0.073;
                        AntData.ris.normal = [0;0;1]; % Normal vector of RIS
                        AntData.ris.InjEff = 1;
        end
end

end

