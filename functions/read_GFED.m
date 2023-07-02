% This function reads timeseries of global mean GFED biomass burning emissions
% ****************************************************************************

function [sources_CH4,s] = read_GFED(sources_CH4,s,syear,eyear,file_GFED,climatology)

if(isstruct(sources_CH4)==0); clear sources_CH4; end

% Read emissions from file that contains GFEDv41s emis from 1997 onwards and IPCC BB emis. before that
data_read_all = importdata(file_GFED,',',2);
years_file = data_read_all.data(:,1);
data_read = data_read_all.data(:,end);

syear_file = years_file(1);
eyear_file = years_file(end);
nyears_file= length(years_file);

% Extract specified years only and convert 1e10 g yr-1 -> Tg yr-1
nyears = eyear-syear+1;
bb_emis_time = zeros(nyears,1);
sindex2 = 1;
sindex = syear-syear_file+1;
eindex = nyears_file - (eyear_file-eyear);
if(climatology>0)
    disp(['WARNING: BB emissions will be set to climatological (1997-2021) value of 15.73 Tg/yr!']);
    if(climatology==1)
        bb_emis_time(sindex2:nyears) = 1573*0.01; % GFED mean between 1997-2021
    else
        disp('***********************************************************************************');
        disp(['*** WARNING: BB clim. emissions will be multiplied by a factor ',num2str(climatology),'! ***']);
        disp('***********************************************************************************');
        bb_emis_time(sindex2:nyears) = 1573*0.01*climatology;
    end

else
    if(syear<syear_file || eyear>eyear_file)
        disp('ERROR in read_GFED.m: Chosen years are outside the range of IPCC-BB/GFED emissions!'); return;
    end
    bb_emis_time(sindex2:nyears) = data_read(sindex:eindex);% these data are already in Tg
end

% Add information to struct (dC13 value and std. from Zhang et al. (2022))
s=s+1; sources_CH4(s) = struct('name','Biomass burning','dC13',-24.6,'dC13_std',2.2,'emissions',bb_emis_time(:));

% Print value for the most recent year
disp(['  - GFED (',int2str(eyear),'): Biomass burning: ',num2str(bb_emis_time(end)),' Tg yr-1']);
