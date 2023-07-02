% This function reads global mean natural emissions
% *************************************************

function [sources_CH4,s,eyear_nat] = read_natural(sources_CH4,s,syear,eyear,nat_emis_fromCLM,file_CLM);

if(isstruct(sources_CH4)==0); clear sources_CH4; end

nyears = eyear-syear+1;


% First wetland emissions
name = 'Wetlands';
dC13 = -61.6; dC13_std = 5.3; % Zhang et al. (2022)

if(nat_emis_fromCLM>0) % use wetland emissions from CLM or Zhang et al. (2023)
    if(nat_emis_fromCLM==1)
        clm_emis_label = 'CLM GSWP3';
        eyear_nat = 2014;
    elseif(nat_emis_fromCLM==2)
        clm_emis_label = 'CLM CRUNCEP';
        eyear_nat = 2016;
    elseif(nat_emis_fromCLM==3)
        clm_emis_label = 'CLM CRUJRA';
        eyear_nat = 2019;
    elseif(nat_emis_fromCLM==10)
        clm_emis_label = 'Zhang23 MERRA2';
        eyear_nat = 2020;
    elseif(nat_emis_fromCLM==11)
        clm_emis_label = 'Zhang23 CRU';
        eyear_nat = 2020;
    else
        disp('ERROR in read_natural.m: Wrong option for nat_emis_fromCLM!'); return;
    end
    column = 1+nat_emis_fromCLM;
    
    data_read_all = importdata(file_CLM,',',2);
    years_file = data_read_all.data(:,1);
    data_read = data_read_all.data(:,column);
    syear_file = years_file(1); eyear_file = years_file(end);
    nyears_file= length(years_file);

    % Extract specified years only
    clm_emis_time = zeros(nyears,1);
    sindex = syear-syear_file+1;
    eindex = nyears_file - (eyear_file-eyear);
    if(syear<syear_file || eyear>eyear_file)
        disp('ERROR in read_natural.m: Chosen years are outside the range of CLM emissions!'); return;
    end
    clm_emis_time(1:nyears) = data_read(sindex:eindex); % these data are already in Tg

    % Subtract other natural emissions (66 Tg/yr in the CH4 box model)
    % and replace the constant wetland emissions with time-varying emissions
    emis(1:nyears) = clm_emis_time(1:nyears) - 66;
    disp(['NOTE: Wetland emissions replaced with ',clm_emis_label]);

else % use constant wetland emissions (same as in CH4 box model)
    eyear_nat = eyear;
    emis(1:nyears) = 149; % Tg/yr
end
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);


% For the remaining emission categories, scale with the ratio
% between other natural emissions from the CH4 box model (215-149 Tg/yr)
% and the bottom-up numbers from the global methane budget
natural_ratio = 66/(222-1); % numbers for "other natural sources" in Saunois et al. (2020)

name = 'Freshwaters';
dC13 = -61.5; dC13_std = 5.0; % Zhang et al. (2022)
emis(1:nyears) = 159*natural_ratio; % Saunois et al. (2020)
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);

name = 'Geological';
dC13 = -42.5; dC13_std = 7.0; % Zhang et al. (2022)
emis(1:nyears) = (38+7)*natural_ratio; % Saunois et al. (2020)
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);

name = 'Wild animals';
dC13 = -65.4; dC13_std = 2.8; % Zhang et al. (2022)
emis(1:nyears) = 2*natural_ratio; % Saunois et al. (2020)
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);

name = 'Termites';
dC13 = -63.4; dC13_std = 5.7; % Zhang et al. (2022)
emis(1:nyears) = 9*natural_ratio; % Saunois et al. (2020)
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);

% Permafrost is skipped for now (unknown dC13 value, and it is only 1 Tg yr-1)

name = 'Oceanic biogenic';% open and coastal';
dC13 = -61.5; dC13_std = 5.0; % assumed same as freshwaters
emis(1:nyears) = 6*natural_ratio; % Saunois et al. (2020)
s=s+1; sources_CH4(s) = struct('name',name,'dC13',dC13,'dC13_std',dC13_std,'emissions',emis(:));
disp(['  - NAT  (',int2str(eyear),'): ',name,': ',num2str(emis(end)),' Tg yr-1']);
