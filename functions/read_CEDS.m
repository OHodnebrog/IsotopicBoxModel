% This function reads timeseries of global mean CEDS anthropogenic emissions
% **************************************************************************

function [sources_CH4,s] = read_CEDS(sources_CH4,s,syear,eyear,file_CEDS)

if(isstruct(sources_CH4)==0); clear sources_CH4; end

% Read emissions from file
data_read = importdata(file_CEDS,',');
nsources = length(data_read.data(:,1));
nyears_file = length(data_read.data(1,:));
years_file_text = data_read.textdata(1,4:end);
syear_file_tmp = char(years_file_text(1));
syear_file = str2double(syear_file_tmp(3:6));
eyear_file_tmp = char(years_file_text(end));
eyear_file = str2double(eyear_file_tmp(3:6));
if(syear<syear_file || eyear>eyear_file)
    disp('ERROR in read_CEDS.m: Chosen years are outside the range of CEDS emissions!'); return;
end

% Define combined sources with separate isotopic ratio
combined_sources = char('Coal/gas/oil/industry','Livestock','Rice','Waste');
combined_dC13    = [     0.5*(-43.7+-44.7),      -65.5,      -62.2, -56.0]; % Zhang2022
combined_dC13_std= [     0.5*(3.8+4.5),          3.3,        3.8,   7.1]; % Zhang2022
ncombined_sources = length(combined_dC13);

% Loop through sources
combined_emis = zeros(ncombined_sources,nyears_file);
for n = 1:nsources
    if(sum(data_read.data(n,:))>0)
        source = char(data_read.textdata(n+1,2));
        
        % Coal/gas/oil/industry
        if(strcmp(source(1),'1')||strcmp(source(1),'2')||strcmp(source(1),'7'))
            combined_emis(1,:) = combined_emis(1,:) + data_read.data(n,:);
            
        % Livestock
        elseif(strcmp(source(1:2),'3B')||strcmp(source(1:2),'3E'))
            combined_emis(2,:) = combined_emis(2,:) + data_read.data(n,:);
            
        % Rice
        elseif(strcmp(source(1:2),'3D'))
            combined_emis(3,:) = combined_emis(3,:) + data_read.data(n,:);
            
        % Waste
        elseif(strcmp(source(1),'5'))
            combined_emis(4,:) = combined_emis(4,:) + data_read.data(n,:);
        
        else
            disp('ERROR in read_CEDS.m: Unknown CEDS source!'); return;
        end
    end
end

% Extract specified years only and convert kt yr-1 -> Tg yr-1
sindex = syear-syear_file+1;
eindex = nyears_file - (eyear_file-eyear);
combined_emis_time = combined_emis(:,sindex:eindex)*1e-3;

% Add information to struct
for n = 1:ncombined_sources
    s=s+1; sources_CH4(s) = struct('name',strtrim(combined_sources(n,:)),'dC13',combined_dC13(n),'dC13_std',combined_dC13_std(n),'emissions',combined_emis_time(n,:));
    
    % Print values for the most recent year
    disp(['  - CEDS (',int2str(eyear),'): ',combined_sources(n,:),': ',num2str(combined_emis_time(n,end)),' Tg yr-1']);
end
