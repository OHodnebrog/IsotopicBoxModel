% This function reads timeseries of global mean EDGAR anthropogenic emissions
% ***************************************************************************

function [sources_CH4,s] = read_EDGAR(sources_CH4,s,syear,eyear,file_EDGAR)

if(isstruct(sources_CH4)==0); clear sources_CH4; end

% Read emissions from file
syear_file = 1970;
eyear_file = 2021;
data_read = importdata(file_EDGAR,';');
read_sources = squeeze(data_read.textdata(2:end,5));
read_data = data_read.data;
read_data(isnan(read_data)) = 0; % replace NaNs with zeros
nsources = length(read_sources);
nyears_file = eyear_file-syear_file+1;
if(syear<syear_file || eyear>eyear_file)
    disp('ERROR in read_EDGAR.m: Chosen years are outside the range of EDGAR emissions!'); return;
end

% Define combined sources with separate isotopic ratio
combined_sources = char('Coal/gas/oil/industry','Livestock','Rice','Waste');
combined_dC13    = [     0.5*(-43.7+-44.7),      -65.5,      -62.2, -56.0]; % Zhang2022
combined_dC13_std= [     0.5*(3.8+4.5),          3.3,        3.8,   7.1]; % Zhang2022
ncombined_sources = length(combined_dC13);

% Loop through sources
combined_emis = zeros(ncombined_sources,nyears_file);
for n = 1:nsources
    if(sum(read_data(n,:))>0)
        source = char(read_sources(n));
        
        % Coal/gas/oil/industry
        if(strcmp(source(1),'1')||strcmp(source(1),'2')||strcmp(source(1:3),'5.B'))
            combined_emis(1,:) = combined_emis(1,:) + read_data(n,:);
            
        % Livestock
        elseif(strcmp(source(1:3),'3.A'));
            combined_emis(2,:) = combined_emis(2,:) + read_data(n,:);
            
        % Rice
        elseif(length(source)>4 && strcmp(source(1:5),'3.C.7'))
            combined_emis(3,:) = combined_emis(3,:) + read_data(n,:);
            
        % Waste
        elseif(strcmp(source(1),'4'))
            combined_emis(4,:) = combined_emis(4,:) + read_data(n,:);

            
        % Skip biomass burning emissions (these are included in GFED)
        elseif(length(source)>4 && strcmp(source(1:5),'3.C.1'))
            % skip this since it is included in GFED
            
        else
            n
            source
            disp('ERROR in read_EDGAR.m: Unknown EDGAR source!'); return;
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
    disp(['  - EDGAR (',int2str(eyear),'): ',combined_sources(n,:),': ',num2str(combined_emis_time(n,end)),' Tg yr-1']);
end
