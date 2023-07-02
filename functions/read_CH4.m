% This function reads timeseries of CH4 (in ppb) from a given file
% ****************************************************************

function [ data_CH4 ] = read_CH4( file_CH4_SSP,file_CH4_NOAA,syear,eyear )

nyears = eyear-syear+1;

% Read SSP file based on Meinshausen et al. (2017, GMD) first
data_read = importdata(file_CH4_SSP,' ',2);
if(syear<data_read.data(1,1)||syear>data_read.data(end,1))
    disp('ERROR in read_CH4.m: Check syear!'); return;
end
index_start = syear-data_read.data(1,1)+1;
index_end   = index_start+nyears-1;
index_max = length(squeeze(data_read.data(:,1)));
if(index_end>index_max)
    index_end = index_max;
    disp(['ERROR in read_CH4.m: CH4 dataset ends before ',int2str(eyear),'!']);
    return;
end

% Then read observed CH4 vmr based on NOAA
data_read_noaa = importdata(file_CH4_NOAA,',',4);
syear_noaa = floor(data_read_noaa.data(1,1));
eyear_noaa = floor(data_read_noaa.data(end,1));
data_noaa_CH4 = data_read_noaa.data(:,3);

% Overwrite time series with observed dataset for years that overlap
data_combined = data_read.data(:,4);
data_combined(syear_noaa-data_read.data(1,1)+1:eyear_noaa-data_read.data(1,1)+1)=data_noaa_CH4(:);

data_CH4 = data_combined(index_start:index_end);
disp(['- CH4 dataset read from file for years ',int2str(syear),'-',int2str(eyear)]);
