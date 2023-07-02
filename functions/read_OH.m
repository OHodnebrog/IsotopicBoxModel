% This function reads timeseries of OH from a given file
% ******************************************************

function [ data_OH,eyear_OH,syear_OH ] = read_OH( file_OH,syear,eyear,init_OH )

nyears = eyear-syear+1;

% Read OH values from file
data_read = importdata(file_OH,',',1);
if(syear>data_read.data(end,1))
    disp('ERROR in read_OH: Check syear!'); return;
end
syear_OH = syear;
if(syear<data_read.data(1,1))
    syear_OH = data_read.data(1,1);
end
index_start = syear_OH-data_read.data(1,1)+1;
index_end   = syear-data_read.data(1,1)+nyears;
index_max = length(squeeze(data_read.data(:,1)));
if(index_end>index_max)
    index_end = index_max;
    eyear_OH = data_read.data(end,1);
else
    eyear_OH = eyear;
end
data_dOH = data_read.data(index_start:index_end,2); clear data_read;
disp(['- OH change dataset for years ',int2str(syear_OH),'-',int2str(eyear_OH),' from file ',file_OH]);

% Convert relative changes in OH to actual OH concentrations
data_OH = nan(nyears,1);
for y = syear_OH:eyear_OH
    data_OH(y-syear+1) = data_dOH(y-syear_OH+1)/data_dOH(1)*init_OH;
end

% Set OH concentration to constant for the latest years, if necessary
if(eyear_OH<eyear)
    data_OH(eyear_OH-syear+2:nyears) = data_OH(eyear_OH-syear+1);
end
