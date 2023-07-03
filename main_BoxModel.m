% ******************************************************
% * Simple tropospheric box model for CH4 and isotopes *
% ******************************************************
% * Author:  Oivind Hodnebrog (oivinho@cicero.oslo.no) *
% * Version: 2023-07-03                                *
% ******************************************************

clear all;


% -------------
% Input options
% -------------

% General options
syear           = 1970; % start year of simulation
eyear           = 2020; % end year of simulation
syear_plot      = 1988; % start year when plotting
eyear_plot      = 2020; % end year when plotting
plot_aerchemmip = 1; % 1= analyze OH from AerChemMIP models
plot_ccmi       = 1; % 1= analyze OH from CCMI models
plot_osloctm3   = 1; % 1= analyze OH from OsloCTM3 model
idealized_OH    = 0; % 0=OH timeseries read from file, 1=constant OH
anthro_emis     = 2; % 1= use CEDS21, 2= use EDGARv7 anthropogenic emission inventory
anthro_const2007= 0; % >0= set chosen anthro category to constant from year 2007 (1=fossil,2=livestock,3=rice,4=waste)
dC13emis_add_std= 0; % e.g., -4= subtract 1-sigma in the isotopic ratio of emissions in category 4 (waste)
nat_emis        = 0; % 0= constant (215 Tg/yr), 1= CLM GSWP3, 2= CLM CRUNCEP, 3= CLM CRUJRA, 10/11= Zhang et al. (2023)
emisBB_clim     = 1; % >0= use GFED climatology (instead of varying) CH4 BB emis. and multiply by given factor (1= no scaling, i.e., use 15.73 Tg/yr)
init_OH         = 0; % 0= use spin-up value, >0= manually set initial OH conc. (molec cm-3)
init_CH4        = 1550; % initial CH4 vmr (ppb)
save_run2matfile= 0; % 1= save all variables to a .mat file (note that all figures will be closed)

% OH reaction and kinetic isotope effect (KIE), alpha=k13/k12
troptemp = 268.5; % temperature in the troposphere (268.5 K gives tau_CH4 wrt. OH of ~9.6 yrs as in CH4 box model if I assume OH conc. of 1e6 molec/cm3)
alpha_OH = 0.9946; % Cantrell et al., 1990
k_CH4_OH = 2.45e-12*exp(-1775/troptemp); % reaction rate of OH with CH4 (cm3/molec/s) - from NASA JPL 2019

% Cl reaction, abundance and KIE
const_Cl = 620; % tropospheric mean Cl conc. (molecules cm-3) from Wang et al. (2019, ACP)
alpha_Cl = 0.938; % KIE for Cl + CH4, see Strode et al. (2020, ACP)
k_CH4_Cl = 7.1e-12*exp(-1270/troptemp); % reaction rate of Cl with CH4 (cm3/molec/s) - from NASA JPL 2019

% Soil loss and KIE
alpha_soil = 0.978; % KIE for CH4 soil sink, from Tyler et al. (1994) - see Strode et al.
k_CH4_soil = 1/(160*3600*24*365); % 1/s (use 160 yrs as in CH4 box model)

% Plotting options
plot_dC13CH4_anomaly = 0; % 1= plot the anomaly of d13C(CH4) instead of raw data
plot_dC13CH4_twoyaxes= 1; % 1= use two y-axes (one for obs and one for models)
plot_emis = 1; % 1= plot CH4 emissions and their delta13C composition
plot_budget = 0; % 1= plot chemical budget for CH4 (only using OsloCTM3 OH)
plot_spinup = 0; % 1= plot timeseries of delta13C during spin-up


% ----------------------------
% Paths, data, constants, etc.
% ----------------------------
addpath('./functions');

% OH datasets
file_OH_aerchemmip      = './inputdata/AerChemMIP_modelmean_OH.txt';
file_OH_aerchemmip_waccm= './inputdata/AerChemMIP_CESM2-WACCM_OH.txt';
file_OH_aerchemmip_gfdl = './inputdata/AerChemMIP_GFDL_OH.txt';
file_OH_aerchemmip_uk   = './inputdata/AerChemMIP_UK_OH.txt';
file_OH_ccmi            = './inputdata/CCMI_modelmean_OH.txt';
file_OH_osloctm3        = './inputdata/OsloCTM3_OH_histO3_ceds2021.txt';

% CH4 datasets
file_CH4_SSP  = './inputdata/WMGHG_vmr_SSP2-4.5';
file_CH4_NOAA = './inputdata/NOAA_MoleFractions_2022.csv';
file_CEDS     = './inputdata/CH4_global_CEDS_emissions_by_sector_2021_04_21.csv';
file_EDGAR    = './inputdata/EDGAR_CH4_1970-2021.csv';
file_GFED     = './inputdata/anthropogenic_emissions_gfed.txt';
file_CLM      = './inputdata/natemis.csv';

% d13C(CH4) observations, from Table S4 in Schaefer et al. (2016)
years_dC13CH4_obs= 1988:2014;
data_dC13CH4_obs = [-47.43,-47.40,-47.37,-47.40,-47.35,-47.37,-47.31,-47.27,-47.31,-47.23,-47.19,-47.20,-47.21,-47.19,...
                    -47.21,-47.19,-47.18,-47.22,-47.20,-47.21,-47.17,-47.25,-47.31,-47.35,-47.35,-47.37,-47.37];

% d13C(CH4) observations, from WMO/GAW GHG bulletin 2022
years_dC13CH4_obs2 = 1999:2020; nyears_dC13CH4_obs2 = length(years_dC13CH4_obs2);
data_dC13CH4_obs2 = nan(nyears_dC13CH4_obs2,1); % * Note that values are set to NaN because data are not openly available
data_dC13CH4_obs2_unc = nan(nyears_dC13CH4_obs2,1); % but received from Sylvia Englund Michel (personal communication) *
data_dC13CH4_obs2_max = data_dC13CH4_obs2+data_dC13CH4_obs2_unc;
data_dC13CH4_obs2_min = data_dC13CH4_obs2-data_dC13CH4_obs2_unc;

% Filename for storing results from box model simulation (as .mat file)
save_runname = ['run_emisBBclim',int2str(emisBB_clim),'_idealizedOH',int2str(idealized_OH),'_initOH',int2str(init_OH),...
                '_anthroemis',int2str(anthro_emis),'_anthroconst',int2str(anthro_const2007),'_natemis',int2str(nat_emis),...
                '_dC13emisaddstd',int2str(dC13emis_add_std)];

% Constants
radius_earth = 6371e3; % m
mass_atm = 5.1352e18; % kg, dry air mass of atmosphere (Trenberth et al., 2005)
molarmass_air = 28.9647; % g/mol, dry air
molarmass_CH4 = 16.04; % g/mol
avogadro   = 6.022e23; % molecules/mol

% Calculate volume of troposphere, very approximately assuming tropopause height of 13 km
volume_troposphere = (4/3*pi*(radius_earth+13e3)^3) - (4/3*pi*radius_earth^3); % m

% Vienna PDB (VPDB) isotopic standard (Zhang and Li, 1990)
R_std = 0.011183;

% Time definitions
dt = 3600*24*365; % seconds in a year
years = [syear:eyear];
nyears= length(years);

% Set Cl concentrations to constant
data_Cl = nan(nyears,1); data_Cl(:) = const_Cl;


% --------------------------------------------------------------------
% Read time series of CH4 emissions (Tg yr-1) and their isotopic ratio
% --------------------------------------------------------------------

disp('- Reading CH4 emissions');
s = 0;
    
% Read anthropogenic emission
if(anthro_emis==1) % CEDS21
    disp('- Will read anthropogenic emissions from CEDS21 inventory');
    [sources_CH4,s] = read_CEDS(0,s,syear,eyear,file_CEDS);
elseif(anthro_emis==2) % EDGARv7
    disp('- Will read anthropogenic emissions from EDGARv7 inventory');
    [sources_CH4,s] = read_EDGAR(0,s,syear,eyear,file_EDGAR);
else
    disp('ERROR: Unknown option for anthro_emis!'); return;
end

% Set anthropogenic emission category to constant from year 2007 onwards, if relevant
if(anthro_const2007>0 && anthro_const2007<5)
    if(syear>2007 || eyear<=2007)
        disp('ERROR: anthro_const2007 option requires time period to include year 2007!'); return;
    end
    disp(['WARNING: ',sources_CH4(anthro_const2007).name,' emissions set to constant from year 2007!']);
    sources_CH4(anthro_const2007).emissions(2007-syear+1:end) = sources_CH4(anthro_const2007).emissions(2007-syear+1);
end

% Read biomass burning and natural emissions
[sources_CH4,s] = read_GFED(sources_CH4,s,syear,eyear,file_GFED,emisBB_clim);
[sources_CH4,s,~] = read_natural(sources_CH4,s,syear,eyear,nat_emis,file_CLM);
nsources_CH4 = s; nsources_CH4_ant = 5;

% Sum up total CH4 emissions
emissions_CH4 = zeros(nyears,1);
for n = 1:nsources_CH4; emissions_CH4(:) = emissions_CH4(:) + sources_CH4(n).emissions(:); end
disp(['  - Total CH4 emissions (',int2str(eyear),'): ',num2str(emissions_CH4(nyears)),' Tg yr-1']);

% Add or subtract 1-sigma of the dC13 signature of chosen emission category
s = abs(dC13emis_add_std);
if(dC13emis_add_std>0) % add 1-sigma
    sources_CH4(s).dC13 = sources_CH4(s).dC13 + sources_CH4(s).dC13_std;
    disp(['WARNING: Adding 1-sigma in the dC13 signature of ',sources_CH4(s).name,' emissions!']);
elseif(dC13emis_add_std<0) % subtract 1-sigma
    sources_CH4(s).dC13 = sources_CH4(s).dC13 - sources_CH4(s).dC13_std;
    disp(['WARNING: Subtracting 1-sigma in the dC13 signature of ',sources_CH4(s).name,' emissions!']);
end
    
% Calculate a weighted average delta13 value (this is the same as weighing the C13/C12 ratios)
dC13_S_CH4 = zeros(nyears,1);
for y = 1:nyears
    for s = 1:nsources_CH4
        dC13_S_CH4(y) = dC13_S_CH4(y) + sources_CH4(s).dC13*sources_CH4(s).emissions(y)/emissions_CH4(y);
    end
end
disp(['- Average d13C value for CH4 emissions for year ',int2str(eyear),' is ',num2str(dC13_S_CH4(end)),' per mil']);
    
% Calculate CH4 source term for each year
for y = 1:nyears
    S_CH4(y) = emissions_CH4(y)*1e12/molarmass_CH4*avogadro/(volume_troposphere*1e6)/dt; % convert Tg/yr -> molec/cm3/s
end


% ----------------------------------------------------------
% Read CH4 and OH concentration, and spin up OH and dC13_CH4
% ----------------------------------------------------------

% Read CH4 concentrations
[data_CH4_obs] = read_CH4(file_CH4_SSP,file_CH4_NOAA,syear,eyear);

% Crudely convert CH4 from ppb to molec/cm3
mass_troposphere = 0.75*mass_atm; % kg, assume 75% of atmospheric mass below tropopause
n_air = mass_troposphere*1e3/molarmass_air*avogadro/(volume_troposphere*1e6); % molec(air)/cm3, number dens. of air in trop.
init_n_CH4 = init_CH4*1e-9*n_air; % ppb -> molec(CH4)/cm3

% Spin up OH and dC13_CH4 for a given initial concentration and emission of CH4
[init_OH_spinup,init_dC13CH4] = calc_initvals(init_n_CH4,S_CH4(1),dC13_S_CH4(1),R_std,n_air,plot_spinup,...
                                              k_CH4_OH,alpha_OH,...
                                              k_CH4_Cl,alpha_Cl,data_Cl(1),...
                                              k_CH4_soil,alpha_soil);
disp(['- Initial dC13CH4 value from spin-up is ',num2str(init_dC13CH4),' permil']);

% Use manually set OH concentration or the one from spin-up
if(init_OH>0)
    disp(['WARNING: Will use manually set OH conc. of ',num2str(init_OH),' molec/cm3!']);
else
    disp(['- Will use spinned up OH conc. of ',num2str(init_OH_spinup),' molec/cm3']);
    init_OH = init_OH_spinup;
end

% Read OH values from file
data_OH_aerchemmip = nan(nyears,1);
data_OH_aerchemmip_waccm = nan(nyears,1);
data_OH_aerchemmip_gfdl = nan(nyears,1);
data_OH_aerchemmip_uk = nan(nyears,1);
data_OH_ccmi = nan(nyears,1);
data_OH_osloctm3 = nan(nyears,1);
if(plot_aerchemmip)
    [data_OH_aerchemmip,eyear_OH_aerchemmip,~] = read_OH(file_OH_aerchemmip,syear,eyear,init_OH);
    [data_OH_aerchemmip_waccm,~,~] = read_OH(file_OH_aerchemmip_waccm,syear,eyear,init_OH);
    [data_OH_aerchemmip_gfdl,~,~] = read_OH(file_OH_aerchemmip_gfdl,syear,eyear,init_OH);
    [data_OH_aerchemmip_uk,~,~] = read_OH(file_OH_aerchemmip_uk,syear,eyear,init_OH);
end
if(plot_ccmi); [data_OH_ccmi,eyear_OH_ccmi,syear_OH_ccmi] = read_OH(file_OH_ccmi,syear,eyear,init_OH); end
if(plot_osloctm3); [data_OH_osloctm3,eyear_OH_osloctm3,~] = read_OH(file_OH_osloctm3,syear,eyear,init_OH); end

% Do an idealized run with constant OH instead of using the OH timeseries
if(idealized_OH==1)
    disp('* WARNING: Replacing OH time series with a constant value! *');
    data_OH_aerchemmip(:) = init_OH;
    data_OH_aerchemmip_waccm(:) = init_OH;
    data_OH_aerchemmip_gfdl(:) = init_OH;
    data_OH_aerchemmip_uk(:) = init_OH;
    data_OH_ccmi(:) = init_OH;
    data_OH_osloctm3(:) = init_OH;
end


% ------------------------------------------------------------------------------
% Run the box model by looping through the years and calculating CH4 and dC13CH4
% ------------------------------------------------------------------------------

if(plot_aerchemmip)
    [data_dC13CH4_aerchemmip,data_CH4_aerchemmip,~,~,~,~] = ...
        calc_CH4_evolution(data_OH_aerchemmip,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'AerChemMIP');
    [data_dC13CH4_aerchemmip_waccm,data_CH4_aerchemmip_waccm,~,~,~,~] = ...
        calc_CH4_evolution(data_OH_aerchemmip_waccm,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'AerChemMIP_WACCM');
    [data_dC13CH4_aerchemmip_gfdl,data_CH4_aerchemmip_gfdl,~,~,~,~] = ...
        calc_CH4_evolution(data_OH_aerchemmip_gfdl,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'AerChemMIP_GFDL');
    [data_dC13CH4_aerchemmip_uk,data_CH4_aerchemmip_uk,~,~,~,~] = ...
        calc_CH4_evolution(data_OH_aerchemmip_uk,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'AerChemMIP_UK');

    % Find the range (min/max) based on the AerChemMIP models
    data_OH_aerchemmip_min = min(data_OH_aerchemmip_waccm,min(data_OH_aerchemmip_gfdl,data_OH_aerchemmip_uk));
    data_OH_aerchemmip_max = max(data_OH_aerchemmip_waccm,max(data_OH_aerchemmip_gfdl,data_OH_aerchemmip_uk));
    data_CH4_aerchemmip_min = min(data_CH4_aerchemmip_waccm,min(data_CH4_aerchemmip_gfdl,data_CH4_aerchemmip_uk));
    data_CH4_aerchemmip_max = max(data_CH4_aerchemmip_waccm,max(data_CH4_aerchemmip_gfdl,data_CH4_aerchemmip_uk));
    data_dC13CH4_aerchemmip_min = min(data_dC13CH4_aerchemmip_waccm,min(data_dC13CH4_aerchemmip_gfdl,data_dC13CH4_aerchemmip_uk));
    data_dC13CH4_aerchemmip_max = max(data_dC13CH4_aerchemmip_waccm,max(data_dC13CH4_aerchemmip_gfdl,data_dC13CH4_aerchemmip_uk));
end
if(plot_ccmi)
    [data_dC13CH4_ccmi,data_CH4_ccmi,~,~,~,~] = calc_CH4_evolution(data_OH_ccmi,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'CCMI');
end
if(plot_osloctm3)
    [data_dC13CH4_osloctm3,data_CH4_osloctm3,budget_S_CH4_osloctm3,...
     budget_l_CH4_OH_osloctm3,budget_l_CH4_Cl_osloctm3,budget_l_CH4_soil_osloctm3] = ...
        calc_CH4_evolution(data_OH_osloctm3,init_n_CH4,S_CH4,init_dC13CH4,dC13_S_CH4,...
                           k_CH4_OH,alpha_OH,k_CH4_Cl,alpha_Cl,data_Cl,...
                           k_CH4_soil,alpha_soil,data_CH4_obs,...
                           R_std,n_air,years,'OsloCTM3');

    % Convert budget unit, molec/cm3/s -> Tg yr-1
    budget_S_CH4_osloctm3 = budget_S_CH4_osloctm3*1e-12*molarmass_CH4/avogadro*volume_troposphere*1e6*dt;
    budget_l_CH4_OH_osloctm3 = budget_l_CH4_OH_osloctm3*1e-12*molarmass_CH4/avogadro*volume_troposphere*1e6*dt;
    budget_l_CH4_Cl_osloctm3 = budget_l_CH4_Cl_osloctm3*1e-12*molarmass_CH4/avogadro*volume_troposphere*1e6*dt;
    budget_l_CH4_soil_osloctm3 = budget_l_CH4_soil_osloctm3*1e-12*molarmass_CH4/avogadro*volume_troposphere*1e6*dt;
end

% Calculate anomaly, if needed
if(plot_dC13CH4_anomaly)
    if(plot_aerchemmip)
        data_dC13CH4_aerchemmip_min = data_dC13CH4_aerchemmip_min - mean(data_dC13CH4_aerchemmip);
        data_dC13CH4_aerchemmip_max = data_dC13CH4_aerchemmip_max - mean(data_dC13CH4_aerchemmip);
        data_dC13CH4_aerchemmip = data_dC13CH4_aerchemmip - mean(data_dC13CH4_aerchemmip);
    end
    if(plot_ccmi); data_dC13CH4_ccmi = data_dC13CH4_ccmi - mean(data_dC13CH4_ccmi); end
    if(plot_osloctm3); data_dC13CH4_osloctm3 = data_dC13CH4_osloctm3 - mean(data_dC13CH4_osloctm3); end
    
    % Normalize the dC13CH4 observation time series to give anomalies with reference to the common time period of Schaefer et al. and WMO/GAW
    syear_obs = max(years_dC13CH4_obs(1),years_dC13CH4_obs2(1));
    eyear_obs = min(years_dC13CH4_obs(end),years_dC13CH4_obs2(end));
    data_dC13CH4_obs = data_dC13CH4_obs - mean(data_dC13CH4_obs(syear_obs-years_dC13CH4_obs(1)+1:eyear_obs-years_dC13CH4_obs(1)+1));
    data_dC13CH4_obs2_min = data_dC13CH4_obs2_min - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
    data_dC13CH4_obs2_max = data_dC13CH4_obs2_max - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
    data_dC13CH4_obs2= data_dC13CH4_obs2 - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
end


% ---------------------------------------
% Plot time evolution of OH, CH4 and d13C
% ---------------------------------------
p = figure; set(p,'PaperPositionMode','auto'); set(p,'Position',[0 0 1000 1200]);

% Colors
c_aerchemmip = [77,175,74]/255;
c_ccmi = [55,126,184]/255;
c_osloctm3 = [228,26,28]/255;
if(idealized_OH==1); c_osloctm3 = [106,61,154]/255; end

% Plot OH
subplot(3,1,1); hold on;
lgd = legend('FontSize',8,'Orientation','horizontal','Location','NorthWest');%'SouthEast');
if(plot_aerchemmip)
    e_index_aerchemmip = min(eyear,eyear_OH_aerchemmip)-syear+1;
    set(lgd,'AutoUpdate','off');
    patch([years(1:e_index_aerchemmip) fliplr(years(1:e_index_aerchemmip))],[data_OH_aerchemmip_max(1:e_index_aerchemmip)' fliplr(data_OH_aerchemmip_min(1:e_index_aerchemmip)')],c_aerchemmip,'EdgeColor','none','FaceAlpha',0.333);
    if(idealized_OH==0); set(lgd,'AutoUpdate','on'); end
    plot(years(1:eyear_OH_aerchemmip-syear+1),data_OH_aerchemmip(1:eyear_OH_aerchemmip-syear+1),'Color',c_aerchemmip,'LineWidth',2,'DisplayName','AerChemMIP');
    if(eyear_OH_aerchemmip<eyear) % then do the latter years with constant OH, if relevant
        set(lgd,'AutoUpdate','off');
        plot(years(eyear_OH_aerchemmip-syear+1:nyears),data_OH_aerchemmip(eyear_OH_aerchemmip-syear+1:nyears),':','Color',c_aerchemmip,'LineWidth',2);
        set(lgd,'AutoUpdate','on');
    end
end
if(plot_ccmi)
    s_index_ccmi = max(syear,syear_OH_ccmi)-syear+1;
    e_index_ccmi = min(eyear,eyear_OH_ccmi)-syear+1;
    if(idealized_OH>0); set(lgd,'AutoUpdate','off'); end
    plot(years(1:eyear_OH_ccmi-syear+1),data_OH_ccmi(1:eyear_OH_ccmi-syear+1),'Color',c_ccmi,'LineWidth',2,'DisplayName','CCMI');
    if(eyear_OH_ccmi<eyear) % then do the latter years with constant OH, if relevant
        set(lgd,'AutoUpdate','off');
        plot(years(eyear_OH_ccmi-syear+1:nyears),data_OH_ccmi(eyear_OH_ccmi-syear+1:nyears),':','Color',c_ccmi,'LineWidth',2);
        set(lgd,'AutoUpdate','on');
    end
end
if(plot_osloctm3)
    osloctm3_label = 'OsloCTM3 (CEDS21+COVID)';
    if(idealized_OH==1); osloctm3_label = 'Constant'; end
    set(lgd,'AutoUpdate','on');
    plot(years(1:eyear_OH_osloctm3-syear+1),data_OH_osloctm3(1:eyear_OH_osloctm3-syear+1),'Color',c_osloctm3,'LineWidth',2,'DisplayName',osloctm3_label);
    if(eyear_OH_osloctm3<eyear) % then do the latter years with constant OH, if relevant
        set(lgd,'AutoUpdate','off');
        plot(years(eyear_OH_osloctm3-syear+1:nyears),data_OH_osloctm3(eyear_OH_osloctm3-syear+1:nyears),':','Color',c_osloctm3,'LineWidth',2);
        set(lgd,'AutoUpdate','on');
    end
end
set(lgd,'AutoUpdate','off');
hold off;
xlim([syear_plot eyear_plot]);
grid on;
title('OH concentrations'); ylabel('molec cm^{-3}');

% Plot CH4
subplot(3,1,2); hold on;
if(plot_aerchemmip)
    patch([years(1:e_index_aerchemmip) fliplr(years(1:e_index_aerchemmip))],[data_CH4_aerchemmip_max(1:e_index_aerchemmip)' fliplr(data_CH4_aerchemmip_min(1:e_index_aerchemmip)')],c_aerchemmip,'EdgeColor','none','FaceAlpha',0.333);
    plot(years(1:eyear_OH_aerchemmip-syear+1),data_CH4_aerchemmip(1:eyear_OH_aerchemmip-syear+1),'Color',c_aerchemmip,'LineWidth',2);
    if(eyear_OH_aerchemmip<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_aerchemmip-syear+1:nyears),data_CH4_aerchemmip(eyear_OH_aerchemmip-syear+1:nyears),':','Color',c_aerchemmip,'LineWidth',2);
    end
end
if(plot_ccmi)
    plot(years(1:eyear_OH_ccmi-syear+1),data_CH4_ccmi(1:eyear_OH_ccmi-syear+1),'Color',c_ccmi,'LineWidth',2);
    if(eyear_OH_ccmi<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_ccmi-syear+1:nyears),data_CH4_ccmi(eyear_OH_ccmi-syear+1:nyears),':','Color',c_ccmi,'LineWidth',2);
    end
end
if(plot_osloctm3)
    plot(years(1:eyear_OH_osloctm3-syear+1),data_CH4_osloctm3(1:eyear_OH_osloctm3-syear+1),'Color',c_osloctm3,'LineWidth',2);
    if(eyear_OH_osloctm3<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_osloctm3-syear+1:nyears),data_CH4_osloctm3(eyear_OH_osloctm3-syear+1:nyears),':','Color',c_osloctm3,'LineWidth',2);
    end
end
plot(years,data_CH4_obs,'k','LineWidth',2);
hold off;
xlim([syear_plot eyear_plot]);
grid on;
title('CH_4 mixing ratios'); ylabel('ppb');

% Plot dC13CH4
subplot(3,1,3); hold on;
if(plot_dC13CH4_twoyaxes); yyaxis left; set(gca,'YColor','k'); end
patch([years_dC13CH4_obs2 fliplr(years_dC13CH4_obs2)],[data_dC13CH4_obs2_max fliplr(data_dC13CH4_obs2_min)],[0.25,0.25,0.25],'EdgeColor','none','FaceAlpha',0.333);
plot(years_dC13CH4_obs2,data_dC13CH4_obs2,'-','Color',[0.25,0.25,0.25],'LineWidth',2,'DisplayName','Obs. WMO/GAW');
plot(years_dC13CH4_obs,data_dC13CH4_obs,'-k','LineWidth',2,'DisplayName','Obs. Schaefer16');
if(plot_dC13CH4_twoyaxes)
    ylabel(char(8240));
    yyaxis right; set(gca,'YColor',c_osloctm3);
end
if(plot_aerchemmip)
    patch([years(1:e_index_aerchemmip) fliplr(years(1:e_index_aerchemmip))],[data_dC13CH4_aerchemmip_max(1:e_index_aerchemmip)' fliplr(data_dC13CH4_aerchemmip_min(1:e_index_aerchemmip)')],c_aerchemmip,'EdgeColor','none','FaceAlpha',0.333);
    plot(years(1:eyear_OH_aerchemmip-syear+1),data_dC13CH4_aerchemmip(1:eyear_OH_aerchemmip-syear+1),'-','Color',c_aerchemmip,'LineWidth',2);
    if(eyear_OH_aerchemmip<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_aerchemmip-syear+1:nyears),data_dC13CH4_aerchemmip(eyear_OH_aerchemmip-syear+1:nyears),':','Color',c_aerchemmip,'LineWidth',2);
    end
end
if(plot_ccmi)
    plot(years(1:eyear_OH_ccmi-syear+1),data_dC13CH4_ccmi(1:eyear_OH_ccmi-syear+1),'-','Color',c_ccmi,'LineWidth',2);
    if(eyear_OH_ccmi<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_ccmi-syear+1:nyears),data_dC13CH4_ccmi(eyear_OH_ccmi-syear+1:nyears),':','Color',c_ccmi,'LineWidth',2);
    end
end
if(plot_osloctm3)
    plot(years(1:eyear_OH_osloctm3-syear+1),data_dC13CH4_osloctm3(1:eyear_OH_osloctm3-syear+1),'-','Color',c_osloctm3,'LineWidth',2);
    if(eyear_OH_osloctm3<eyear) % then do the latter years with constant OH, if relevant
        plot(years(eyear_OH_osloctm3-syear+1:nyears),data_dC13CH4_osloctm3(eyear_OH_osloctm3-syear+1:nyears),':','Color',c_osloctm3,'LineWidth',2);
    end
end
hold off;
xlim([syear_plot eyear_plot]);
grid on;
if(plot_dC13CH4_anomaly); title('\Delta \delta^{13}C-CH_4'); else title('\delta^{13}C-CH_4'); end
ylabel(char(8240));
set(gcf,'Color','w');
savename = 'OH_CH4_dC13CH4';
saveas(gcf,['./plots/',savename,'.fig']);
print(gcf,['./plots/',savename,'.png'],'-dpng','-r300');


% -------------------------------------------
% Plot CH4 emissions and their isotopic ratio
% -------------------------------------------
if(plot_emis)
    colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
              0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330;...
              0.6350 0.0780 0.1840];
    colors(1:4,:) = [217,95,2; 117,112,179; 231,41,138; 230,171,2]/255; % overwrite colors for anthropogenic sectors
    colors(6,:) = [31,120,180]/255; % overwrite color for wetland
    p = figure; set(p,'PaperPositionMode','auto'); set(p,'Position',[0 0 1200 900]);

    % Plot time evolution of CH4 emissions
    subplot(2,1,1); hold on;
    yyaxis left; set(gca,'YColor','k');
    for n = 1:4 %nsources_CH4
        linetype = '-'; icolor=n;
        if(n>nsources_CH4_ant); linetype=':'; icolor=n-nsources_CH4_ant;end
        plot([syear:eyear],sources_CH4(n).emissions,linetype,'Color',colors(icolor,:),'LineWidth',2,'DisplayName',sources_CH4(n).name);
    end
    legend('Location','EastOutside'); hold off; ylim([0 200]); ylabel(['Tg yr^{-1}']);
    grid on;
    title('CH_4 emissions');
    yyaxis right; set(gca,'YColor',[0.5,0.5,0.5]);
    plot([syear:eyear],dC13_S_CH4,'-','Color',[0.5,0.5,0.5],'LineWidth',2,'DisplayName','\delta^{13}C-CH_4 signature'); ylabel(char(8240));

    % Plot bar chart with isotopic signature and uncertainty
    subplot(2,1,2); hold on;
    for n = 1:6 %nsources_CH4
        linetype = '-'; icolor=n;
        bar(n,sources_CH4(n).dC13,'FaceColor',colors(icolor,:),'LineStyle','none','DisplayName',sources_CH4(n).name);
    end
    lgd = legend('Location','EastOutside');
    set(lgd,'AutoUpdate','off');
    for n = 1:6 %nsources_CH4
        errorbar(n,sources_CH4(n).dC13,sources_CH4(n).dC13_std,sources_CH4(n).dC13_std,'Color',[0,0,0],'LineWidth',2);%,'LineStyle',linetype,'DisplayName',sources_CH4(n).name);
    end
    hold off;
    title('\delta^{13}C-CH_4 isotopic signature of emissions'); ylabel(char(8240));
    xlim([0.5 6.5]);
    grid on;
    set(gcf,'Color','w');
    set(gca,'xtick',0);
    set(gca,'xticklabel',{[]});
    savename = 'CH4_emissions';
    saveas(gcf,['./plots/',savename,'.fig']);
    print(gcf,['./plots/',savename,'.png'],'-dpng','-r300');
end


% ------------------------------------
% Plot chemical budget (OsloCTM3 only)
% ------------------------------------
if(plot_budget && plot_osloctm3)
    figure; hold on;
    plot(years(2:nyears),budget_S_CH4_osloctm3(2:nyears),'LineWidth',2,'DisplayName','Emissions');
    plot(years(2:nyears),budget_l_CH4_OH_osloctm3(2:nyears),'LineWidth',2,'DisplayName','Loss through OH reaction');
    plot(years(2:nyears),budget_l_CH4_Cl_osloctm3(2:nyears),'LineWidth',2,'DisplayName','Loss through Cl reaction');
    plot(years(2:nyears),budget_l_CH4_soil_osloctm3(2:nyears),'LineWidth',2,'DisplayName','Loss through soil uptake');
    set(gcf,'Color','w'); ylabel(['Tg yr^{-1}']); legend('Location','East');
    title('Chemical budget of CH_4');
end


% ---------------------------------
% Save all variables to a .mat file
% ---------------------------------
if(save_run2matfile)
    disp('- Closing all figures before saving run to .mat file');
    close all; % close all figures
    save_runname_version = 2;
    save_runname_file = ['./matfiles/',save_runname,'.mat'];
    disp(['NOTE: Saving to mat file ',save_runname_file]);
    save(save_runname_file);
end
