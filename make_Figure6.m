
% Script to make plot based on combining data from different box model runs.
% The script requires main_BoxModel.m to be run first using save_run2matfile=1.
% *****************************************************************************
clear all;


% Input options
INP_emisBB_clim     = 1; % >0= use GFED climatology (instead of varying) CH4 BB emis. and multiply by given factor (1= no scaling, i.e., use 15.73 Tg/yr)
INP_idealized_OH    = 0; % 0=OH timeseries read from file, 1=constant OH
INP_init_OH         = 0; % 0= use spin-up value, >0= manually set initial OH conc. (molec cm-3)
INP_anthro_emis     = 2; % 1= use CEDS21, 2= use EDGARv7 anthropogenic emission inventory
INP_anthro_const2007= 0; % >0= set chosen anthro category to constant from year 2007 (1=fossil,2=livestock,3=rice,4=waste)
INP_nat_emis        = 0; % 0= constant (215 Tg/yr), 1= CLM GSWP3, 2= CLM CRUNCEP, 3= CLM CRUJRA, 10/11= Zhang et al. (2023)
INP_inc_dC13emis_unc= 0; % e.g., 6= plot +/-1-sigma in the isotopic ratio of emissions in category 6 (wetlands) as shading


% Do some checking
if(INP_inc_dC13emis_unc>0 && INP_anthro_const2007==0 && INP_nat_emis==0)
    disp('ERROR: Need to have either INP_anthro_const2007>0 or INP_nat_emis>0 when INP_inc_dC13emis_unc>0!');
    return;
end
if((INP_idealized_OH>0 && INP_anthro_const2007>0) || (INP_idealized_OH>0 && INP_nat_emis>0) || (INP_anthro_const2007>0 && INP_nat_emis>0))
    disp('ERROR: Need to have only one of the parameters INP_idealized_OH, INP_anthro_const2007 or INP_nat_emis be >0!');
    return;
end
if(INP_nat_emis>0)
    disp('NOTE: Difference will be taken as run2-run1 instead of run1-run2!');
    reverse_diff = -1;
else
    reverse_diff = 1;
end

% Load data
run1 = ['./matfiles/run_emisBBclim',int2str(INP_emisBB_clim),'_idealizedOH0_initOH',int2str(INP_init_OH),'_anthroemis',int2str(INP_anthro_emis),'_anthroconst0_natemis0_dC13emisaddstd0.mat'];
run2 = ['./matfiles/run_emisBBclim',int2str(INP_emisBB_clim),'_idealizedOH',int2str(INP_idealized_OH),'_initOH',int2str(INP_init_OH),'_anthroemis',int2str(INP_anthro_emis),'_anthroconst',int2str(INP_anthro_const2007),'_natemis',int2str(INP_nat_emis),'_dC13emisaddstd0.mat'];
load(run1);            % run1 data will be the default
run2data = load(run2); % run2 data will be loaded into a struct
if(plot_dC13CH4_anomaly==1 || run2data.plot_dC13CH4_anomaly==1)
    disp('ERROR: Need to have plot_dC13CH4_anomaly=0 in input data!'); return;
end

% Load additional mat files with +/-1-sigma uncertainty
if(INP_inc_dC13emis_unc>0)
    run1data_high = load(strrep(run1,'std0',['std',int2str(INP_inc_dC13emis_unc)]));
    run1data_low  = load(strrep(run1,'std0',['std-',int2str(INP_inc_dC13emis_unc)]));
    run2data_high = load(strrep(run2,'std0',['std',int2str(INP_inc_dC13emis_unc)]));
    run2data_low  = load(strrep(run2,'std0',['std-',int2str(INP_inc_dC13emis_unc)]));
    if(run1data_high.plot_dC13CH4_anomaly==1 || run1data_low.plot_dC13CH4_anomaly==1 || run2data_high.plot_dC13CH4_anomaly==1 || run2data_low.plot_dC13CH4_anomaly==1)
        disp('ERROR: Need to have plot_dC13CH4_anomaly=0 in input data!'); return;
    end
end


% Overwrite plotting options
plot_dC13CH4_anomaly = 1; % 1= plot the anomaly of d13C(CH4) instead of raw data for the observations
set_maxmin_d13C= 1;
if(INP_idealized_OH>0); minval_d13C= -0.4; maxval_d13C= 0.15; end % a) OH plot
if(INP_anthro_const2007>0); minval_d13C= -0.4; maxval_d13C= 0.3; end % b) anthro emis plot
if(INP_nat_emis>0); minval_d13C= -0.4; maxval_d13C= 0.25; end % c) natural emis plot
if(INP_anthro_const2007>0 || INP_nat_emis>0) % only plot OsloCTM3 for the emission perturbations
    plot_aerchemmip = 0; plot_ccmi = 0;
end

% Take the difference in dC13CH4 between run 1 and 2
if(plot_aerchemmip)
    data_dC13CH4_aerchemmip_waccm = data_dC13CH4_aerchemmip_waccm-run2data.data_dC13CH4_aerchemmip_waccm;
    data_dC13CH4_aerchemmip_gfdl = data_dC13CH4_aerchemmip_gfdl-run2data.data_dC13CH4_aerchemmip_gfdl;
    data_dC13CH4_aerchemmip_uk = data_dC13CH4_aerchemmip_uk-run2data.data_dC13CH4_aerchemmip_uk;
    data_dC13CH4_aerchemmip = data_dC13CH4_aerchemmip-run2data.data_dC13CH4_aerchemmip;
    
    % Find the range (min/max) based on the AerChemMIP models
    data_dC13CH4_aerchemmip_min = min(data_dC13CH4_aerchemmip_waccm,min(data_dC13CH4_aerchemmip_gfdl,data_dC13CH4_aerchemmip_uk));
    data_dC13CH4_aerchemmip_max = max(data_dC13CH4_aerchemmip_waccm,max(data_dC13CH4_aerchemmip_gfdl,data_dC13CH4_aerchemmip_uk));
end
if(plot_ccmi); data_dC13CH4_ccmi = data_dC13CH4_ccmi-run2data.data_dC13CH4_ccmi; end
if(plot_osloctm3)
    data_dC13CH4_osloctm3 = reverse_diff*(data_dC13CH4_osloctm3-run2data.data_dC13CH4_osloctm3);
    if(INP_inc_dC13emis_unc>0)
        data_dC13CH4_osloctm3_high = reverse_diff*(run1data_high.data_dC13CH4_osloctm3-run2data_high.data_dC13CH4_osloctm3);
        data_dC13CH4_osloctm3_low  = reverse_diff*(run1data_low.data_dC13CH4_osloctm3-run2data_low.data_dC13CH4_osloctm3);
    end
end

% Normalize the dC13CH4 observation time series to give anomalies with reference to the common time period of Schaefer et al. and WMO/GAW
if(plot_dC13CH4_anomaly)
    if(idealized_OH<2)
        if(exist('save_runname_version','var') && save_runname_version==2)
            syear_obs = max(years_dC13CH4_obs(1),years_dC13CH4_obs2(1));
            eyear_obs = min(years_dC13CH4_obs(end),years_dC13CH4_obs2(end));
            data_dC13CH4_obs = data_dC13CH4_obs - mean(data_dC13CH4_obs(syear_obs-years_dC13CH4_obs(1)+1:eyear_obs-years_dC13CH4_obs(1)+1));
            data_dC13CH4_obs2_min = data_dC13CH4_obs2_min - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
            data_dC13CH4_obs2_max = data_dC13CH4_obs2_max - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
            data_dC13CH4_obs2= data_dC13CH4_obs2 - mean(data_dC13CH4_obs2(syear_obs-years_dC13CH4_obs2(1)+1:eyear_obs-years_dC13CH4_obs2(1)+1));
        else
            syear_obs = max(years_dC13CH4_obs(1),syear_WMO);
            eyear_obs = min(years_dC13CH4_obs(end),eyear_WMO);
            data_dC13CH4_obs = data_dC13CH4_obs - mean(data_dC13CH4_obs(syear_obs-years_dC13CH4_obs(1)+1:eyear_obs-years_dC13CH4_obs(1)+1));
            data_dC13CH4_obs2_min = data_dC13CH4_obs2_min - mean(data_dC13CH4_obs2(syear_obs-syear+1:eyear_obs-syear+1));
            data_dC13CH4_obs2_max = data_dC13CH4_obs2_max - mean(data_dC13CH4_obs2(syear_obs-syear+1:eyear_obs-syear+1));
            data_dC13CH4_obs2= data_dC13CH4_obs2 - mean(data_dC13CH4_obs2(syear_obs-syear+1:eyear_obs-syear+1));
        end
    end
end


% Plotting
% ********
% Colors
c_aerchemmip = [77,175,74]/255;
c_ccmi = [55,126,184]/255;
c_osloctm3 = [228,26,28]/255;
cat_colors = [217,95,2; 117,112,179; 231,41,138; 230,171,2]/255;


% Time evolution of d13C
p = figure; set(p,'PaperPositionMode','auto'); set(p,'Position',[0 0 1000 600]);
hold on;
lgd = legend('FontSize',8,'Location','South');%NorthWest');
if(plot_dC13CH4_anomaly) % plot zero line
    set(lgd,'AutoUpdate','off');
    plot(years,zeros(nyears,1),'-k');%,'LineWidth',1);
    set(lgd,'AutoUpdate','on');
end
if(idealized_OH<2)
    set(lgd,'AutoUpdate','off');
    if(exist('save_runname_version','var') && save_runname_version==2)
        patch([years_dC13CH4_obs2 fliplr(years_dC13CH4_obs2)],[data_dC13CH4_obs2_max fliplr(data_dC13CH4_obs2_min)],[0.25,0.25,0.25],'EdgeColor','none','FaceAlpha',0.333);
        set(lgd,'AutoUpdate','on');
        plot(years_dC13CH4_obs2,data_dC13CH4_obs2,'-','Color',[0.25,0.25,0.25],'LineWidth',2,'DisplayName','Obs. WMO/GAW');
    else
        patch([years(s_index_WMO:e_index_WMO) fliplr(years(s_index_WMO:e_index_WMO))],[data_dC13CH4_obs2_max(s_index_WMO:e_index_WMO)' fliplr(data_dC13CH4_obs2_min(s_index_WMO:e_index_WMO)')],[0.25,0.25,0.25],'EdgeColor','none','FaceAlpha',0.333);
        set(lgd,'AutoUpdate','on');
        plot(years,data_dC13CH4_obs2,'-','Color',[0.25,0.25,0.25],'LineWidth',2,'DisplayName','Obs. WMO/GAW');
    end
    plot(years_dC13CH4_obs,data_dC13CH4_obs,'-k','LineWidth',2,'DisplayName','Obs. Schaefer et al. (2016)');
end
if(plot_aerchemmip)
    e_index_aerchemmip = min(eyear,eyear_OH_aerchemmip)-syear+1;
    set(lgd,'AutoUpdate','off');
    patch([years(1:e_index_aerchemmip) fliplr(years(1:e_index_aerchemmip))],[data_dC13CH4_aerchemmip_max(1:e_index_aerchemmip)' fliplr(data_dC13CH4_aerchemmip_min(1:e_index_aerchemmip)')],c_aerchemmip,'EdgeColor','none','FaceAlpha',0.333);
    set(lgd,'AutoUpdate','on');
    plot(years(1:eyear_OH_aerchemmip-syear+1),data_dC13CH4_aerchemmip(1:eyear_OH_aerchemmip-syear+1),'Color',c_aerchemmip,'LineWidth',2,'DisplayName','AerChemMIP');
end
if(plot_ccmi)
    s_index_ccmi = max(syear,syear_OH_ccmi)-syear+1;
    e_index_ccmi = min(eyear,eyear_OH_ccmi)-syear+1;
    plot(years(1:eyear_OH_ccmi-syear+1),data_dC13CH4_ccmi(1:eyear_OH_ccmi-syear+1),'Color',c_ccmi,'LineWidth',2,'DisplayName','CCMI');
end
if(plot_osloctm3)
    if(INP_anthro_const2007>0 || INP_nat_emis>0) % plot solid line all the way for the emission perturbations (possibly except natural)
        if(INP_anthro_const2007>0)
            cat_label = sources_CH4(INP_anthro_const2007).name;
            cat_color = cat_colors(INP_anthro_const2007,:);
            eyear_emis = eyear;
        elseif(INP_nat_emis==1)
            cat_label = 'CLM GSWP3'; cat_color = [102,194,164]/255;
            eyear_emis = 2014;
        elseif(INP_nat_emis==2)
            cat_label = 'CLM CRUNCEP'; cat_color = [44,162,95]/255;
            eyear_emis = 2016;
        elseif(INP_nat_emis==3)
            cat_label = 'CLM CRUJRA'; cat_color = [0,109,44]/255;
            eyear_emis = 2019;
        elseif(INP_nat_emis==10)
            cat_label = 'LPJ\_wsl (MERRA2)'; cat_color = [106,61,154]/255;
            eyear_emis = 2020;
        elseif(INP_nat_emis==11)
            cat_label = 'LPJ\_wsl (CRU)'; cat_color = [31,120,180]/255;
            eyear_emis = 2020;
        end
        if(INP_inc_dC13emis_unc>0) % plot uncertainty as shading
            set(lgd,'AutoUpdate','off');
            patch([years(1:eyear_emis-syear+1) fliplr(years(1:eyear_emis-syear+1))],[data_dC13CH4_osloctm3_high(1:eyear_emis-syear+1)' fliplr(data_dC13CH4_osloctm3_low(1:eyear_emis-syear+1)')],cat_color,'EdgeColor','none','FaceAlpha',0.333);
            set(lgd,'AutoUpdate','on');
        end
        plot(years(1:eyear_emis-syear+1),data_dC13CH4_osloctm3(1:eyear_emis-syear+1),'Color',cat_color,'LineWidth',2,'DisplayName',cat_label);
    else
        plot(years(1:eyear_OH_osloctm3-syear+1),data_dC13CH4_osloctm3(1:eyear_OH_osloctm3-syear+1),'Color',c_osloctm3,'LineWidth',2,'DisplayName','OsloCTM3 (CEDS21+COVID)');
    end
end
set(lgd,'AutoUpdate','off');
if(plot_aerchemmip && eyear_OH_aerchemmip<eyear) % then do the latter years with constant OH, if relevant
    plot(years(eyear_OH_aerchemmip-syear+1:nyears),data_dC13CH4_aerchemmip(eyear_OH_aerchemmip-syear+1:nyears),':','Color',c_aerchemmip,'LineWidth',2);
end
if(plot_ccmi && eyear_OH_ccmi<eyear)
    plot(years(eyear_OH_ccmi-syear+1:nyears),data_dC13CH4_ccmi(eyear_OH_ccmi-syear+1:nyears),':','Color',c_ccmi,'LineWidth',2);
end
if(plot_osloctm3 && eyear_OH_osloctm3<eyear && INP_anthro_const2007==0 && INP_nat_emis==0)
    plot(years(eyear_OH_osloctm3-syear+1:nyears),data_dC13CH4_osloctm3(eyear_OH_osloctm3-syear+1:nyears),':','Color',c_osloctm3,'LineWidth',2);
end
hold off;
xlim([syear_plot eyear_plot]);
if(set_maxmin_d13C==1); ylim([minval_d13C maxval_d13C]); end
grid on;
if(plot_dC13CH4_anomaly); title('\Delta \delta^{13}C-CH_4'); else title('\delta^{13}C-CH_4'); end
ylabel(char(8240));
set(gcf,'Color','w');
savename = 'dC13CH4_combinedruns';
saveas(gcf,['./plots/',savename,'.fig']);
print(gcf,['./plots/',savename,'.png'],'-dpng','-r300');
