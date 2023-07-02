% This function spins up the model to calculate initial OH and dC13_CH4 values for a given CH4 vmr and emission
% *************************************************************************************************************

function [ init_OH,init_dC13CH4 ] = calc_initvals( init_n_CH4,init_S_CH4,init_dC13_S_CH4,R_std,n_air,plot_spinup,...
                                                  k_CH4_OH,alpha_OH,...
                                                  k_CH4_Cl,alpha_Cl,const_Cl,...
                                                  k_CH4_soil,alpha_soil )


% Initial values for the spin-up
nyears = 300; % number of spin-up years
init_dC13CH4 = -50; % per mil

% Constants
dt = 3600*24*365; % seconds in a year


% Calculate the OH concentration that will give mass balance
init_OH = ((init_S_CH4/init_n_CH4)-k_CH4_Cl*const_Cl-k_CH4_soil)/k_CH4_OH; % molec/cm3

% Calculate loss of CH4 (molec/cm3/s)
loss_CH4_OH = k_CH4_OH*init_OH*init_n_CH4; % CH4+OH
loss_CH4_Cl = k_CH4_Cl*const_Cl*init_n_CH4; % CH4+Cl
loss_CH4_soil = k_CH4_soil*init_n_CH4; % CH4 loss to soil
loss_CH4 = loss_CH4_OH + loss_CH4_Cl + loss_CH4_soil;

% Calculate the reaction rates for each of 12CH4 and 13CH4
% Note that the k_12CH4 rates are assumed same as k_CH4 (loss will be scaled later in the code)
k_12CH4_OH  = k_CH4_OH; % reaction rate of OH with 12CH4 (cm3/molec/s)
k_13CH4_OH  = k_12CH4_OH*alpha_OH; % reaction rate of OH with 13CH4 (cm3/molec/s)
k_12CH4_Cl  = k_CH4_Cl; % reaction rate of Cl with 12CH4 (cm3/molec/s)
k_13CH4_Cl  = k_12CH4_Cl*alpha_Cl; % reaction rate of Cl with 13CH4 (cm3/molec/s)
k_12CH4_soil= k_CH4_soil; % soil sink rate for 12CH4 (1/s)
k_13CH4_soil= k_12CH4_soil*alpha_soil; % soil sink rate for 13CH4 (1/s)

% Initial fractions of C12 and C13
init_ratio_C13_C12 = R_std*(init_dC13CH4/1000+1); % see also https://ciaaw.org/carbon.htm and Table 5 in https://ciaaw.org/pubs/EXER-2000.pdf
init_C12_frac = 1/(1+init_ratio_C13_C12); % because C13_frac+C12_frac=1
init_C13_frac = 1-init_C12_frac;

% Initialize and loop through all years
data_dC13CH4 = nan(nyears,1); data_dC13CH4(1) = init_dC13CH4;
data_C13 = nan(nyears,1); data_C13(1) = init_C13_frac*init_n_CH4;
data_C12 = nan(nyears,1); data_C12(1) = init_C12_frac*init_n_CH4;
for y = 2:nyears

    % Calculate loss of C13 and C12 (molec/cm3/s)
    loss_C13_OH  = k_13CH4_OH*init_OH*data_C13(y-1);
    loss_C12_OH  = k_12CH4_OH*init_OH*data_C12(y-1);
    loss_C13_Cl  = k_13CH4_Cl*const_Cl*data_C13(y-1);
    loss_C12_Cl  = k_12CH4_Cl*const_Cl*data_C12(y-1);
    loss_C13_soil= k_13CH4_soil*data_C13(y-1);
    loss_C12_soil= k_12CH4_soil*data_C12(y-1);
    
    % Scale so that total loss of C13 and C12 to OH/Cl/soil equals loss of total CH4 to OH/Cl/soil
    loss_C13_OH = loss_C13_OH*loss_CH4_OH/(loss_C13_OH+loss_C12_OH);
    loss_C12_OH = loss_C12_OH*loss_CH4_OH/(loss_C13_OH+loss_C12_OH);
    loss_C13_Cl = loss_C13_Cl*loss_CH4_Cl/(loss_C13_Cl+loss_C12_Cl);
    loss_C12_Cl = loss_C12_Cl*loss_CH4_Cl/(loss_C13_Cl+loss_C12_Cl);
    loss_C13_soil = loss_C13_soil*loss_CH4_soil/(loss_C13_soil+loss_C12_soil);
    loss_C12_soil = loss_C12_soil*loss_CH4_soil/(loss_C13_soil+loss_C12_soil);
        
    % Sum the various loss terms for C13 and C12
    loss_C13 = loss_C13_OH + loss_C13_Cl + loss_C13_soil;
    loss_C12 = loss_C12_OH + loss_C12_Cl + loss_C12_soil;
    
    % Calculate a "net source term" by using S_CH4 but scaled by the fraction of C13 and C12
    % in CH4 emissions (delta13 values for CH4 emissions is calculated further up)
    ratio_C13_C12_S_CH4 = R_std*(init_dC13_S_CH4/1000+1);
    frac_C12_S_CH4 = 1/(1+ratio_C13_C12_S_CH4); % because frac_C13+frac_C12_frac=1
    frac_C13_S_CH4 = 1-frac_C12_S_CH4;
    S_C13 = init_S_CH4*frac_C13_S_CH4;
    S_C12 = init_S_CH4*frac_C12_S_CH4;

    % Calculate change in C13 and C12 using mass balance
    data_C13(y) = data_C13(y-1) + (S_C13 - loss_C13)*dt;
    data_C12(y) = data_C12(y-1) + (S_C12 - loss_C12)*dt;

    % Calculate the new delta13 isotopic ratio
    data_dC13CH4(y) = (data_C13(y)/data_C12(y)/R_std-1)*1000; % (per mil)
end

init_CH4 = init_n_CH4/(1e-9*n_air); % molec(CH4)/cm3 -> ppb

init_dC13CH4 = data_dC13CH4(end); % set new initial value

if(plot_spinup)
    figure; plot([1:nyears],data_dC13CH4); ylabel(char(8240));
    title(['Spin-up of \delta^{13}C-CH_4 for CH_4=',num2str(init_CH4,4),'  ppb, S_{CH4}=',num2str(init_S_CH4,6),' molec/cm3/s']);
end
