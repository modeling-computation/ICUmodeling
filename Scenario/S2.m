%% Comprehensive Scenario 2-1 threshold 0.7, delay 7
clear all; clc;

total_day = 296; % 2021-02-26 ~ 2021-12-18
age_num   = 16; 

% States per age group: S V E I Hm Hs C R D
state    = 9; 
statenum = age_num*state;

% parameters
h       = 1;
num     = 1/h;
W_R     = 1/180; % waning recover
W_V     = 1/180; % waning vaccination
alpha_E = 1/3.5; % Progression exposed to infection
alpha_I = 1/6.8; % Progression infection to Hospitalized
r_UR    = 1/14;  % Unreported cases to Recovery
r_mR    = 1/14;  % Progression Mild to Recovery
r_CR    = 1/21;  % Progression Critical to Recovery
sigma   = 0.911; % Vaccine efficacy

% load Data
initial     = xlsread('../Data/initial.xlsx'); % intital
rho         = xlsread('../Data/unreported_rate.xlsx')/100; % Unreported rate
cm          = xlsread('../Data/contact_matrix.xlsx'); % contact matrix
population  = xlsread('../Data/population.xlsx'); % popultion by age group
NN          = population(1:end-1);

params  = struct('age_num',age_num,'h',h, 'W_R',W_R, 'W_V',W_V, 'alpha_E',alpha_E, 'alpha_I',alpha_I,...
                 'r_UR',r_UR,'r_mR',r_mR,'r_CR',r_CR,'sigma',sigma, 'cm',cm, 'rho',rho, 'NN',NN);

% confirmed cases
infection_num = xlsread('../Data/covid19_cases.xlsx');
infection_num = infection_num(1:total_day,2:age_num+1);

% death
death_num = xlsread('../Data/covid19_death.xlsx');
death_num = death_num(1:total_day,:);

% ICU
ICU = xlsread('../Data/current_icu.xlsx');
ICU = ICU(1:total_day,2:age_num+1);

% Daily vaccination
vaccine = xlsread('../Data/vaccine.xlsx'); 
vaccine = vaccine(1:total_day,:);

% severe_rate
severe_rate = xlsread('../Data/severe_rate_day.xlsx');
severe_rate = severe_rate(1:total_day,:);

% initial value : If it does not exist previous result, then 0
seg = @(k) ((k-1)*age_num+1) : (k*age_num);

initial_x         = zeros(1,statenum);                                                    
initial_x(seg(1)) = initial(seg(1));                                            % S
initial_x(seg(2)) = zeros(1,age_num);                                           % V
initial_x(seg(3)) = initial(seg(2));                                            % E
initial_x(seg(4)) = initial(seg(3));                                            % I
initial_x(seg(5)) = initial(seg(4)).*(1-severe_rate(1,:));                      % Hm
initial_x(seg(6)) = initial(seg(4)).*(severe_rate(1,:));                        % Hs
initial_x(seg(7)) = [0,0,0,0,0,0,0.5,0.5,0.5,0.5,7.5,7.5,19.5,19.5,27,61];      % C
initial_x(seg(8)) = initial(seg(5));                                            % R
initial_x(seg(9)) = initial(seg(6));                                            % D

for i=1:length(death_num)
    if i==1
        death=initial(seg(6));
    else
        death(i,:)=death(i-1,:)+death_num(i,seg(1)+1);
    end

end

fit_beta    = xlsread('../Result/fitting/beta_fitting.xlsx');
fit_alpha_s = xlsread('../Result/fitting/alpha_s_fitting.xlsx');
fit_f       = xlsread('../Result/fitting/f_fitting.xlsx');

result_f       = xlsread('../Result/fitting/result_f.xlsx');
result_alpha_s = xlsread('../Result/fitting/result_alpha_s.xlsx');

scenario_f = result_f(end,:);
scenario_alpha_s = result_alpha_s(end,:);

level_1=xlsread('../Result/beta_level1.xlsx');

npi      = zeros(4,age_num);
npi(1,:) = level_1*1.0;
npi(2,:) = level_1*0.7;
npi(3,:) = level_1*0.5;
npi(4,:) = level_1*0.3;
%%
% No scenario 2021-02-26 ~ 2021-11-30
scenario_start = 249; 

total_x=initial_x;
for i = 1:scenario_start-1
    phi=vaccine(i,:); P_s=severe_rate(i,:);

    beta_hat    = fit_beta(i,:);
    alpha_s_hat = fit_alpha_s(i,:);
    f_s_hat     = fit_f(i,:);
    f_c_hat     = 0.72*f_s_hat;

    temp_x=odeoperation_ode(@odef, statenum, [i, i+1], h, total_x, age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi, P_s, beta_hat, alpha_s_hat, f_s_hat, f_c_hat);

    total_x=temp_x(end,:);

    if i==1
        xx=temp_x(1:end-1,:);
    else
        xx=[xx;temp_x(1:end-1,:)];
    end

end

scenario_initial = xx(end,:);

%%
% Scenario 2021.11.01-2021.12.18 
scenario_day = total_day-scenario_start+1;

beta_step=4; bed_step=5;

S  = zeros(scenario_day,age_num,beta_step,bed_step);
V  = zeros(scenario_day,age_num,beta_step,bed_step);
E  = zeros(scenario_day,age_num,beta_step,bed_step);
I  = zeros(scenario_day,age_num,beta_step,bed_step);
Hm = zeros(scenario_day,age_num,beta_step,bed_step);
Hs = zeros(scenario_day,age_num,beta_step,bed_step);
C  = zeros(scenario_day,age_num,beta_step,bed_step);
R  = zeros(scenario_day,age_num,beta_step,bed_step);
D  = zeros(scenario_day,age_num,beta_step,bed_step);

Bed = ones(scenario_day,beta_step,bed_step)*800;
IOR = zeros(scenario_day,beta_step,bed_step); % ior

dp = inf; beta_delay = 7; bed_delay = 7;

re_bed = [0,0.05,0.1,0.15,0.2];

for k = 1 : bed_step
    for j = 1 : beta_step
        total_x = scenario_initial;
        beta_dp = inf;
        bed_dp  = inf;

        for i = 1:scenario_day
            phi=vaccine(i+scenario_start-1,:);P_s=severe_rate(i+scenario_start-1,:);
            
            IOR(i,j,k)=sum(total_x(seg(7)))/Bed(i,j,k);
    
            if (IOR(i,j,k)>=0.7)
                beta_dp = min(i,beta_dp);
                bed_dp  = min(i,bed_dp);
            end
            
            if (i>=beta_dp+beta_delay)
                beta = npi(j,:);
            else
                beta = npi(1,:);
            end
             
            if (i==bed_dp+bed_delay-1)
                Bed(i+1:end,j,k) = Bed(i,j,k)*(1+re_bed(k));
                bed_dp           = inf;
            end
    
            alpha_s = scenario_alpha_s;
            f_s     = scenario_f;
            f_c     = 0.72*f_s;
            
            temp_x=odeoperation_ode(@odef, statenum, [i,i+1], h, total_x, age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi, P_s, beta, alpha_s, f_s, f_c);

            total_x=temp_x(end,:);
    
            if i==1
                xx=temp_x(1:end-1,:);
            else
                xx=[xx;temp_x(1:end-1,:)];
            end
    
        end

        S(:,:,j,k)  = xx(:,seg(1));
        V(:,:,j,k)  = xx(:,seg(2));
        E(:,:,j,k)  = xx(:,seg(3));
        I(:,:,j,k)  = xx(:,seg(4));
        Hm(:,:,j,k) = xx(:,seg(5));
        Hs(:,:,j,k) = xx(:,seg(6));
        C(:,:,j,k)  = xx(:,seg(7));
        R(:,:,j,k)  = xx(:,seg(8));
        D(:,:,j,k)  = xx(:,seg(9));
    end
end

beta_colors=["blue","green","red","yellow","magenta"];
bed_colors=["cyan","black","#4DBEEE","#A2142F","#EDB120"];

% Scenario 3 ICU and Bed
for k=1:bed_step
    for j=1:beta_step
        subplot(5,5,5*(j-1)+k);
        plot(1:scenario_day,sum(C(:,:,j,k),2),'Color',beta_colors(j),'LineWidth',2.5);hold on
        plot(1:scenario_day,Bed(:,j,k),'Color',bed_colors(k),'LineWidth',2.5);hold on
        xlim([1 scenario_day])
        ylim([0 2000])
    end
end
%% Save result
save_path='..\Result\S2\';
g_name=string(1:20);

% Case
S_case=squeeze(sum(alpha_I*(1-rho)'.*I(:,:,:),2));
S_case=array2table(S_case,'VariableNames',g_name);
writetable(S_case, strcat(save_path,'case.xlsx'));

% IOR
S_IOR=array2table(IOR(:,:),'VariableNames',g_name);
writetable(S_IOR, strcat(save_path,'IOR.xlsx'));

% new ICU
S_new_icu=squeeze(sum((1-scenario_f).*scenario_alpha_s.*Hs(:,:,:),2));
S_new_icu=array2table(S_new_icu,'VariableNames',g_name);
writetable(S_new_icu, strcat(save_path,'new_icu.xlsx'));


% death
S_death=squeeze(sum(D(:,:,:),2));
S_death=array2table(S_death,'VariableNames',g_name);
writetable(S_death, strcat(save_path,'death.xlsx'));

% icu
S_icu=squeeze(sum(C(:,:,:),2));
S_icu=array2table(S_icu,'VariableNames',g_name);
writetable(S_icu, strcat(save_path,'icu.xlsx'));

% Bed
S_bed=array2table(Bed(:,:),'VariableNames',g_name);
writetable(S_bed, strcat(save_path,'bed.xlsx'));