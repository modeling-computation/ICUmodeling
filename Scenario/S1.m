%% NPI Scenario. NPI level control
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

level_1  = xlsread('../Result/beta_level1.xlsx');
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

scenario_initial=xx(end,:);
%%
% Scenario 2021.11.01-2021.12.18 
step = 4;

scenario_day=total_day-scenario_start+1;
S  = zeros(scenario_day,age_num,step);
V  = zeros(scenario_day,age_num,step);
E  = zeros(scenario_day,age_num,step);
I  = zeros(scenario_day,age_num,step);
Hm = zeros(scenario_day,age_num,step);
Hs = zeros(scenario_day,age_num,step);
C  = zeros(scenario_day,age_num,step);
R  = zeros(scenario_day,age_num,step);
D  = zeros(scenario_day,age_num,step);

Bed = ones(scenario_day,step)*800;
IOR = zeros(scenario_day,step); 

dp  = inf; delay = 7;
for j = 1 : step

    total_x=scenario_initial;

    for i = 1:scenario_day

        phi=vaccine(i+scenario_start-1,:);P_s=severe_rate(i+scenario_start-1,:);

        IOR(i,j)=sum(total_x(seg(7)))/Bed(i,j);
        
        if (IOR(i,j)>=0.4)
            dp=min(i,dp);
        end

        if (i>=dp+delay)
            beta = npi(j,:);
        else
            beta = npi(1,:);
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

    S(:,:,j)  = xx(:,seg(1));
    V(:,:,j)  = xx(:,seg(2));
    E(:,:,j)  = xx(:,seg(3));
    I(:,:,j)  = xx(:,seg(4));
    Hm(:,:,j) = xx(:,seg(5));
    Hs(:,:,j) = xx(:,seg(6));
    C(:,:,j)  = xx(:,seg(7));
    R(:,:,j)  = xx(:,seg(8));
    D(:,:,j)  = xx(:,seg(9));
end

ap=dp+delay;
% Scenario 1-1 ICU patient
figure(1)
for k=1:step
    plot(1:scenario_day,sum(C(:,:,k),2),'LineWidth',2.5);hold on
end

xlim([1 scenario_day])
ylim([0 1500])
xline(dp,'--r', 'LineWidth',2);hold on
yline(800,'--r','LineWidth',2);hold on


% Scenario 1-1 IOR after 2021-11-01
figure(2)
for k=1:step
    plot(1:scenario_day,IOR(:,k),'LineWidth',2.5);hold on
end
xlim([1 scenario_day])
ylim([0 2])
h1=fill([ap ap scenario_day scenario_day], [0 2 2 0],'k');
set(h1,'facealpha',0.1)
xline(dp,'--r', 'LineWidth',2);hold on
yline(1,'--r','LineWidth',2);hold on
yline(0.8,'--b','LineWidth',2);hold on
yline(0.6,'--b','LineWidth',2);hold on


set(gca, 'fontsize',20)
ylabel('IOR')

%% Save result
save_path='..\Result\S1\';
g_name=["no","con1","con2","con3"];
g_name2=["g1","g2","g3","g4","g5","g6","g7","g8","g9","g10","g11","g12","g13","g14","g15","g16"];

% Case
S1_case=squeeze(sum(alpha_I*(1-rho)'.*I(:,:,:),2));
S1_case=array2table(S1_case,'VariableNames',g_name);
writetable(S1_case, strcat(save_path,'case.xlsx'));

% Case by age
for i=1:step
    S1_case=alpha_I*(1-rho)'.*I(:,:,i);
    S1_case=array2table(S1_case,'VariableNames',g_name2);
    writetable(S1_case, strcat(save_path,'case_',g_name(i),'.xlsx'));
end

% IOR
S1_IOR=array2table(IOR(:,:),'VariableNames',g_name);
writetable(S1_IOR, strcat(save_path,'IOR.xlsx'));

% new ICU
S1_new_icu=squeeze(sum((1-scenario_f).*scenario_alpha_s.*Hs,2));
S1_new_icu=array2table(S1_new_icu,'VariableNames',g_name);
writetable(S1_new_icu, strcat(save_path,'new_icu.xlsx'));

for i=1:step
    S1_new_icu=squeeze((1-scenario_f).*scenario_alpha_s.*Hs(:,:,i));
    S1_new_icu=array2table(S1_new_icu,'VariableNames',g_name2);
    writetable(S1_new_icu, strcat(save_path,'new_icu_',g_name(i),'.xlsx'));
end

% death
S1_death=squeeze(sum(D,2));
S1_death=array2table(S1_death,'VariableNames',g_name);
writetable(S1_death, strcat(save_path,'death.xlsx'));

for i=1:step
    S1_death=D(:,:,i);
    S1_death=array2table(S1_death,'VariableNames',g_name2);
    writetable(S1_death, strcat(save_path,'death_',g_name(i),'.xlsx'));
end

% icu
S1_icu=squeeze(sum(C(:,:,:),2));
S1_icu=array2table(S1_icu,'VariableNames',g_name);
writetable(S1_icu, strcat(save_path,'icu.xlsx'));

for i=1:4
    S1_icu=squeeze(C(:,:,i));
    S1_icu=array2table(S1_icu,'VariableNames',g_name2);
    writetable(S1_icu, strcat(save_path,'icu_',g_name(i),'.xlsx'));
end

% Bed
S1_bed=array2table(Bed,'VariableNames',g_name);
writetable(S1_bed, strcat(save_path,'bed.xlsx'));