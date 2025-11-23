%% ========================================================================
%  SEIUHCRD-V age-structured fitting 
%  - Date: 2021-02-26 ~ 2021-12-18 (total_day = 296)
% ========================================================================
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
population  = xlsread('../Data/Population.xlsx'); % popultion by age group
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
%%
starttime = [0, 69, 114, 136, 150, 192, 211, 234, 248];
endtime   = [starttime(2:end), total_day-1]; 
nPeriods  = numel(starttime);

result_beta     = zeros(nPeriods, age_num);
result_alpha_s  = zeros(nPeriods, age_num);
result_f        = zeros(nPeriods, age_num);

beta_initial = [0.02*ones(1,age_num),0.001*ones(1,age_num),0.001*ones(1,age_num)]; 

total_x = initial_x;
xx      = [];

container=cell(1, 5);
container{5} = params;

for i=1:nPeriods
    t0 = starttime(i);t1 = endtime(i);
    time_interval = [t0, t1];

    phi = vaccine(t0+1:t1+1,:); P_s = severe_rate(t0+1:t1+1,:);

    container{1} = time_interval;
    container{2} = total_x;
    container{3} = phi;
    container{4} = P_s;

    % Fit targets: [cases_by_age, ICU_by_age, death_cum_by_age]
    fit_cases = infection_num(t0+1:t1+1,:);
    fit_icu   = ICU(t0+1:t1+1,:);
    fit_death = death(t0+1:t1+1,:);
    fitdata   = [fit_cases, fit_icu, fit_death];
    
    % Bounds
    lb = [zeros(1,age_num), zeros(1,age_num), zeros(1,age_num)];
    ub = [100*ones(1,age_num), ones(1,age_num), ones(1,age_num)];
    
    % Fitting
    [x,resnorm,residual,exitflag,output,lambda,jacobian]=...
        lsqcurvefit(@odeoperation, beta_initial, container, fitdata, lb, ub);

    % Estimated parameters
    beta_hat    = x(seg(1));
    alpha_s_hat = x(seg(2));
    f_s_hat     = x(seg(3)); 
    f_c_hat     = 0.72*f_s_hat;

    result_beta(i,:)     = beta_hat;
    result_alpha_s(i,:)  = alpha_s_hat;
    result_f(i,:)        = f_s_hat;
    
    totalnum=[t0,t1];
    temp_x=odeoperation_ode(@odef, statenum, totalnum, h, total_x, age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi, P_s, beta_hat, alpha_s_hat, f_s_hat, f_c_hat);
    
    total_x=temp_x(end,:);

    if i==1
        xx=temp_x(1:end-1,:);
    else
        xx=[xx;temp_x(1:end-1,:)];
    end 
    
end

xx = [xx;temp_x(end,:)];

S  = xx(:,seg(1));
V  = xx(:,seg(2));
E  = xx(:,seg(3));
I  = xx(:,seg(4));
Hm = xx(:,seg(5));
Hs = xx(:,seg(6));
C  = xx(:,seg(7));
R  = xx(:,seg(8));
D  = xx(:,seg(9));

beta_fitting     = zeros(total_day, age_num);
alpha_s_fitting  = zeros(total_day, age_num);
f_fitting        = zeros(total_day, age_num);

% inclusive starts; fill until next start-1
starttime2 = [starttime, total_day-1];

for i = 1:(numel(starttime2)-1)
    t0 = starttime2(i)   + 1; t1 = starttime2(i+1) + 1;   
    beta_fitting(t0:t1, :)   = repmat(result_beta(i, :), t1 - t0 + 1, 1);
    alpha_s_fitting(t0:t1,:) = repmat(result_alpha_s(i,   :), t1 - t0 + 1, 1);
    f_fitting(t0:t1,    :)   = repmat(result_f(i,    :), t1 - t0 + 1, 1);
end
%% Plotting
plot_time=[0,total_day-1];
time=plot_time(1):plot_time(2);
time_plot=plot_time(1)+1:plot_time(2)+1;


figure(1)
plot(time,sum(alpha_I*(1-rho)'.*I(time_plot,:),2),'LineWidth',3);hold on
plot(time,sum(infection_num(time_plot,:),2),'*r');hold on
for j=2:length(starttime)
    xline(starttime(j),'--r');hold on
end
set(gca, 'FontSize', 20)
set(gca, 'FontWeight', 'bold')
xticks([0 total_day-1])
xticklabels({'2021-02-26','2021-12-18'})
legend({'Fitting','Case (Data)'})
axis tight

figure(2)
for i=1:age_num
    subplot(6,3,i);
    plot(time,alpha_I*(1-rho(i))*I(time_plot,i),'LineWidth',3);hold on
    plot(time,infection_num(time_plot,i),'*r','MarkerSize',1);hold on
    for j=2:length(starttime)
        xline(starttime(j),'--r');hold on
    end
end


figure(3)
plot(time,sum(C(time_plot,:),2),'Color',[1 0.8 0],'LineWidth',3);hold on
plot(time,sum(ICU(time_plot,:),2),'*r');hold on
for j=2:length(starttime)
    xline(starttime(j),'--r');hold on
end
set(gca, 'FontSize', 20)
set(gca, 'FontWeight', 'bold')
xticks([0 total_day-1])
xticklabels({'2021-02-26','2021-12-18'})
legend({'Fitting','ICU (Data)'})
axis tight


figure(4)
for i=1:age_num
    subplot(6,3,i);
    plot(time,C(time_plot,i),'LineWidth',3);hold on
    plot(time,ICU(time_plot,i),'*r','MarkerSize',1);hold on
    for j=2:length(starttime)
        xline(starttime(j),'--r');hold on
    end
end

figure(5)
plot(time,sum(D(time_plot,:),2),'Color','k','LineWidth',3);hold on
plot(time,sum(death(time_plot,:),2),'*r');hold on
for j=2:length(starttime)
    xline(starttime(j),'--r');hold on
end
set(gca, 'FontSize', 20)
set(gca, 'FontWeight', 'bold')
xticks([0 total_day-1])
xticklabels({'2021-02-26','2021-12-18'})
legend({'Fitting','Death (Data)'})
axis tight

beds=xlsread('..\Data\beds.xlsx'); 

figure(6)
for i=1:age_num
    subplot(6,3,i);
    plot(time,D(time_plot,i),'LineWidth',3);hold on
    plot(time,death(time_plot,i),'*r','MarkerSize',1);hold on
    for j=2:length(starttime)
        xline(starttime(j),'--r');hold on
    end

end

%%
cases_mae  = zeros(1, age_num);
cases_rmse = zeros(1, age_num);
icu_mae  = zeros(1, age_num);
icu_rmse = zeros(1, age_num);
death_mae  = zeros(1, age_num);
death_rmse = zeros(1, age_num);


cases_error  = infection_num(:,:)-alpha_I*(1-rho)'.*I(:,:);
icu_error   = ICU(:,:)-C(:,:);
death_error = death(:,:)-D(:,:);

cases_mae(1,:) = mean(abs(cases_error));
cases_rmse(1,:) = sqrt(mean(abs(cases_error).^2));

icu_mae(1,:) = mean(abs(icu_error));
icu_rmse(1,:) = sqrt(mean(abs(icu_error).^2));

death_mae(1,:) = mean(abs(death_error));
death_rmse(1,:) = sqrt(mean(abs(death_error).^2));


save_path='.\Result\';

writematrix(cases_mae,strcat(save_path,'cases_mae.xlsx'));
writematrix(cases_rmse,strcat(save_path,'cases_rmse.xlsx'));

writematrix(icu_mae,strcat(save_path,'icu_mae.xlsx'));
writematrix(icu_rmse,strcat(save_path,'icu_rmse.xlsx'));

writematrix(death_mae,strcat(save_path,'death_mae.xlsx'));
writematrix(death_rmse,strcat(save_path,'death_rmse.xlsx'));

%% Save result

save_path='.\Result\';

writematrix(beta_fitting,strcat(save_path,'beta_fitting.xlsx'));
writematrix(alpha_s_fitting,strcat(save_path,'alpha_s_fitting.xlsx'));
writematrix(f_fitting,strcat(save_path,'f_fitting.xlsx'));

writematrix(result_beta,strcat(save_path,'result_beta.xlsx'));
writematrix(result_alpha_s,strcat(save_path,'result_alpha_s.xlsx'));
writematrix(result_f,strcat(save_path,'result_f.xlsx'));

% fitting - cases by age
g_name=["g1","g2","g3","g4","g5","g6","g7","g8","g9","g10","g11","g12","g13","g14","g15","g16"];
fitcase=[sum(infection_num(time_plot,:),2),sum(alpha_I*(1-rho)'.*I(time_plot,:),2),alpha_I*(1-rho)'.*I(time_plot,:)];
fitcase=array2table(fitcase,'VariableNames',['data','case',g_name]);
writetable(fitcase, strcat(save_path,'fit_case.xlsx'));

% fitting - death by age
fitdeath=[sum(death(time_plot,:),2), sum(D(time_plot,:),2), D(time_plot,:)];
fitdeath=array2table(fitdeath,'VariableNames',['data','death',g_name]);
writetable(fitdeath, strcat(save_path,'fit_death.xlsx'));

% fitting - ICU by age
fitdeath=[sum(ICU(time_plot,:),2), sum(C(time_plot,:),2),C(time_plot,:)];
fitdeath=array2table(fitdeath,'VariableNames',['data','icu',g_name]);
writetable(fitdeath, strcat(save_path,'fit_icu.xlsx'));


% Compartment
total_compartment=array2table(xx);
writetable(total_compartment, strcat(save_path,'total_compartment.xlsx'));

