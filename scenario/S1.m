%% Scenario 1. NPI scenario 
clear all
clc

%% Setting
total_day=296;
age_num=17;

state=9;
statenum=age_num*state;

h = 1; 
num=1/h;

W_R=1/180; % waning recover
W_V=1/180; % waning vaccination
eta=1/3.5; % Progression exposed to infection
alpha=1/6.8; % Progression infection to Hospitalized
r_mR=1/14; % Progression Mild to Recovery
r_CR=1/21; % Progression Critical to Recovery

level_1=xlsread('..\parameters\beta_level1.xlsx');
level_2=xlsread('..\parameters\beta_level2.xlsx');
level_3=xlsread('..\parameters\beta_level3.xlsx');
level_4=xlsread('..\parameters\beta_level4.xlsx');

fit_beta=zeros(total_day,age_num,4);
fit_beta(:,:,1)=level_1;
fit_beta(:,:,2)=level_2;
fit_beta(:,:,3)=level_3;
fit_beta(:,:,4)=level_4;

fit_pc = xlsread('..\parameters\icu_admission_rate.xlsx'); % ICU admission rate
fit_f = xlsread('..\parameters\severe_case_fatality_rate.xlsx'); % Severe case-fatality rate

% intital
initial=xlsread('..\parameters\initial.xlsx'); 

% contact matrix
cm=xlsread('..\parameters\contact_matrix_all.xlsx');

% popultion by age group
population=xlsread('..\data\Population.xlsx');
NN=population(1:end-1);

% infection case
infection_num=xlsread('..\data\covid19_cases.xlsx');
infection_num=infection_num(1:total_day,:);

% death
death_num=xlsread('..\data\covid19_death.xlsx');
death_num=death_num(1:total_day,:);

for i=1:length(death_num)
    if i==1
        death=initial(5*age_num+1:6*age_num);
    else
        death(i,:)=death(i-1,:)+death_num(i,2:18);
    end
end

% vaccine effect data
vaccine1=xlsread('..\data\vaccine1.xlsx'); 
vaccine2=xlsread('..\data\vaccine2.xlsx'); 
vaccine3=xlsread('..\data\vaccine3.xlsx'); 

vaccine1=vaccine1(1:total_day,:);
vaccine2=vaccine2(1:total_day,:);
vaccine3=vaccine3(1:total_day,:);
vaccine=vaccine2;

%
severe_rate=xlsread('..\parameters\severe_rate_day.xlsx');

% initial value : If it does not exist previous result, then 0
initial_x=zeros(1,statenum);                                                    
initial_x(1:1*age_num)=initial(1:1*age_num)-[0,0,0,0,0,0,0,1,0,1,8,7,20,19,27,27,34]; % S
initial_x(1*age_num+1:2*age_num)=zeros(1,age_num); %V
initial_x(2*age_num+1:3*age_num)=initial(1*age_num+1:2*age_num); %E
initial_x(3*age_num+1:4*age_num)=initial(2*age_num+1:3*age_num); %I
initial_x(4*age_num+1:5*age_num)=initial(3*age_num+1:4*age_num).*(1-severe_rate(1,:)); %Hm
initial_x(5*age_num+1:6*age_num)=initial(4*age_num+1:5*age_num); %R
initial_x(6*age_num+1:7*age_num)=initial(3*age_num+1:4*age_num).*(severe_rate(1,:)); %Hs
initial_x(7*age_num+1:8*age_num)=[0,0,0,0,0,0,0,1,0,1,8,7,20,19,27,27,34]; %C
initial_x(8*age_num+1:9*age_num)=initial(5*age_num+1:6*age_num); %D

% ICU
ICU=xlsread('..\data\current_icu.xlsx');

%% Scenario
% Scenario 2021.11.01-2021.12.18 
step=4;

delay=7;
threshold=0.4;

scenario_start=249;
scenario_end=296;

S=zeros(total_day,age_num,step);
V=zeros(total_day,age_num,step);
E=zeros(total_day,age_num,step);
I=zeros(total_day,age_num,step);
Hm=zeros(total_day,age_num,step);
R=zeros(total_day,age_num,step);
Hs=zeros(total_day,age_num,step);
C=zeros(total_day,age_num,step);
D=zeros(total_day,age_num,step);

Bed=ones(total_day,step)*800;
IOR=zeros(total_day,step); 

dp=inf; 

for j = 1 : step

    total_x=initial_x;

    for i = 1:total_day
        phi=vaccine(i,:);
        
        IOR(i,j)=sum(total_x(7*age_num+1:8*age_num))/Bed(i,j);
        
        if (IOR(i,j)>=threshold) && (i>=scenario_start)
            dp=min(i,dp);
        end

        if (i>=dp+delay)
            beta = fit_beta(i,:,j);
        else
            beta = fit_beta(i,:,1);
        end

        mu=fit_pc(i,:);

        f_s=fit_f(i,:);
        f_c=0.72*f_s;

        P_s=severe_rate(i,:);
        
        temp_x=odeoperation_ode(@odef,statenum,[i,i+1],h,total_x, age_num, W_R, W_V, NN, beta, cm, eta, P_s, alpha, r_mR, f_s, r_CR, f_c, mu,phi);
        
        total_x=temp_x(end,:);

        if i==1
            xx=temp_x(1:end-1,:);
        else
            xx=[xx;temp_x(1:end-1,:)];
        end

    end

    S(:,:,j)=xx(:,1:1*age_num);
    V(:,:,j)=xx(:,1*age_num+1:2*age_num);
    E(:,:,j)=xx(:,2*age_num+1:3*age_num);
    I(:,:,j)=xx(:,3*age_num+1:4*age_num);
    Hm(:,:,j)=xx(:,4*age_num+1:5*age_num);
    R(:,:,j)=xx(:,5*age_num+1:6*age_num);
    Hs(:,:,j)=xx(:,6*age_num+1:7*age_num);
    C(:,:,j)=xx(:,7*age_num+1:8*age_num);
    D(:,:,j)=xx(:,8*age_num+1:9*age_num);

end

ap=dp+delay;
% Scenario 1 ICU patient total period
figure(1)
for k=1:step
    plot(1:total_day,sum(C(1:total_day,:,k),2),'LineWidth',2.5);hold on
end

xlim([1 total_day])
ylim([0 2000])
h1=fill([ap ap scenario_end scenario_end], [0 2000 2000 0],'k');
set(h1,'facealpha',0.1)

xline(dp,'--r', 'LineWidth',2);hold on
yline(800,'--b','LineWidth',2);hold on

legend('NPI Level 1','NPI Level 2','NPI Level 3','NPI Level 4')

% Scenario 1 ICU patient
figure(2)
for k=1:step
    plot(scenario_start:scenario_end,sum(C(scenario_start:scenario_end,:,k),2),'LineWidth',2.5);hold on
end

xlim([scenario_start scenario_end])
ylim([0 1500])

xline(dp,'--r', 'LineWidth',2);hold on
yline(800,'--r','LineWidth',2);hold on

legend('NPI Level 1','NPI Level 2','NPI Level 3','NPI Level 4')


% Scenario 1 IOR after 2021-11-01
figure(3)
for k=1:step
    plot(scenario_start:scenario_end,IOR(scenario_start:scenario_end,k),'LineWidth',2.5);hold on
end
xlim([scenario_start scenario_end])
ylim([0 2])
h1=fill([ap ap scenario_end scenario_end], [0 2 2 0],'k');
set(h1,'facealpha',0.1)
xline(dp,'--r', 'LineWidth',2);hold on
yline(1,'--r','LineWidth',2);hold on
yline(0.8,'--b','LineWidth',2);hold on
yline(0.6,'--b','LineWidth',2);hold on

legend('NPI Level 1','NPI Level 2','NPI Level 3','NPI Level 4')

set(gca, 'fontsize',20)
ylabel('IOR')



%% Save result
save_path='..\result\S1\';
g_name=["no","con1","con2","con3"];
g_name2=["g1","g2","g3","g4","g5","g6","g7","g8","g9","g10","g11","g12","g13","g14","g15","g16","g17"];

% Case
S1_case=squeeze(alpha*sum(I(scenario_start:scenario_end,:,:),2));
S1_case=array2table(S1_case,'VariableNames',g_name);
writetable(S1_case, strcat(save_path,'case.xlsx'));

% Case by age
for i=1:4
    S1_case=alpha*I(scenario_start:scenario_end,:,i);
    S1_case=array2table(S1_case,'VariableNames',g_name2);
    writetable(S1_case, strcat(save_path,'case_',g_name(i),'.xlsx'));
end

% IOR
S1_IOR=array2table(IOR(scenario_start:scenario_end,:),'VariableNames',g_name);
writetable(S1_IOR, strcat(save_path,'IOR.xlsx'));

% new ICU
S1_new_icu=squeeze(sum(fit_pc(scenario_start:scenario_end,1:17).*Hs(scenario_start:scenario_end,:,:),2));
S1_new_icu=array2table(S1_new_icu,'VariableNames',g_name);
writetable(S1_new_icu, strcat(save_path,'new_icu.xlsx'));

for i=1:4
    S1_new_icu=squeeze(fit_pc(scenario_start:scenario_end,1:17).*Hs(scenario_start:scenario_end,:,i));
    S1_new_icu=array2table(S1_new_icu,'VariableNames',g_name2);
    writetable(S1_new_icu, strcat(save_path,'new_icu_',g_name(i),'.xlsx'));
end

% death
S1_death=squeeze(sum(D(scenario_start:scenario_end,1:17,:),2));
S1_death=array2table(S1_death,'VariableNames',g_name);
writetable(S1_death, strcat(save_path,'death.xlsx'));

for i=1:4
    S1_death=D(scenario_start:scenario_end,1:17,i);
    S1_death=array2table(S1_death,'VariableNames',g_name2);
    writetable(S1_death, strcat(save_path,'death_',g_name(i),'.xlsx'));
end


% icu
S1_icu=squeeze(sum(C(scenario_start:scenario_end,:,:),2));
S1_icu=array2table(S1_icu,'VariableNames',g_name);
writetable(S1_icu, strcat(save_path,'icu.xlsx'));

for i=1:4
    S1_icu=squeeze(C(scenario_start:scenario_end,:,i));
    S1_icu=array2table(S1_icu,'VariableNames',g_name2);
    writetable(S1_icu, strcat(save_path,'icu_',g_name(i),'.xlsx'));
end

% Bed
S1_bed=array2table(Bed(scenario_start:scenario_end,:),'VariableNames',g_name);
writetable(S1_bed, strcat(save_path,'bed.xlsx'));