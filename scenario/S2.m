%% Scenario 2 Comprehensive scenario
clear all
clc

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
scenario_start=249;
scenario_end=296;

beta_delay=7; 
bed_delay=7;

threshold=0.7;

re_bed=[0,0.05,0.1,0.15,0.2];

beta_step=4;
bed_step=5;

S=zeros(total_day,age_num,beta_step,bed_step);
V=zeros(total_day,age_num,beta_step,bed_step);
E=zeros(total_day,age_num,beta_step,bed_step);
I=zeros(total_day,age_num,beta_step,bed_step);
Hm=zeros(total_day,age_num,beta_step,bed_step);
R=zeros(total_day,age_num,beta_step,bed_step);
Hs=zeros(total_day,age_num,beta_step,bed_step);
C=zeros(total_day,age_num,beta_step,bed_step);
D=zeros(total_day,age_num,beta_step,bed_step);

Bed=ones(total_day,beta_step,bed_step)*800;
IOR=zeros(total_day,beta_step,bed_step); % ior

dp=inf; 



for k=1:bed_step

    for j = 1 : beta_step

    total_x=initial_x;
    beta_dp=inf;
    bed_dp=inf;

    for i = 1:total_day
        phi=vaccine(i,:);
        
        IOR(i,j,k)=sum(total_x(7*age_num+1:8*age_num))/Bed(i,j,k);

        if (IOR(i,j,k)>=threshold) && (i>=scenario_start)
            beta_dp=min(i,beta_dp);
            bed_dp=min(i,bed_dp);
        end
        
        if (i>=beta_dp+beta_delay)
            beta = fit_beta(i,:,j);
        else
            beta = fit_beta(i,:,1);
        end
         
        if (i==bed_dp+bed_delay-1)
            Bed(i+1:end,j,k)=Bed(i,j,k)*(1+re_bed(k));
            bed_dp=inf;
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


    S(:,:,j,k)=xx(:,1:1*age_num);
    V(:,:,j,k)=xx(:,1*age_num+1:2*age_num);
    E(:,:,j,k)=xx(:,2*age_num+1:3*age_num);
    I(:,:,j,k)=xx(:,3*age_num+1:4*age_num);
    Hm(:,:,j,k)=xx(:,4*age_num+1:5*age_num);
    R(:,:,j,k)=xx(:,5*age_num+1:6*age_num);
    Hs(:,:,j,k)=xx(:,6*age_num+1:7*age_num);
    C(:,:,j,k)=xx(:,7*age_num+1:8*age_num);
    D(:,:,j,k)=xx(:,8*age_num+1:9*age_num);
    end
end

beta_colors=["blue","green","red","yellow","magenta"];
bed_colors=["cyan","black","#4DBEEE","#A2142F","#EDB120"];

% Scenario 3 ICU and Bed
for k=1:bed_step
    for j=1:beta_step
        subplot(5,5,5*(j-1)+k);
        plot(scenario_start:scenario_end,sum(C(scenario_start:scenario_end,:,j,k),2),'Color',beta_colors(j),'LineWidth',2.5);hold on
        plot(scenario_start:scenario_end,Bed(scenario_start:scenario_end,j,k),'Color',bed_colors(k),'LineWidth',2.5);hold on
        xlim([scenario_start scenario_end])
        ylim([0 2000])
    end
end


%% Save result
save_path='..\Result\S2\';
g_name=string(1:20);

% Case
S2_case=squeeze(alpha*sum(I(scenario_start:scenario_end,:,:),2));
S2_case=array2table(S2_case,'VariableNames',g_name);
writetable(S2_case, strcat(save_path,'case.xlsx'));

% IOR
S2_IOR=array2table(IOR(scenario_start:scenario_end,:),'VariableNames',g_name);
writetable(S2_IOR, strcat(save_path,'IOR.xlsx'));

% new ICU
S2_new_icu=squeeze(sum(fit_pc(scenario_start:scenario_end,1:17).*Hs(scenario_start:scenario_end,:,:),2));
S2_new_icu=array2table(S2_new_icu,'VariableNames',g_name);
writetable(S2_new_icu, strcat(save_path,'new_icu.xlsx'));


% death
S2_death=squeeze(sum(D(scenario_start:scenario_end,1:17,:),2));
S2_death=array2table(S2_death,'VariableNames',g_name);
writetable(S2_death, strcat(save_path,'death.xlsx'));

% icu
S2_icu=squeeze(sum(C(scenario_start:scenario_end,:,:),2));
S2_icu=array2table(S2_icu,'VariableNames',g_name);
writetable(S2_icu, strcat(save_path,'icu.xlsx'));

% Bed
S2_bed=array2table(Bed(scenario_start:scenario_end,:),'VariableNames',g_name);
writetable(S2_bed, strcat(save_path,'bed.xlsx'));