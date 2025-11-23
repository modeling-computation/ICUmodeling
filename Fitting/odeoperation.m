function C1_t = odeoperation(est,container)

time_interval = container{1};
x             = container{2};
phi           = container{3};
P_s           = container{4};
params        = container{5};

% defined parameters
age_num       = params.age_num;
NN            = params.NN;
cm            = params.cm;
h             = params.h;
W_R           = params.W_R; % waning recover
W_V           = params.W_V; % waning vaccination
alpha_E       = params.alpha_E; % Progression exposed to infection
alpha_I       = params.alpha_I; % Progression infection to Hospitalized
r_UR          = params.r_UR;  % Unreported cases to Recovery
r_mR          = params.r_mR;  % Progression Mild to Recovery
r_CR          = params.r_CR;  % Progression Critical to Recovery
rho           = params.rho;   % unreported rate
sigma         = params.sigma; % vaccine efficacy

% estimated parameters
beta          = est(1:age_num);
alpha_s       = est(age_num+1:2*age_num);
f_s           = est(2*age_num+1:3*age_num);
f_c           = 0.72*f_s;

% time interval
data_time     = time_interval(1):h:time_interval(2);

for i = 1:length(data_time)-1

   k1 = odef((i-1)*h,x(i,:), age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi(i,:), P_s(i,:), beta, alpha_s, f_s, f_c);
   k2 = odef((i-1)*h+h/2,x(i,:)+h*k1'/2, age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi(i,:), P_s(i,:), beta, alpha_s, f_s, f_c);
   k3 = odef((i-1)*h+h/2,x(i,:)+h*k2'/2, age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi(i,:), P_s(i,:), beta, alpha_s, f_s, f_c);
   k4 = odef(i*h,x(i,:)+h*k3', age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi(i,:), P_s(i,:), beta, alpha_s, f_s, f_c);
   x(i+1,:) = x(i,:) + h*(k1' + 2*k2' + 2*k3' + k4')/6;

end

C1_t=[alpha_I*(1-rho)'.*x(:,3*age_num+1:4*age_num), x(:,6*age_num+1:7*age_num), x(:,8*age_num+1:9*age_num)];