function dx=odef(t,x,age_num, W_R, W_V, NN, cm, alpha_E, alpha_I, rho, r_UR, r_mR, r_CR, sigma, phi, P_s, beta, alpha_s, f_s, f_c)
                        
dS  = zeros(age_num,1);
dV  = zeros(age_num,1);
dE  = zeros(age_num,1);
dI  = zeros(age_num,1);
dHm = zeros(age_num,1);
dHs = zeros(age_num,1);
dC  = zeros(age_num,1);
dR  = zeros(age_num,1);
dD  = zeros(age_num,1);

S  = x(1          :1*age_num);
V  = x(1*age_num+1:2*age_num);
E  = x(2*age_num+1:3*age_num);
I  = x(3*age_num+1:4*age_num);
Hm = x(4*age_num+1:5*age_num);
Hs = x(5*age_num+1:6*age_num);
C  = x(6*age_num+1:7*age_num);
R  = x(7*age_num+1:8*age_num);
D  = x(8*age_num+1:9*age_num);


for i=1:age_num

    dS(i)  = - beta(i)*sum(cm(i,:).*I./NN')*S(i) + W_R*R(i) + W_V*V(i) -phi(i);
    dV(i)  = phi(i) - W_V*V(i) - (1-sigma)*beta(i)*sum(cm(i,:).*I./NN')*V(i);
    dE(i)  = beta(i)*sum(cm(i,:).*I./NN')*S(i) + (1-sigma)*beta(i)*sum(cm(i,:).*I./NN')*V(i) - alpha_E*E(i);
    dI(i)  = alpha_E*E(i) - (1-rho(i))*alpha_I*I(i) - rho(i)*r_UR*I(i);
    dHm(i) = (1-rho(i))*(1-P_s(i))*alpha_I*I(i) - r_mR*Hm(i);
    dHs(i) = (1-rho(i))*P_s(i)*alpha_I*I(i) - f_s(i)*Hs(i) - (1-f_s(i))*alpha_s(i)*Hs(i);
    dC(i)  = (1-f_s(i))*alpha_s(i)*Hs(i) - f_c(i)*C(i) - (1-f_c(i))*r_CR*C(i);
    dR(i)  = r_mR*Hm(i) + (1-f_c(i))*r_CR*C(i) + rho(i)*r_UR*I(i) - W_R*R(i);
    dD(i)  = f_s(i)*Hs(i) + f_c(i)*C(i);

end

dx = [dS;dV;dE;dI;dHm;dHs;dC;dR;dD];
