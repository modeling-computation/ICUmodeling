function dx=odef(t,x,age_num, W_R, W_V, NN, beta, cm, eta, P_s, alpha, r_mR, f_s, r_CR, f_c, alpha_S, phi)

dS = zeros(age_num,1);
dV = zeros(age_num,1);
dE = zeros(age_num,1);
dI = zeros(age_num,1);

dHm = zeros(age_num,1);
dR = zeros(age_num,1);
dHs = zeros(age_num,1);
dC = zeros(age_num,1);
dD = zeros(age_num,1);


S=x(1:1*age_num);
V=x(1*age_num+1:2*age_num);
E=x(2*age_num+1:3*age_num);
I=x(3*age_num+1:4*age_num);
Hm=x(4*age_num+1:5*age_num);
R=x(5*age_num+1:6*age_num);
Hs=x(6*age_num+1:7*age_num);
C=x(7*age_num+1:8*age_num);
D=x(8*age_num+1:9*age_num);



for i=1:age_num

    dS(i)  = R(i)*W_R + V(i)*W_V - S(i)*beta(i)*sum(cm(i,:).*(I(1:age_num))./NN')-phi(i);
 
    dV(i)  = phi(i)-V(i)*W_V ;

    dE(i)  = S(i)*beta(i)*sum(cm(i,:).*(I(1:age_num))./NN') - E(i)*eta;

    dI(i)  = E(i)*eta - I(i)*alpha;

    dHm(i) = I(i)*(1-P_s(i))*alpha - Hm(i)*r_mR;

    dR(i)  = Hm(i)*r_mR + C(i)*(1-f_c(i))*r_CR - R(i)*W_R ;

    dHs(i) = I(i)*P_s(i)*alpha - Hs(i)*f_s(i) - Hs(i)*(1-f_s(i))*alpha_S(i);
    
    dC(i)  = Hs(i)*(1-f_s(i))*alpha_S(i) - C(i)*f_c(i) - C(i)*(1-f_c(i))*r_CR;

    dD(i)  = Hs(i)*f_s(i) + C(i)*f_c(i);

end

dx = [dS;dV;dE;dI;dHm;dR;dHs;dC;dD];
