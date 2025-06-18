function x=odeoperation_ode(FN,statenum,totalnum,h,initialx, age_num, W_R, W_V, NN, beta, cm, eta, P_s, alpha, r_mR, f_s, r_CR, f_c, mu,phi)

odefuctionname=FN;
total_iter=(totalnum(2)-totalnum(1))/h;
x=zeros(total_iter+1,statenum);
x(1,:)=initialx;

for i = 1:total_iter
   k1 = odefuctionname((i-1)*h,x(i,:), age_num, W_R, W_V, NN, beta, cm, eta, P_s(i,:), alpha, r_mR, f_s, r_CR, f_c, mu, phi(i,:));
   k2 = odefuctionname((i-1)*h+h/2,x(i,:)+h*k1'/2, age_num, W_R, W_V, NN, beta, cm, eta, P_s(i,:), alpha, r_mR, f_s, r_CR, f_c, mu, phi(i,:));
   k3 = odefuctionname((i-1)*h+h/2,x(i,:)+h*k2'/2, age_num, W_R, W_V, NN, beta, cm, eta, P_s(i,:), alpha, r_mR, f_s, r_CR, f_c, mu, phi(i,:));
   k4 = odefuctionname(i*h,x(i,:)+h*k3', age_num, W_R, W_V, NN, beta, cm, eta, P_s(i,:), alpha, r_mR, f_s, r_CR, f_c, mu, phi(i,:));
   x(i+1,:) = x(i,:) + h*(k1' + 2*k2' + 2*k3' + k4')/6;
end

