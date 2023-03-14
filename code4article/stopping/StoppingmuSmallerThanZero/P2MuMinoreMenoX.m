function p = P2MuMinoreMenoX(h,mu,s,x)
    alpha_10 = exp((mu-h)/s)/2;
   
    % CORRETTA
    %p = h/(2*s)*alpha_10*exp((+x+mu)/s) + (1-1/2*exp((+x+mu)/s))*alpha_10;

    alpha_20 = (h-s)*exp((+mu)/s)*alpha_10/(2*s);
    eta_2 = alpha_10;

    p = alpha_20*exp((+x)/s) + eta_2;