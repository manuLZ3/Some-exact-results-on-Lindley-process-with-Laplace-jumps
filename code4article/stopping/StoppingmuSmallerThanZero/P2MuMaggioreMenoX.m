function p = P2MuMaggioreMenoX(h,mu,s,x)
    alpha_10 = exp((mu-h)/s)/2;
    
    %% formula (4.48b)
    %p = alpha_10*(exp((x+mu)/s)-exp(-(x+mu)/s))/8 + alpha_10*exp((x+mu)/s)/(4*s)*(h-x-mu) +alpha_10*exp((-x-mu)/s)/4;

    %% formula (4.48c)
    %p = exp(x/s)/(4*s)*exp(mu/s)*alpha_10*(s/2+h-(x+mu)) + exp(-x/s)/(4*s)*exp(-mu/s)*alpha_10*s*0.5;

    
    alpha_220 = exp(mu/s)*alpha_10*(s*0.5+h)/(2*s);
    alpha_221 = -exp(mu/s)*alpha_10/(2*s);
    beta_220 = exp(-mu/s)*alpha_10*s*0.5/(2*s);

    p = exp(x/s)*(alpha_220+alpha_221*(x+mu)) + exp(-x/s)*beta_220;
end