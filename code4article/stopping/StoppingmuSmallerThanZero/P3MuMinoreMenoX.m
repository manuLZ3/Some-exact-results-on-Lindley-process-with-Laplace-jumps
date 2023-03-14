function p = P3MuMinoreMenoX(h,mu,s,x)
    
    %% coefficienti passo 1
    alpha_10 = exp((mu-h)/s)/2;
    
    %% coefficienti passo 2 partendo dal primo bin
    alpha_210 = (h-s)*exp((+mu)/s)*alpha_10/(2*s);
    eta_210 = alpha_10;
    
    %% coefficienti passo 2 partendo dal secondo bin
    alpha_220 = exp(mu/s)*alpha_10*(s*0.5+h)/(2*s);
    alpha_221 = -exp(mu/s)*alpha_10/(2*s);
    beta_220 = exp(-mu/s)*alpha_10*s*0.5/(2*s);

    %% coefficienti passo 3 partendo dal primo bin
    eta_31 = alpha_210+eta_210;
    alpha_310 = ( -mu*alpha_210 + (1-exp(mu/s))*s*eta_210 + (h+mu)*alpha_220 + (h+mu)^2/2*alpha_221 + s/2*(exp(2*mu/s)-exp(-2*h/s))*beta_220 );
    alpha_310 = exp(mu/s)/(2*s)*alpha_310-exp(mu/s)/2*(alpha_210+eta_210);


    %p = exp((x+mu)/s)/(2*s)*( -mu*alpha_210 + (1-exp(mu/s))*s*eta_210 + (h+mu)*alpha_220 + (h+mu)^2/2*alpha_221 + s/2*(exp(2*mu/s)-exp(-2*h/s))*beta_220 )...
     %   + (1-exp((mu+x)/s)/2)*(alpha_210+eta_210);
    p = alpha_310*exp(x/s)+eta_31;