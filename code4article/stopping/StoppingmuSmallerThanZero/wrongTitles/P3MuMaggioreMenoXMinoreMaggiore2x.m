function p = P3MuMaggioreMenoXMaggiore2x(h,mu,s,x)
   
    % coefficienti primo passo
    alpha_10 = exp((mu-h)/s)/2;
    
    % coefficienti secondo passo, primo bin
    alpha_210 = (h-s)*exp((+mu)/s)*alpha_10/(2*s);
    eta_21 = alpha_10;

    % coefficienti secondo passo, secondo bin
    alpha_220 = exp(mu/s)*alpha_10*(s*0.5+h)/(2*s);
    alpha_221 = -exp(mu/s)*alpha_10/(2*s);
    beta_220 = exp(-mu/s)*alpha_10*s*0.5/(2*s);
    
    % coefficienti terzo passo, secondo bin
    alpha_330 = s/2*alpha_220 + (h+mu)*alpha_220 -(s/2)^2*alpha_221 + (h+mu)^2/2*alpha_221 - s/2*exp(-2*h/s)*beta_220;
    alpha_330 = alpha_330 *exp(mu/s)/(2*s);
    alpha_331 = s/2*alpha_221 - alpha_220;
    alpha_331 = alpha_331 *exp(mu/s)/(2*s);
    alpha_332 = -alpha_221/2;
    alpha_332 = alpha_332 *exp(mu/s)/(2*s);
    beta_320  = s/2*(exp(-2*mu/s)-1)*alpha_210+s*(exp(-mu/s)-1)*eta_21 - s/2*exp(-2*mu/s)*alpha_220 +(s/2)^2*exp(-2*mu/s)*alpha_221 + s/2*beta_220;
    beta_320  = beta_320*exp(-mu/s)/(2*s) + exp(-mu/s)/2*(alpha_210 + eta_21);  % l'ultima è la prob di fermarsi in 2 partendo da 0
    beta_321 = beta_220;
    beta_321 = beta_321*exp(-mu/s)/(2*s);

    p = exp(x/s)*( alpha_330+alpha_331*(x+2*mu)+alpha_332*(x+2*mu)^2 ) + exp(-x/s)*(beta_320 + beta_321*(x+2*mu));


end