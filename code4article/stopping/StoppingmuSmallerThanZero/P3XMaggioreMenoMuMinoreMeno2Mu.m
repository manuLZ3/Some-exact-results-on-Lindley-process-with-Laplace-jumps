%function p = P3MuMaggioreMenoXMinoreMeno2x(h,mu,s,x)
function p = P3XMaggioreMenoMuMinoreMeno2Mu(h,mu,s,x)
    %{  
                            OLD -WRONG
    alpha_10 = exp((mu-h)/s);
    
    alpha_20 = exp(mu/s)*alpha_10*(s/2+h);
    alpha_21 = -exp(mu/s)*alpha_10;
    beta_20 = exp(-mu/s)*alpha_10*s*0.5;
    
    %p = exp((-mu-x)/s)/(8*s^2)*(  alpha_20*s/2*( exp(2*(x+mu)/s)-1 ) + alpha_21*( s/2*(exp(2*(x+mu)/s)*(x+2*mu)-mu) - (s/2)^2*(exp(2*(x+mu)/s)-1)) + beta_20*(x+mu)  ) + exp((x+mu)/s)/(8*s^2)*( alpha_20*(h-x-mu) + alpha_21*((h+mu)^2/2)-(x+2*mu)^2/2) + s/2*beta_20*(exp(-2*(x+mu)/s)-exp(-2*h/s) )+ P2MuMaggioreMenoX(h,mu,s,0)*exp((-x-mu)/s)/2;
      
    alpha_320 = alpha_20*((h+mu)+(s/2)) + alpha_21 *((h+mu)^2/2-(s/2)^2) + beta_20*(-(s/2)*exp(-2*h/s));
    alpha_321 = -alpha_20 + alpha_21 *(s/2);
    alpha_322 = - 0.5*alpha_21;
    %beta_32 = 0;
    beta_321 = beta_20;
    beta_320 = beta_20*(-mu+s/2) + alpha_20*(-s/2) + alpha_21*(-mu*(s/2)+(s/2)^2);

    p = exp((x+mu)/s)/(8*s^2)* (alpha_320 + alpha_321*(x+2*mu)+ alpha_322*(x+2*mu)^2) + exp((-x-mu)/s)/(8*s^2)* (beta_320 + beta_321*(x+2*mu)) + P2MuMaggioreMenoX(h,mu,s,0)*exp((-x-mu)/s)/2;
    %}

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
    alpha_320 = (h+mu)*alpha_220 + s/2*(exp(2*mu/s)-exp(-2*h/s))*beta_220 + (h+mu)^2/2*alpha_221 + s/2*alpha_210 - s*exp(mu/s)*eta_21;
    alpha_320 = alpha_320*exp(mu/s)/(2*s);
    alpha_321 = -alpha_210;
    alpha_321 = alpha_321*exp(mu/s)/(2*s);
    beta_320  = -s/2*alpha_210-s*eta_21;
    beta_320  = beta_320*exp(-mu/s)/(2*s) + exp(-mu/s)/2*(alpha_210 + eta_21); % l'ultima è la prob di fermarsi in 2 partendo da 0
    eta_32 = eta_21;

    p = exp(x/s)*( alpha_320+alpha_321*(x+2*mu) ) + exp(-x/s)*beta_320 + eta_32;


end