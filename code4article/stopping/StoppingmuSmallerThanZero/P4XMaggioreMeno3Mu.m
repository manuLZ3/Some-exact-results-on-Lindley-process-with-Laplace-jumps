function p = P4XMaggioreMeno3Mu(h,mu,s,x)
   
    %% coefficienti primo passo
    alpha_10 = exp((mu-h)/s)/2;
    
    %% coefficienti secondo passo, primo bin
    alpha_210 = (h-s)*exp((+mu)/s)*alpha_10/(2*s);
    eta_21 = alpha_10;

    %% coefficienti secondo passo, secondo bin
    alpha_220 = exp(mu/s)*alpha_10*(s*0.5+h)/(2*s);
    alpha_221 = -exp(mu/s)*alpha_10/(2*s);
    beta_220 = exp(-mu/s)*alpha_10*s*0.5/(2*s);

    %% coefficienti terzo passo, primo bin
    eta_31 = alpha_210+eta_21;
    alpha_310 = ( -mu*alpha_210 + (1-exp(mu/s))*s*eta_21 + (h+mu)*alpha_220 + (h+mu)^2/2*alpha_221 + s/2*(exp(2*mu/s)-exp(-2*h/s))*beta_220 );
    alpha_310 = exp(mu/s)/(2*s)*alpha_310-exp(mu/s)/2*(alpha_210+eta_21);

    %% coefficienti terzo passo, secondo bin
    alpha_320 = (h+mu)*alpha_220 + s/2*(exp(2*mu/s)-exp(-2*h/s))*beta_220 + (h+mu)^2/2*alpha_221 + s/2*alpha_210 - s*exp(mu/s)*eta_21;
    alpha_320 = alpha_320*exp(mu/s)/(2*s);
    alpha_321 = -alpha_210;
    alpha_321 = alpha_321*exp(mu/s)/(2*s);
    beta_320  = -s/2*alpha_210-s*eta_21;
    beta_320  = beta_320*exp(-mu/s)/(2*s) + exp(-mu/s)/2*(alpha_210 + eta_21); % l'ultima è la prob di fermarsi in 2 partendo da 0
    eta_32 = eta_21;

    %% coefficienti terzo passo, terzo bin
    alpha_330 = s/2*alpha_220 + (h+mu)*alpha_220 -(s/2)^2*alpha_221 + (h+mu)^2/2*alpha_221 - s/2*exp(-2*h/s)*beta_220;
    alpha_330 = alpha_330 *exp(mu/s)/(2*s);
    alpha_331 = s/2*alpha_221 - alpha_220;
    alpha_331 = alpha_331 *exp(mu/s)/(2*s);
    alpha_332 = -alpha_221/2;
    alpha_332 = alpha_332 *exp(mu/s)/(2*s);
    beta_330  = s/2*(exp(-2*mu/s)-1)*alpha_210+s*(exp(-mu/s)-1)*eta_21 - s/2*exp(-2*mu/s)*alpha_220 +(s/2)^2*exp(-2*mu/s)*alpha_221 + s/2*beta_220;
    beta_330  = beta_330*exp(-mu/s)/(2*s) + exp(-mu/s)/2*(alpha_210 + eta_21);  % l'ultima è la prob di fermarsi in 2 partendo da 0
    beta_331 = beta_220;
    beta_331 = beta_331*exp(-mu/s)/(2*s);
    
    % senza spezzare l'integrale in bin
    %p =-exp(-(mu + x)/s)*eta_31/2 - exp(-(mu + x)/s)*alpha_310/4 + eta_31*exp(-(2*mu + x)/s)/2 + alpha_310*exp(-(3*mu + x)/s)/4 + (-alpha_321*s^2*exp(-(5*mu + x)/s) - 2*alpha_321*mu*s*exp(-(3*mu + x)/s) + alpha_321*s^2*exp(-(3*mu + x)/s) + 2*alpha_320*s*exp(-(5*mu + x)/s) - 4*beta_320*exp(-(mu + x)/s)*mu - 4*eta_32*s*exp(-(2*mu + x)/s) + 4*eta_32*s*exp(-(3*mu + x)/s) - 2*alpha_320*s*exp(-(3*mu + x)/s))/(8*s) - (3*exp(-(5*mu + x)/s)*s^3*alpha_332 - 3*exp(-(5*mu + x)/s)*s^2*alpha_331 + 6*exp(-(5*mu + x)/s)*s*alpha_330 + 3*alpha_331*s^2*exp((mu + x)/s) - 3*alpha_332*s^3*exp((mu + x)/s) - 6*alpha_330*s*exp((mu + x)/s) + 12*alpha_330*exp((mu + x)/s)*mu + 12*alpha_330*exp((mu + x)/s)*x + 30*exp((mu + x)/s)*mu^2*alpha_331 + 6*exp((mu + x)/s)*x^2*alpha_331 + 76*exp((mu + x)/s)*mu^3*alpha_332 + 4*exp((mu + x)/s)*x^3*alpha_332 - 6*h^2*alpha_331*exp((mu + x)/s) - 12*alpha_330*h*exp((mu + x)/s) - 4*h^3*alpha_332*exp((mu + x)/s) - 54*beta_331*mu^2*exp(-(mu + x)/s) - 36*beta_331*mu*exp(-(mu + x)/s)*x - 18*beta_331*s*mu*exp(-(mu + x)/s) - 6*beta_331*s*x*exp(-(mu + x)/s) + 108*exp((mu + x)/s)*mu^2*x*alpha_332 - 36*alpha_332*mu*s*x*exp((mu + x)/s) - 54*mu^2*alpha_332*s*exp((mu + x)/s) + 18*mu*s^2*alpha_332*exp((mu + x)/s) + 6*s^2*x*alpha_332*exp((mu + x)/s) - 6*s*x^2*alpha_332*exp((mu + x)/s) - 18*alpha_331*s*mu*exp((mu + x)/s) - 6*alpha_331*s*x*exp((mu + x)/s) + 36*exp((mu + x)/s)*mu*x^2*alpha_332 + 36*alpha_331*mu*exp((mu + x)/s)*x - 24*h^2*mu*alpha_332*exp((mu + x)/s) - 48*mu^2*alpha_332*h*exp((mu + x)/s) - 24*mu*alpha_331*h*exp((mu + x)/s) + 6*beta_331*s*h*exp(-(2*h - mu - x)/s) + 12*beta_331*mu*s*exp(-(2*h - mu - x)/s) - 36*beta_330*exp(-(mu + x)/s)*mu - 12*beta_330*exp(-(mu + x)/s)*x - 6*exp(-(mu + x)/s)*x^2*beta_331 - 6*beta_330*s*exp(-(mu + x)/s) - 3*beta_331*s^2*exp(-(mu + x)/s) + 3*beta_331*s^2*exp(-(2*h - mu - x)/s) + 6*beta_330*s*exp(-(2*h - mu - x)/s))/(24*s) + exp(-(mu + x)/s)*(alpha_310 + eta_31)/2;
    % spezzandolo ma assumendo x >-h*mu>h ...assurdo
    %p=-exp(-(x + mu)/s)*eta_31/2 - exp(-(x + mu)/s)*alpha_310/4 + eta_31*exp(-(2*mu + x)/s)/2 + alpha_310*exp(-(3*mu + x)/s)/4 + (-s^2*alpha_321*exp(-(5*mu + x)/s) - 2*s*mu*alpha_321*exp(-(3*mu + x)/s) + s^2*alpha_321*exp(-(3*mu + x)/s) + 2*s*alpha_320*exp(-(5*mu + x)/s) - 4*s*eta_32*exp(-(2*mu + x)/s) + 4*s*eta_32*exp(-(3*mu + x)/s) - 2*s*alpha_320*exp(-(3*mu + x)/s) - 4*beta_320*exp(-(x + mu)/s)*mu)/(8*s) - (exp(-(5*mu + x)/s)*s^3*alpha_332 - 2*s*mu^2*alpha_332*exp(-(7*mu + x)/s) - 2*s^2*mu*alpha_332*exp(-(7*mu + x)/s) - s^3*alpha_332*exp(-(7*mu + x)/s) - 2*beta_331*mu^2*exp(-(x + mu)/s) - exp(-(5*mu + x)/s)*s^2*alpha_331 + 2*s*mu*alpha_331*exp(-(7*mu + x)/s) + s^2*alpha_331*exp(-(7*mu + x)/s) + 4*beta_330*exp(-(x + mu)/s)*mu + 2*exp(-(5*mu + x)/s)*s*alpha_330 - 2*s*alpha_330*exp(-(7*mu + x)/s))/(8*s) + (-3*s^3*alpha_332*exp(-(7*mu + x)/s) + 3*s^2*alpha_331*exp(-(7*mu + x)/s) - 6*s*alpha_330*exp(-(7*mu + x)/s) - 36*exp((x + mu)/s)*mu*x^2*alpha_332 - 76*exp((x + mu)/s)*mu^3*alpha_332 - 30*alpha_331*mu^2*exp((x + mu)/s) - 18*s^2*mu*alpha_332*exp((x + mu)/s) - 6*s^2*x*alpha_332*exp((x + mu)/s) + 6*s*x^2*alpha_332*exp((x + mu)/s) + 18*s*alpha_331*mu*exp((x + mu)/s) + 6*s*alpha_331*x*exp((x + mu)/s) + 48*beta_331*mu^2*exp(-(x + mu)/s) + 36*beta_331*mu*exp(-(x + mu)/s)*x - 6*s*beta_331*h*exp(-(2*h - mu - x)/s) - 12*s*beta_331*mu*exp(-(2*h - mu - x)/s) + 24*h^2*mu*alpha_332*exp((x + mu)/s) + 48*mu^2*alpha_332*h*exp((x + mu)/s) + 24*mu*alpha_331*h*exp((x + mu)/s) + 6*s*beta_331*x*exp(-(x + mu)/s) + 18*s*beta_331*mu*exp(-(x + mu)/s) - 36*exp((x + mu)/s)*mu*x*alpha_331 - 6*exp((x + mu)/s)*x^2*alpha_331 - 4*exp((x + mu)/s)*x^3*alpha_332 - 12*alpha_330*exp((x + mu)/s)*mu - 12*alpha_330*exp((x + mu)/s)*x - 3*s^2*alpha_331*exp((x + mu)/s) + 3*s^3*alpha_332*exp((x + mu)/s) + 6*s*alpha_330*exp((x + mu)/s) - 3*s^2*beta_331*exp(-(2*h - mu - x)/s) - 6*s*beta_330*exp(-(2*h - mu - x)/s) + 4*h^3*alpha_332*exp((x + mu)/s) + 6*h^2*alpha_331*exp((x + mu)/s) + 12*alpha_330*h*exp((x + mu)/s) + 48*beta_330*exp(-(x + mu)/s)*mu + 12*beta_330*exp(-(x + mu)/s)*x + 6*exp(-(x + mu)/s)*x^2*beta_331 + 3*s^2*beta_331*exp(-(x + mu)/s) + 6*s*beta_330*exp(-(x + mu)/s) + 54*s*mu^2*alpha_332*exp((x + mu)/s) + 36*s*mu*alpha_332*x*exp((x + mu)/s) - 108*exp((x + mu)/s)*mu^2*x*alpha_332 - 6*s*mu^2*alpha_332*exp(-(7*mu + x)/s) - 6*s^2*mu*alpha_332*exp(-(7*mu + x)/s) + 6*s*mu*alpha_331*exp(-(7*mu + x)/s))/(24*s) + exp(-(x + mu)/s)*(alpha_310 + eta_31)/2;
    p = -exp(-(x + mu)/s)*eta_31/2 - exp(-(x + mu)/s)*alpha_310/4 + eta_31*exp(-(2*mu + x)/s)/2 + alpha_310*exp(-(3*mu + x)/s)/4 - (2*alpha_321*mu*s*exp(-(3*mu + x)/s) - alpha_321*s^2*exp(-(3*mu + x)/s) + alpha_321*s^2*exp(-(5*mu + x)/s) + 4*eta_32*s*exp(-(2*mu + x)/s) - 4*eta_32*s*exp(-(3*mu + x)/s) + 2*alpha_320*s*exp(-(3*mu + x)/s) - 2*alpha_320*s*exp(-(5*mu + x)/s) + 4*beta_320*exp(-(x + mu)/s)*mu)/(8*s) - (3*exp((7*mu + x)/s)*beta_331*s^2 + 6*exp((7*mu + x)/s)*beta_330*s + 6*exp(-(5*mu + x)/s)*s*alpha_330 + 4*exp((x + mu)/s)*x^3*alpha_332 - 3*exp(-(5*mu + x)/s)*s^2*alpha_331 + 112*exp((x + mu)/s)*mu^3*alpha_332 + 6*exp((x + mu)/s)*x^2*alpha_331 + 3*exp(-(5*mu + x)/s)*s^3*alpha_332 + 12*alpha_330*exp((x + mu)/s)*x + 48*alpha_330*exp((x + mu)/s)*mu + 48*exp((x + mu)/s)*mu^2*alpha_331 - 12*beta_330*exp(-(x + mu)/s)*x - 6*exp(-(x + mu)/s)*x^2*beta_331 - 6*beta_330*s*exp(-(x + mu)/s) - 36*beta_330*exp(-(x + mu)/s)*mu - 3*beta_331*s^2*exp(-(x + mu)/s) - 3*alpha_332*s^3*exp((x + mu)/s) + 3*alpha_331*s^2*exp((x + mu)/s) - 6*alpha_330*s*exp((x + mu)/s) - 54*exp(-(x + mu)/s)*mu^2*beta_331 + 108*mu^2*alpha_332*exp((x + mu)/s)*x + 36*alpha_331*mu*exp((x + mu)/s)*x + 36*exp((x + mu)/s)*mu*x^2*alpha_332 - 6*beta_331*s*x*exp(-(x + mu)/s) - 18*beta_331*mu*s*exp(-(x + mu)/s) + 6*s^2*x*alpha_332*exp((x + mu)/s) - 36*exp(-(x + mu)/s)*mu*x*beta_331 + 18*alpha_332*mu*s^2*exp((x + mu)/s) - 54*mu^2*alpha_332*s*exp((x + mu)/s) - 6*alpha_331*s*x*exp((x + mu)/s) - 6*s*x^2*alpha_332*exp((x + mu)/s) - 18*alpha_331*mu*s*exp((x + mu)/s) - 6*exp((7*mu + x)/s)*beta_331*s*mu - 36*alpha_332*mu*s*x*exp((x + mu)/s))/(24*s) + (4*h^3*alpha_332*exp((x + mu)/s) + 24*h^2*mu*alpha_332*exp((x + mu)/s) + 48*mu^2*alpha_332*h*exp((x + mu)/s) + 36*exp((x + mu)/s)*mu^3*alpha_332 - 6*exp((7*mu + x)/s)*beta_331*s*mu + 3*exp((7*mu + x)/s)*beta_331*s^2 + 6*h^2*alpha_331*exp((x + mu)/s) + 24*mu*alpha_331*h*exp((x + mu)/s) + 18*exp((x + mu)/s)*mu^2*alpha_331 - 6*beta_331*s*h*exp(-(2*h - mu - x)/s) - 12*beta_331*mu*s*exp(-(2*h - mu - x)/s) - 3*beta_331*s^2*exp(-(2*h - mu - x)/s) + 6*exp((7*mu + x)/s)*beta_330*s + 12*alpha_330*h*exp((x + mu)/s) + 36*alpha_330*exp((x + mu)/s)*mu - 6*beta_330*s*exp(-(2*h - mu - x)/s))/(24*s) + exp(-(x + mu)/s)*(alpha_310 + eta_31)/2;
end