function a_ni0s = alpha_ni_muSmallerThan0(h,n,H,mu,s, A_n_pre,B_n_pre,C_n_pre)
    %{ ==============================================================
    %{ Calcolca i coefficienti a_ni0 che in P_i[N=n] moltiplicano e^u
    %{ ==============================================================
    
    % prob di fermarsi in n-1 passi partendo da 0
    m_nMeno1_1 = (n-2);%min(n-1,H-1+1)-1; %min(n,L-i+1)-1;
    %% monomi
    p = 0:m_nMeno1_1;
    monomi = transpose((0+(n-2)*mu).^p);
    K0 = ( mtimes(A_n_pre(1,:), monomi) + mtimes(B_n_pre(1,:), monomi) ) + C_n_pre(1);
    
    a_ni0s = zeros(H,1);
    %% calcolo dei 'a_{n,i,0}'
    % il ciclo va fino a 'L-1' perchè 'a_{n,L,0}' si calcola diversamente
    for i = 2:H
      A = 0;
      i_star = i-1;
      for j = (i_star+1):H
         % estremi di integrazione
         Ij_l = -(j-1)*mu;    
         Ij_r =  min(-j*mu,h);     
         m_nmeno1_j = (n-2); %min(n-1,H-j+1)-1;
         % coefficienti in P_j[N=n+1]
         b_nMeno1j = B_n_pre(j,:);
         a_nMeno1j = A_n_pre(j,:);
         c_nMeno1j = C_n_pre(j);
         %% Calcolo dei K_j^(2). Formula (4.74)
         for p = 0:m_nmeno1_j 
               %integrand = generalizedGamma(p,s/2,(n-1)*mu);
               %gamma_jp = integrand.evaluate_integral(Ij_l,Ij_r); 
               %fun = @(y,p,n,k,s) (y+(n-2)*k).^p.*exp(-2.*y./s);
               %gamma_jp = integral(@(y) fun(y,p,n,mu,s), Ij_l, Ij_r);
               gamma_jp = gammaMia(Ij_l, Ij_r ,p, s/2, (n-2)*mu);
               power_jp = ((Ij_r+(n-2)*mu)^(p+1)-(Ij_l+(n-2)*mu)^(p+1))/(p+1);
               A = A + a_nMeno1j(p+1)*power_jp + b_nMeno1j(p+1)*gamma_jp;   
         end
         %A = A - c_nj*(2*s)^(n-1)*(exp(-Ij_r/s)-exp(-Ij_l/s));
         A = A - c_nMeno1j*s*(exp(-Ij_r/s)-exp(-Ij_l/s));
      end
       
     if H > 1 %&& i+1 <= H
          %estremi di integrazione
          Iistar_l = -(i_star-1)*mu;  
          Iistar_r = min(-(i_star)*mu,h);        
          %m_ipiu1n = min(L,n)-1;
          m_ipiu1_nmeno1 = (n-2); %min((n-1),H-i_star+1)-1; %max(min(n-1,L-(i+1)+1)-1,0);
          % coefficienti in P_j[N=n+1]
          b_nMeno1istar = B_n_pre(i_star,:);
          a_nMeno1istar = A_n_pre(i_star,:);
          c_nMeno1istar = C_n_pre(i_star);  
          for p = 0:m_ipiu1_nmeno1                                                                   
              power_istarp = ((Iistar_r+(n-2)*mu)^(p+1))/(p+1);
              gamma_istarp = 0;
              for d=0:p
                    % WRONG!
                    %gamma_ip = gamma_ip - (Iipiu1_r+(n-1)*mu)^(p-d)*exp(-2*(Iipiu1_r+(n-1)*mu)/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1);
                    gamma_istarp = gamma_istarp - (Iistar_r+(n-2)*mu)^(p-d)*exp(-2*(Iistar_r)/s)*factorial(p)/factorial(p-d)*(s/2)^(d+1);
              end                   
              A = A + a_nMeno1istar(p+1)*( power_istarp + (-1)^(p)*factorial(p)*(s/2)^(p+1) ) + b_nMeno1istar(p+1)*gamma_istarp;
          end
          %A = A - c_nistar*(2*s)^(n-1)*exp(-Iistar_r/s);
          A = A - c_nMeno1istar*s*exp(-Iistar_r/s);
     end
    %A = A*exp(mu/s); 
    A = A*exp(mu/s)/(2*s);
 
    a_ni0s(i) = A;
    end
    
    %% formula (4.77g)
    A = 0;
    for j = 1:H
        % estremi di integrazione
        Ij_l = -(j-1)*mu;    
        Ij_r =  min(-j*mu,h);     
        m_nmeno1_j = (n-2); %min(n-1,H-j+1)-1;
        % coefficienti in P_j[N=n+1]
        b_nMeno1j = B_n_pre(j,:);
        a_nMeno1j = A_n_pre(j,:);
        c_nMeno1j = C_n_pre(j);
        % il ciclo va fino a '-1' perchè si indicizza da 'p+1': 'b_nj(1)' è il coefficiente della potenza 0
        for p = 0:m_nmeno1_j 
              %integrand = generalizedGamma(p,s/2,(n-1)*mu);
              %gamma_jp = integrand.evaluate_integral(Ij_l,Ij_r); 
              %fun = @(y,p,n,k,s) (y+(n-2)*k).^p.*exp(-2.*y./s);
              %gamma_jp = integral(@(y) fun(y,p,n,mu,s), Ij_l, Ij_r);
              gamma_jp = gammaMia(Ij_l, Ij_r ,p, s/2, (n-2)*mu);
              power_jp = ((Ij_r+(n-2)*mu)^(p+1)-(Ij_l+(n-2)*mu)^(p+1))/(p+1);
              A = A + a_nMeno1j(p+1)*power_jp + b_nMeno1j(p+1)*gamma_jp;   
         end
         %A = A - c_nj*(2*s)^(n-1)*(exp(-Ij_r/s)-exp(-Ij_l/s));
         A = A - c_nMeno1j*s*(exp(-Ij_r/s)-exp(-Ij_l/s));
    end
    a_ni0s(1) = A*exp(mu/s)/(2*s)-K0*exp(mu/s)/(2);
end