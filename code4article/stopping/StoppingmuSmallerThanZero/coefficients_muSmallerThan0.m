function [A_n,B_n,C_n] = coefficients_muSmallerThan0(h,n,mu,s)
   

   %% Matrici dei coefficienti:
   % -------------------------
   %       i° riga : intervallo i nella partizione
   %    p° colonna : coefficiente di (y-nk)^p
   %    dimensione : (L X \bar{m}_{n,1})
    
   % smaller integer such that ´L*k>h´. Number of bins in the partition
   H = ceil(h/-mu);
    
   % max degree of polynomial among the functions P_i(n,w) e P_i(n-1,w)
   m_n_1 = (n-1);%min(n,H)-1;  % min(n,L+1)-1 
   bar_m_n_1 = m_n_1+1; % si aggiunge uno per il termine noto
   
  
   %% caso base
   if n == 1
       A_n = transpose(repelem(exp((-h+mu)/s)/2,H));
       B_n = transpose(repelem(0,               H));
       C_n = transpose(repelem(0,               H));
       
   else
   
       % crea coefficienti al passo precedente 
       [A_n_pre,B_n_pre,C_n_pre] = coefficients_muSmallerThan0(h,n-1,mu,s);
       
       % prob di fermarsi in n-1 passi partendo da 0
       % grado massimo
       m_nMeno1_1 = (n-2); %min(n-1,H-1+1)-1; %min(n,L-i+1)-1;
       %% monomi
       p = 0:m_nMeno1_1;
       monomi = transpose((0+(n-2)*mu).^p);
       K0 = ( mtimes(A_n_pre(1,:), monomi) + mtimes(B_n_pre(1,:), monomi) ) + C_n_pre(1);
      
       % inizializza matrici dei coefficienti
       A_n = zeros([H,bar_m_n_1]);
       B_n = zeros([H,bar_m_n_1]);
       C_n = zeros([H,1]);
       for i = 2:H
           i_star = i-1;
           m_nMeno1_i_star = n-1;%min(n-1,H-i_star+1);
           for j = 1:bar_m_n_1-1
               a=0;
               b=0;
               for k = (j-1):m_nMeno1_i_star-1
                   %% formula (4.5)
                   a = a + A_n_pre(i_star,k+1)*(-s/2)^(k-j+1)*factorial(k)/factorial(j);
                   %% formula (4.7)
                   b = b + B_n_pre(i_star,k+1)*(s/2)^(k-j+1)*factorial(k)/factorial(j);
               end
               %A_n(i,j+1) = -exp(mu/s)*a;
               %B_n(i,j+1) = exp(-mu/s)*b;
               A_n(i,j+1) = -exp(mu/s)*a/(2*s);
               B_n(i,j+1) = exp(-mu/s)*b/(2*s);
           end
           %if i < H
           C_n(i) = C_n_pre(i_star);
           %end
       end
       % il caso ´i=1´ è trattato separatamente 
       C_n(1) = K0;
        

       %  calcola i coefficienti per i=0 (prima colonna)
       a_ni0s = alpha_ni_muSmallerThan0(h,n,H,mu,s, A_n_pre,B_n_pre,C_n_pre);
       b_ni0s = beta_ni_muSmallerThan0(h,n,H,mu,s, A_n_pre,B_n_pre,C_n_pre);
       % colloca i termini noti a_ni0 sulla prima colonna
       A_n(:,1) = a_ni0s;
       B_n(:,1) = b_ni0s;
   end
end




   