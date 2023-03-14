function p = Prob_N_conLuguale1(h,n,mu,s,w)
    %{      
    ------------ Old verison ------- WRONG
    if n==1
       p = 1-exp((h-mu-w)/s)/2;
    elseif n==2
       b_210 = exp((h-mu)/s)/2 - exp(-mu/s)/2 - h/(4*s)*exp((h-2*mu)/s) + exp(-h/s)/4 ;%+exp(-mu/s)/2 + exp((h-2*mu)/s)/4
       p = exp(-w/s)*b_210;
    else
       
       %{
       %% ricorsivo
       b_nmeno10 = P_N_L1(h,n-1,mu,s,0);
       b_n10 = b_nmeno10*exp(-mu/s)*( h/(2*s)+1/2 );
       %}

       %% froma chiusa
       b_210 = exp((h-mu)/s)/2 - exp(-mu/s)/2 - h/(4*s)*exp((h-2*mu)/s) + exp(-h/s)/4; %+exp(-mu/s)/2 + exp((h-2*mu)/s)/4
       b_n10 = b_210*(exp(-mu/s)*( h/(2*s)+1/2 ))^(n-2);
       %b_n10 = b_210*(exp(-(n-2)*mu/s)*( h/(2*s)+1/2 ))^(n-2);
       p = exp(-w/s)*b_n10;
    %}

       if n==1
            p = 1-exp((h-mu-w)/s)/2;
       elseif n==2
            b_2 = exp((h-mu)/s)/2 - ( h/(4*s) + 1/4 )*exp((h-2*mu)/s);
            p = exp(-w/s)*b_2;
       else 
            b_2 = exp((h-mu)/s)/2 - ( h/(4*s) + 1/4 )*exp((h-2*mu)/s);
            b_n = b_2*exp(-(n-2)*mu/s)*( h/(2*s)+1/2 )^(n-2);
            p = b_n*exp(-w/s);
       end
    
end