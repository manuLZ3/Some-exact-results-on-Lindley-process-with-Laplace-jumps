function  P = Prob_N_symbolicoConLuguale1(h, n, k, s)


    %% Definizione variabili
    w = sym('w');
    x = sym('x');
    assume(w>0);
    assume(x>0); 

    %% Calcolo integrale delle probabilità
    P = [];
    P_n = [];

    if n == 1
        P_n = [1-exp((h-k-w)/s)/2];
        P = P_n;
    else
        P_pre = Prob_N_symbolicoConLuguale1(h, n-1, k, s);
        %P_nmeno1 = P_pre(length(P_pre));
        P_nmeno1 = P_pre(end);
        % funzione del passo precedente da integrare
        Pnmeno1_integranda = subs(P_nmeno1,w,x);
        % funzione del passo precedente valutata in 0
        Pnmeno1_da0 = subs(P_nmeno1,w,0);

        % estremi di integrazione 
        I_l = 0;
        I_r = h;
        P_n = int(Pnmeno1_integranda*exp((x-k-w)/s)/(2*s), x, I_l, I_r)+ Pnmeno1_da0*exp((-k-w)/s)/2;

        P = [P_pre;P_n];
    end
end
