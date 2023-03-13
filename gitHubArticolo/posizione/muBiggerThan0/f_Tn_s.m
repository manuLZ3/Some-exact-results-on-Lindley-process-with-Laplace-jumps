function y = f_Tn_s(u,n,k,x,s, coefficients) 

    %% create coefficients
    if ~exist('coefficients','var')
        [A_n,B_n,c_n] = f_Tn_coefficients_s(n,k,x,s);
    else 
        [A_n,B_n,c_n] = coefficients;
    end

    %% identify the bin of partition at which `u` belongs to
    found=false;
    if u==0
        y = c_n;
    else
        for i=1:n
            if ((i-1)*k<u) && (u<=i*k)
                found=true;
                break
            end
        end
        if not(found)
            if ((n-1)*k<u) && (u<=n*k+x)
                i = n+1;
            else
                i = n+2;
            end
        end
        
        %% compute output
        % max degree of monomila
        m = min(n,i)-1;
        % create monomial functions
        monomi = [];
        for p=0:m
            monomi = [monomi, (u-n*k)^p];
        end
        if n>i
            monomi = [monomi, repelem(0,n-i)];
        end
        % output
        y = dot(monomi, A_n(:,i))*exp(u/s) + dot(monomi, B_n(:,i))*exp(-u/s);
    end
end

%{
function y = f_T1(u,x,k,s)
if u < 0
   y = 0;
elseif u == 0
   y = exp((-x-k)/s)/2;
elseif u < k+x
   y = exp((u-x-k)/s)/(2*s);
else
   y = exp((-u+x+k)/s)/(2*s);
end
end
%}