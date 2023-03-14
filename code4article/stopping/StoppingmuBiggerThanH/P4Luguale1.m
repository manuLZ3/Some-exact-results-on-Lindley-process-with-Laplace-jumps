function p4 = P4Luguale1(h,mu,s,x)
     b_2 = exp((h-mu)/s)/2 - ( h/(4*s) + 1/4 )*exp((h-2*mu)/s);
     p4 = b_2*( h/(2*s) + 1/2 )^2*exp((-2*mu)/s)*exp(-x/s);
end