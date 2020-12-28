function [Lzerot]=Mdu(a,N,n_in,Nc)
N_pa = sum(N);
M_du1 = zeros(n_in,N_pa);
k0=1;
[Al,L0]=LaguerreNetwork(a(k0),N(k0));
M_du1(1,1:N(1))=L0';
cc=N(1);
for k0=2:n_in;
    [Al,L0]=LaguerreNetwork(a(k0),N(k0));
    M_du1(k0,cc+1:cc+N(k0))=L0';    
    cc=cc+N(k0);
end
Lzerot=M_du1;
end