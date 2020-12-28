%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function:  LaguerreNetwork(a,N)
%%% Details:   Laguerre functions to resolve the parameter matrix A and L(0)
%%% Input:     scaling factor a and scaling number N
%%% Output:    parametric matrix A and the initial Laguerre vector L(0)
%%% Details:   for example N=5 the matrix A is represented as follow.    
%%%            L(0)=sqrt(beta)*[1 -a a^2 -a^3 a^4]';
%%%            A=[  a          0         0       0    0;
%%%                beta        a         0       0    0
%%%               -a*beta     beta       a       0    0
%%%               a^2*beta   -a*beta    beta     a    0
%%%              -a^3*beta   a^2*beta  -a*beta  beta  a ]
%%% Date:      12-Dec-2020
%%% By:        Hongqian WEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,L0]=LaguerreNetwork(a,N)
v(1,1)=a;
L0(1,1)=1;
for k=2:N
v(k,1)=(-a).^(k-2)*(1-a*a);
L0(k,1)=(-a).^(k-1);
end
L0=sqrt((1-a*a))*L0;
A(:,1)=v;
for i=2:N
A(:,i)=[zeros(i-1,1);v(1:N-i+1,1)];
end

