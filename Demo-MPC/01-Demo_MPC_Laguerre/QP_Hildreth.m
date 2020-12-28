%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function:  QP_Hildreth(E,F,M,Gama)
%%% Details:   Solver for quadratic objective functions with linear
%%%            constraints. QP_Hildreth finds a minimum for a problem
%%%            specified by CostFunction and its retraints with the
%%%            Hildreth method.Specially, the input F and Gama are the
%%%            vector of doubles.The solving speed is similar to the
%%%            quadprog function in the Matlab tools.
%%% CostFunc:  J = (1/2)x'Ex+F'x
%%% Subj to:   Mx<=Gama
%%% Date:      13-Nov-2020
%%% By:        Hongqian WEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, fval, Lambda] = QP_Hildreth(E,F,M,Gama)
% M=A_cons;
[n1,m1]=size(M);
X=-E\F;
fval = 0.5*X'*E*X+X'*F;
Lambda = 0;
kk=0;
for i=1:n1
if (M(i,:)*X>Gama(i)) kk=kk+1;
else
kk=kk+0;
end
end
if (kk==0) return; end

P = M*(E\M');
d = M*(E\F)+Gama;
[n, m] = size(d);
x_ini = zeros(n,m);
Lambda = x_ini;
al = 10;

for km=1:60
    lambda_p = Lambda;
    for i = 1:n
        w = P(i,:)*Lambda-P(i,i)*Lambda(i,1);
        w = w+d(i,1);
        la = -w/P(i,i);
        Lambda(i,1) = max(0,la);
    end
    al = (Lambda-lambda_p)'*(Lambda-lambda_p);
    if al<10e-8
        break;
    end
end
X = -E\F-E\M'*Lambda;
fval = 0.5*X'*E*X+X'*F; %Update the final cost value with the constraints
end
