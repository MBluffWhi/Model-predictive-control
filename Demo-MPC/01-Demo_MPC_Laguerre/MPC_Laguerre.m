function [sys,x0,str,ts] = MPC_Laguerre(t,x,u,flag)
%%%         This is a demo of model predictive control based on LaguerreNetwork.
%%%         The spring-damping system is used to be the demo plant. The argument model
%%%         contains the data structures needed by the controller. The input 
%%%         to the S-function block is a vector signal consisting of the
%%%         measured outputs and the reference values. The output of the 
%%%         S-function block is a vector signal consisting of the control 
%%%         variables and the estimated state vector,potentially including estimated 
%%%         disturbance states.
%%% Date:   13-Nov-2020
%%% By:     Hongqian WEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch flag,
 case 0
  [sys,x0,str,ts] = mdlInitializeSizes; % Initialization
  
 case 2
  sys = mdlUpdates(t,x,u); % Update discrete states
  
 case 3
  sys = mdlOutputs(t,x,u); % Calculate outputs
 


 case {1,4,9} % Unused flags
  sys = [];
  
 otherwise
  error(['unhandled flag = ',num2str(flag)]); % Error handling
end
% End of dsfunc.

%==============================================================
% Initialization
%==============================================================

function [sys,x0,str,ts] = mdlInitializeSizes

% Call simsizes for a sizes structure, fill it in, and convert it 
% to a sizes array.

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 3;
sizes.NumOutputs     = 1;
sizes.NumInputs      = 3; %S-function的输入包含状态变量(2个)和输出的参考(1个)
sizes.DirFeedthrough = 1; % Matrix D is non-empty.
sizes.NumSampleTimes = 1;
sys = simsizes(sizes); 
x0 =[0;0;0];   
global U StateVecTemp;
U = 0;
StateVecTemp = [0; 0];
% Initialize the discrete states.
str = [];             % Set str to an empty matrix.
ts  = [0.05 0];       % sample time: [period, offset]
%End of mdlInitializeSizes
		      
%==============================================================
% Update the discrete states
%==============================================================
function sys = mdlUpdates(t,x,u)
  
sys = x;
%End of mdlUpdate.

%==============================================================
% Calculate outputs
%==============================================================
function sys = mdlOutputs(t,x,u)
    global delta_u;
    global U;
    global StateVecTemp;
    global Dis
    tic
    Ts=0.05;
    Np =40; % predictive horizon
    %% 连续模型的离散化处理
    MtxAc = [0 1;-1 -0.5]; MtxBc = [0; 1];MtxCc = [1 0]; MtxDc=zeros(1,1);
%     [MtxAd,MtxBd,MtxCd,MtxDd]=c2dm(MtxAc,MtxBc,MtxCc,MtxDc,Ts) % Eula
%     front differential equation
    [m1,n1]=size(MtxCc);
    [n1,n_in]=size(MtxBc);
    Nx = n1; %Nx: number of state variable 
    Nu = n_in; %Nu: number of control variable 
    Ns = m1; %Ns: number of observed/system output variable 
    MtxAd = MtxAc*Ts+eye(Nx);
    MtxBd = MtxBc*Ts; MtxDd=MtxDc*Ts; MtxCd = MtxCc;
    
    %%Augmented model
    MtxA=eye(n1+m1,n1+m1);  % MtxA represented Matrix A
    MtxA(1:n1,1:n1)=MtxAd;
    MtxA(n1+1:n1+m1,1:n1)=MtxCd*MtxAd;
    MtxB=zeros(n1+m1,n_in);
    MtxB(1:n1,:)=MtxBd;
    MtxB(n1+1:n1+m1,:)=MtxCd*MtxBd;
    MtxC=zeros(m1,n1+m1);
    MtxC(:,n1+1:n1+m1)=eye(m1,m1);
    [n, n_in] = size(MtxB); % n: number of state variable in agumented model
    
    q = 200;
    a=0.8;
    N=4;
    R = 0.5;
    C_bar = eye(Ns).*q;
    C_td = [zeros(Ns,Nx) C_bar];
    Q_bar = sqrt(C_td'*C_td);     
%     R_bar = kron(R,eye(N));         
    R_bar = 0.5;
    %%processing state variable
    StateVec = zeros(Nx,1);
    StateVec(1) = u(1);
    StateVec(2) = u(2);
    r = u(3); %referece value R_vect=[r1(k) r2(k) ... rq(k)]^(T) 
    y = MtxCc*StateVec; % system output value
    % x(ki) = [delta_x(k);y_td]  ;%y_td = y(ki)-r(ki)  ki=instant(k)
    Xf = [StateVec-StateVecTemp;y-r];  
    
    [Omega,Psi] = Objective_MPC_Laguerre(MtxA,MtxB,MtxC,MtxDd,Xf,a,N,Np,Q_bar,R_bar);
    Nc = 2;        %Nc represented the number of constaints on control input
    Lzerot = Mdu(a,N,n_in,Nc);  % obtain the first control increment
    %% Quadratic programing construction
    %cost function: J = (1/2)*(Eta^T)*H*Eta+(Eta^T)*f;
    M = 0; Gama=0;    %constraints of demo
    eta = QP_Hildreth(Omega,Psi,M,Gama);  %QP solver developed according to Hildreth method
   %% Output of s-function and previous state saving
    StateVecTemp = StateVec;
    delta_u = Lzerot*eta;
    U(1) = U(1)+delta_u;
    sys = U(1);
    toc
% End of mdlOutputs.