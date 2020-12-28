function [sys,x0,str,ts] = MPC_PulseCtl(t,x,u,flag)
%%%         This is a demo of model predictive control based on pulse-operator-based input sequense.
%%%         The spring-damping system is used to be the demo plant. The argument model
%%%         contains the data structures needed by the controller. The input 
%%%         to the S-function block is a vector signal consisting of the
%%%         measured outputs and the reference values. The output of the 
%%%         S-function block is a vector signal consisting of the control 
%%%         variables and the estimated state vector,potentially including estimated 
%%%         disturbance states.
%%% Date:   9-Sep-2020
%%% By:     Hongqian WEI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
sizes.NumInputs      = 3; %S-function input contains 2state variables and 1 reference value.
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
    %% 连续模型的离散化处理
    MtxAc = [0 1;-4 0]; MtxBc = [1; 0];MtxCc = [0 1]; MtxDc=zeros(1,1);
    [MtxAd,MtxBd,MtxCd,MtxDd]=c2dm(MtxAc,MtxBc,MtxCc,MtxDc,Ts) 
    [m1,n1]=size(MtxCd);
    [n1,n_in]=size(MtxBd);
    Nx = n1; 
    Nu = n_in; 
    Ns = m1; 
    %%Augmented model
    MtxA=eye(n1+m1,n1+m1);  % MtxA represented Matrix A
    MtxA(1:n1,1:n1)=MtxAd;
    MtxA(n1+1:n1+m1,1:n1)=MtxCd*MtxAd;
    MtxB=zeros(n1+m1,n_in);
    MtxB(1:n1,:)=MtxBd;
    MtxB(n1+1:n1+m1,:)=MtxCd*MtxBd;
    MtxC=zeros(m1,n1+m1);
    MtxC(:,n1+1:n1+m1)=eye(m1,m1);
    [n, n_in] = size(MtxB) 
    
    %%处理状态变量
    StateVec = zeros(Nx,1);
    StateVec(1) = u(1);
    StateVec(2) = u(2);
    r = u(3); % outpur reference R_vect=[r1(k) r2(k) ... rq(k)]^(T) 
    y = MtxCd*StateVec; % system output
    Xf = [StateVec-StateVecTemp;y];  %Xf is the key to solve cost function
    
    
    Np =20; %predictive horizon
    Nc=6;   %control horizon

    F_cell = cell(Np,1);
    PHI_cell = cell(Np,Nc);
    W_cell = cell(Np,1);
 
    for i=1:Np
      F_cell{i,1} = MtxC*(MtxA^i);
      
       %{ 
      if i==1
          W_cell{i,1} = MtxC*MtxA^(i-1)*MtxD;
      else
          W_cell{i,1} = W_cell{i-1,1}+MtxC*MtxA^(i-1)*MtxD;  %recursive approach
      end
      %}
      
      for j = 1:Nc
         if j<=i
              PHI_cell{i,j} = MtxC*(MtxA^(i-j))*MtxB;
         else
              PHI_cell{i,j} = zeros(Ns,Nu);
         end
      end
    end
    PHI = cell2mat(PHI_cell);   
    F= cell2mat(F_cell); 
%     W = cell2mat(W_cell);
    Rs = ones(Np,1);
    Rsbar = kron(Rs,r);     % Expand reference r to total predictive horizon with Kronecker production 
    Rbar = 0.5*eye(Nu*Nc);  % weight coefficient for control increment
    H = 0.5*PHI'*PHI+Rbar;   %#ok<MHERM>
    f = -2*PHI'*(Rsbar-F*Xf);
    %% Constraints
    % inequality constraints: the demo is not added.
    %% Solver
    %cost function: J = (1/2)*(Eta^T)*H*Eta+(Eta^T)*f;
    M = 0; Gama=0;
    X = QP_Hildreth(H,f,M,Gama);  %Hildreth's QP solver.
   %% S-function output and previous state saving
    StateVecTemp = StateVec;
    delta_u = X(1);
    U(1) = U(1)+delta_u;
    sys = U(1);
    toc
% End of mdlOutputs.