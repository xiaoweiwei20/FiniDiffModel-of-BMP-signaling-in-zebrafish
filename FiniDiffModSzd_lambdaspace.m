% function name: FiniDiffModSzd_lambdaspace.m
% author: Xu Wang and David Umulis 

% functionality:
%    This function uses a 1D finite difference method to solve the ODE for
%    Zebrafish BMP gradient formation during embryogenesis including sizzled, bmp1a and tolloid. 
%    Based on the primary function FiniteTest.m for the evaluation of differential 
%    equations using the 4th order runge-kutta technique with a fixed step-size.
%
% function input:
%    n - number of nodes to evaluate finite difference equations
%    tRange - time interval to evaluate differential equations, eg. [0 5000]
%    parameters - a stucture storing input parameter values like diffusion
%                 rates, and reaction rates, etc. 
%        .Ltot - length of the embryo (1D assumption)       microns
%        .Lven - length of ventral region (1D assumption)   microns
%        .LdorC - length of dorsal Chd expression region    microns
%        .LvenB - length of dorsal bmp expression region    microns
%        .LdorN - length of dorsal Nog expression region    microns
%        .LvenS - length of dorsal sizzled expression region    microns
%        .LvenT - length of dorsal tolloid expression region    microns
%        .DB - diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
%        .DC - diffusion rate of Chordin       (microns^2*s^-1)*60s*m^-1
%        .DN - diffusion rate of Noggin          (microns^2*s^-1)*60s*m^-1
%        .DS - diffusion rate of sizzled          (microns^2*s^-1)*60s*m^-1
%        .DBC - diffusion rate of [BMP,Chd]    (microns^2*s^-1)*60s*m^-1
%        .DBN - diffusion rate of [BMP,Nog]      (microns^2*s^-1)*60s*m^-1
%        .k1 - binding rates for BMP ligand and Chordin      nM^-1*m^-1
%        .k_1 - unbinding rates for BMP ligand and Chordin   m^-1
%        .k2 - binding rates for BMP ligand and Noggin       nM^-1*m^-1
%        .k_2 - unbinding rates for BMP ligand and Noggin    m^-1
%        .k3 - binding rates for sizzled and tolloid nM^-1*m^-1
%        .k_3 - unbinding rates for sizzled and tolloid  m^-1 
%        .decB - decay rate of Ligand (BMP)    m^-1  
%        .decC - decay rate of Chd             m^-1  
%        .decN - decay rate of Nog             m^-1  
%        .decS - decay rate of Sizzled         m^-1 
%        .decT - decay rate of Tolloid         m^-1 
%        .j1 - production rate of BMP          nM*m^-1  
%        .j2 - production rate of Chordin      nM*m^-1
%        .j3 - production rate of Noggin       nM*m^-1
%        .j4 - production rate of sizzled       nM*m^-1
%        .j5 - production rate of tolloid       nM*m^-1
%        .lambda1 - tld processing rate of BChd     nM^-1*m^-1
%        .lambda2 - tld processing rate of Chd   nM^-1*m^-1
%     save_or_not - whether to save the output into .mat files (0 or 1)
%                   default = 0 (not save)
% NOTE: make n that can exactly divided by Ltot for easy calculation
%
% function output:
%    B - vector storing distribution of BMP 
%    C - vector storing distribution of Chordin
%    N - vector storing distribution of Noggin
%    S - vector storing distribution of Noggin
%    BC - vector storing distribution of BMP-Chordin complex
%    BN - vector storing distribution of BMP-Noggin complex
%% Extract parameter
function [B, C, N, S,BC, BN,time] = FiniDiffModSzd_lambdaspace(n, tRange, parameters)

tic; % start record time

% Create a function that handles faster evaluation of differential equations
fdAdt = @dAdt;

%Capture current time to measure computational speed
%----------------- initialize geometry and vectors --------------------
B0 = zeros(1,n);   % initialize BMP vector
C0 = zeros(1,n);   % initialize Chordin vector
N0 = zeros(1,n);   % initialize Noggin vector
S0 = zeros(1,n);   % initialize Sizzled vector
BC0 = zeros(1,n);  % initialize BMP_Chordin vector
BN0 = zeros(1,n);  % initialize BMP_Noggin vector

initial = [B0, C0, N0, S0, BC0, BN0];    % initialize condition vector

options = odeset('RelTol', 1e-9);    % set solving options relative error tolerance
% Set to force parameters
% solve ode using 15s
[time, D] = ode15s(fdAdt, tRange, initial, options, n,parameters);   

%------------------------- store data in vectors --------------------------
B = D(:, 1:n);      
C = D(:, n+1:2*n);
N = D(:, 2*n+1:3*n);
S = D(:, 3*n+1:4*n);
BC = D(:, 4*n+1:5*n);
BN = D(:, 5*n+1:6*n);

%--------------------------- make x data vector --------------------------- 
start = -parameters.Ltot + parameters.Ltot/n;
X = start:(parameters.Ltot/n):parameters.Ltot;

end

% This fucntion is the set of differential equations
function dY = dAdt(t, Y, n, parameters) 

%----------------------- obtain parameter values --------------------------
  %%% geometry
Ltot = parameters.Ltot;   % length of the embryo (1D assumption)       microns
% Lven = parameters.Lven;  % length of ventral region for Bmp expression         microns
% LvenTld = parameters.LvenTld; % length of ventral region for Tld expression    microns
% LvenSzd = parameters.LvenSzd; % length of ventral region for Sizzled     microns
LdorC = parameters.LdorC;   % length of dorsal Chd expression region   microns
LdorN = parameters.LdorN;   % length of dorsal Nog expression region   microns

  %% diffusion rates
DB = parameters.DB;       % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
DC = parameters.DC;       % diffusion rate of Chordin         (microns^2*s^-1)*60s*m^-1
DN = parameters.DN;       % diffusion rate of Noggin          (microns^2*s^-1)*60s*m^-1
DS = parameters.DS;       % diffusion rate of Sizzled        (microns^2*s^-1)*60s*m^-1
DBC = parameters.DBC;     % diffusion rate of [BMP,Chd]       (microns^2*s^-1)*60s*m^-1
DBN = parameters.DBN;     % diffusion rate of [BMP,Nog]       (microns^2*s^-1)*60s*m^-1

  %% kinetic parameters
k1 = parameters.k1;       % binding rates for BMP ligand and Chordin          nM^-1*m^-1
k_1 = parameters.k_1;     % unbinding rates for BMP ligand and Chordin        m^-1
k2 = parameters.k2;       % binding rates for BMP ligand and Noggin           nM^-1*m^-1
k_2 = parameters.k_2;     % unbinding rates for BMP ligand and Noggin         m^-1
kmt = parameters.kmt;       % Michaelis constant of the proteinase Tolloid  nM
kit = parameters.kit;       % the inhibitor constant(dissociation constants) of the protianse tolloid  nM
kma = parameters.kma;       % Michaelis constant of the proteinase bmp1a  nM
kia = parameters.kia;       % the inhibitor constant(dissociation constants) of the protianse bmp1a nM

  %% Positive Feedback
Vs = parameters.Vs;         % Vmax of the BMP on Sizzled
k = parameters.k;         % parameter for the hill function
nu = parameters.nu;       % cooperative parameter

  %% decay/recycle rates 
decB = parameters.decB;   % decay rate of Ligand (BMP)    nM*m^-1  
decC = parameters.decC;   % decay rate of Chd             nM*m^-1  
decN = parameters.decN;   % decay rate of Nog             nM*m^-1 
decS = parameters.decS;   % decay rate of Sizzled             nM*m^-1
decBC = parameters.decBC;   % decay rate of BC    nM*m^-1  
decBN = parameters.decBN;   % decay rate of BN            nM*m^-1  
  %% production rates 
j1 = parameters.j1;       % production rate of BMP          nM*m^-1  
j2 = parameters.j2;       % production rate of Chordin      nM*m^-1
j3 = parameters.j3;       % production rate of Noggin       nM*m^-1

  %% Tolloid behavior
lambda_tld_C = parameters.lambda_tld_C;     % the max velocity tld processing rate of Chd  nM^-1*m^-1
lambda_tld_BC = parameters.lambda_tld_BC;   % the max velocity tld processing rate of BC   nM^-1*m^-1

lambda_bmp1a_C = parameters.lambda_bmp1a_C;     % the max velocity bmp1a processing rate of Chd  nM^-1*m^-1
lambda_bmp1a_BC = parameters.lambda_bmp1a_BC;   % the max velocity bmp1a processing rate of BC   nM^-1*m^-1

%------------- Convert boundaries to positions of nearest node ------------
% nvenSzd = round(LvenSzd*n/Ltot);
ndorC = round(LdorC*n/Ltot);
ndorN = round(LdorN*n/Ltot);
%--- Increment size for finite difference method,distance between nodes --- 
dx = Ltot/(n-1);
%----------------- Initialize prepatterns for secretion -------------------
etaC = zeros(1,n);
etaN = zeros(1,n);
% etaT = zeros(1,n);
% etaS = zeros(1,n);
x_lig = 0:dx:Ltot;% change x from space to node
yspace = (exp(-((x_lig+dx)/20-11)/5))./(0.1+exp(-((x_lig+dx)/20-11)/5));
etaB = j1 * yspace;      % define ventral region for BMP production
lambda_tld_BC_space = lambda_tld_BC * yspace;
lambda_tld_C_space = lambda_tld_C * yspace;
% plot((1:n),yspace)
p=1;
qa = (t-6480)./(2000+abs(t-6480));
q = qa.*(qa>0);
for i = 1:n
    etaC(i) = j2 * double( i>(n-ndorC) );  % define dorsal region for Chordin production
%     etaC(i) = j2 * chdm(i);
    etaN(i) = j3 * double( i>(n-ndorN) );  % define dorsal region for Noggin production
%    Tld(i) = tld_conc;       % Tolloid uniformly distributed in the embryo
%     etaT(i) = j4 * double(i<=nvenTld);   % Tolloid distributed ventrally
end

%------------ convert Y vector into vectors for each component ------------
Bmp = Y(1:n);
Chd = Y(n+1:2*n);
Nog = Y(2*n+1:3*n);
Szd = Y(3*n+1:4*n);
BC = Y(4*n+1:5*n);
BN = Y(5*n+1:6*n);

%----------------------- Zero out difference vectors ----------------------
dBmp = zeros(1,n);
dChd = zeros(1,n);
dNog = zeros(1,n);
dSzd = zeros(1,n);
dBC = zeros(1,n);
dBN = zeros(1,n);

%--------------------Begin solving for ODEs at each time point--------
%--------------------Index one corresponds to ventral midline---------

dBmp(1) = DB/dx^2*(-2*Bmp(1)+2*Bmp(2)) ...
        - k1*Bmp(1)*Chd(1) + k_1*BC(1) ...
        - k2*Bmp(1)*Nog(1) + k_2*BN(1) ...
        - decB*Bmp(1) + etaB(1) + q*lambda_tld_BC_space(1)*BC(1)/(1+Szd(1)/kit+(Chd(1)+BC(1))/kmt)+p*lambda_bmp1a_BC*BC(1)/(1+Szd(1)/kia+(Chd(1)+BC(1))/kma);
    
dChd(1) = DC/dx^2*(-2*Chd(1)+2*Chd(2)) ...
        - k1*Bmp(1)*Chd(1) + k_1*BC(1) ...
        - decC*Chd(1) + etaC(1) - q*lambda_tld_C_space(1)*Chd(1)/(1+Szd(1)/kit+(Chd(1)+BC(1))/kmt) - p*lambda_bmp1a_C*Chd(1)/(1+Szd(1)/kia+(Chd(1)+BC(1))/kma);
      
dNog(1) = DN/dx^2*(-2*Nog(1)+2*Nog(2)) ...
        - k2*Bmp(1)*Nog(1) + k_2*BN(1) ...
        - decN*Nog(1) + etaN(1);
        
dSzd(1) = DS/dx^2*(-2*Szd(1)+2*Szd(2))...
        - decS*Szd(1) + Vs*Bmp(1)^nu/(k^nu+Bmp(1)^nu);
      
dBC(1) = DBC/dx^2*(-2*BC(1)+2*BC(2)) ...
       + k1*Bmp(1)*Chd(1) - k_1*BC(1) ...
       -q*lambda_tld_BC_space(1)*BC(1)/(1+Szd(1)/kit+(Chd(1)+BC(1))/kmt)-p*lambda_bmp1a_BC*BC(1)/(1+Szd(1)/kia+(Chd(1)+BC(1))/kma)- decBC*BC(1);
       
dBN(1) = DBN/dx^2*(-2*BN(1)+2*BN(2)) ...
       + k2*Bmp(1)*Nog(1) - k_2*BN(1)- decBN*BN(1);
   
%----------------------Internal node points--------------------------    
for i=2:1:n-1
       
dBmp(i) = DB/dx^2*(Bmp(i-1)-2*Bmp(i)+Bmp(i+1)) ...
        - k1*Bmp(i)*Chd(i) + k_1*BC(i) ...
        - k2*Bmp(i)*Nog(i) + k_2*BN(i) ...
        - decB*Bmp(i) + etaB(i) + q*lambda_tld_BC_space(i)*BC(i)/(1+Szd(i)/kit+(Chd(i)+BC(i))/kmt)+p*lambda_bmp1a_BC*BC(i)/(1+Szd(i)/kia+(Chd(i)+BC(i))/kma);
    
dChd(i) = DC/dx^2*(Chd(i-1)-2*Chd(i)+Chd(i+1)) ...
        - k1*Bmp(i)*Chd(i) + k_1*BC(i) ...
        - decC*Chd(i) + etaC(i) - q*lambda_tld_C_space(i)*Chd(i)/(1+Szd(i)/kit+(Chd(i)+BC(i))/kmt) - p*lambda_bmp1a_C*Chd(i)/(1+Szd(i)/kia+(Chd(i)+BC(i))/kma);
      
      
dNog(i) = DN/dx^2*(Nog(i-1)-2*Nog(i)+Nog(i+1)) ...
        - k2*Bmp(i)*Nog(i) + k_2*BN(i) ...
        - decN*Nog(i) + etaN(i);
        
dSzd(i) = DS/dx^2*(Szd(i-1)-2*Szd(i)+Szd(i+1))...
        - decS*Szd(i) + Vs*Bmp(i)^nu/(k^nu+Bmp(i)^nu);
      
dBC(i) = DBC/dx^2*(BC(i-1)-2*BC(i)+BC(i+1)) ...
       + k1*Bmp(i)*Chd(i) - k_1*BC(i) ...
       -q*lambda_tld_BC_space(i)*BC(i)/(1+Szd(i)/kit+(Chd(i)+BC(i))/kmt)-p*lambda_bmp1a_BC*BC(i)/(1+Szd(i)/kia+(Chd(i)+BC(i))/kma)- decBC*BC(i);
       
dBN(i) = DBN/dx^2*(BN(i-1)-2*BN(i)+BN(i+1)) ...
       + k2*Bmp(i)*Nog(i) - k_2*BN(i)- decBN*BN(i);
       
end
%--------------nth node point corresponds to dorsal midline----------------  
dBmp(n) = DB/dx^2*(2*Bmp(n-1)-2*Bmp(n)) ...
        - k1*Bmp(n)*Chd(n) + k_1*BC(n) ...
        - k2*Bmp(n)*Nog(n) + k_2*BN(n) ...
        - decB*Bmp(n) + etaB(n) + q*lambda_tld_BC_space(n)*BC(n)/(1+Szd(n)/kit+(Chd(n)+BC(n))/kmt)+p*lambda_bmp1a_BC*BC(n)/(1+Szd(n)/kia+(Chd(n)+BC(n))/kma);
    
dChd(n) = DC/dx^2*(2*Chd(n-1)-2*Chd(n)) ...
        - k1*Bmp(n)*Chd(n) + k_1*BC(n) ...
        - decC*Chd(n) + etaC(n) - q*lambda_tld_C_space(n)*Chd(n)/(1+Szd(n)/kit+(Chd(n)+BC(n))/kmt) - p*lambda_bmp1a_C*Chd(n)/(1+Szd(n)/kia+(Chd(n)+BC(n))/kma);
      
dNog(n) = DN/dx^2*(2*Nog(n-1)-2*Nog(n)) ...
        - k2*Bmp(n)*Nog(n) + k_2*BN(n) ...
        - decN*Nog(n) + etaN(n);
        
dSzd(n) = DS/dx^2*(2*Szd(n-1)-2*Szd(n))...
        - decS*Szd(n) + Vs*Bmp(n)^nu/(k^nu+Bmp(n)^nu);
      
dBC(n) = DBC/dx^2*(2*BC(n-1)-2*BC(n)) ...
       + k1*Bmp(n)*Chd(n) - k_1*BC(n) ...
       -q*lambda_tld_BC_space(n)*BC(n)/(1+Szd(n)/kit+(Chd(n)+BC(n))/kmt)-p*lambda_bmp1a_BC*BC(n)/(1+Szd(n)/kia+(Chd(n)+BC(n))/kma)- decBC*BC(n);
       
dBN(n) = DBN/dx^2*(2*BN(n-1)-2*BN(n)) ...
       + k2*Bmp(n)*Nog(n) - k_2*BN(n)- decBN*BN(n);
   
%-----------------------Update solution vector---------------------------
dY = [dBmp dChd dNog dSzd dBC dBN]';

end