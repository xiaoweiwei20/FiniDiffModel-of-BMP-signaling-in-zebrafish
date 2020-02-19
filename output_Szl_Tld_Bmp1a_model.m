% Code name: parameter_screening_grid_Nogsink_kinetic.m
% Author: Xu Wang
% Last update: 02/19/2020
%
% Model outputs for each parameter combination will be compared againest 
% experimental results at 5.7 hpf for WT,CLF,NLF,ALF,TLF,SLF and TALF. The total
% error sum of squares for each of these combinations will be calculated.
% Later we will  set up a SSE threshold based on the initial values above,
% and figure out how many combinations (and which combinations) have SSE's
% <= the set threshold and have similar profiles to experiments.
% 
clear; clc;
% 
% %-------------- preparation before screening ----------------
% % %%%%%%%%%% load sparse parameter grid for screening
load para_grid;
load para_grid_jBC;
load para_grid_k;
load para_grid_ki;

%%%%%%%%%% Constant model parameters for model
  %%% nodes and time range
n = 36;   % number of nodes to evaluate finite difference equations
tRange = [0 7900];  % time interval to evaluate differential equations 5.7hpf
  %%% geometry
parameters.Ltot = 700;    % length of the embryo (1D assumption)       microns
parameters.Lven = 400;    % length of ventral region for Bmp            microns
parameters.LvenTld = 700; % length of ventral region for Tld         microns
% parameters.LvenSzd = 230; % length of ventral region for Sizzled     microns
parameters.LdorC = 140;   % length of dorsal Chd expression region   microns
parameters.LdorN = 78;   % length of dorsal Nog expression region   microns
  %%% diffusion rates
parameters.DB = 4.4;       % diffusion rate of ligand (BMP)    (microns^2*s^-1)*60s*m^-1
% parameters.DC = 7;       % diffusion rate of Chd    (microns^2*s^-1)*60s*m^-1
parameters.DS = 10;       % diffusion rate of Sizzled   (microns^2*s^-1)*60s*m^-1
  %%% decay rates
parameters.decB = 8.9*10^(-5);       % decay rate of BMP
parameters.decC = 9.6*10^(-5);       % decay rate of Chd
parameters.nu = 4;  % cooperative parameter
parameters.Vs = 100;
%-------------- screening operations ----------------
MBMP_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MChd_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MNog_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MSzd_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MBC_WT   = zeros(n,1000000);  % initialization of a vecotor storing WT model outputs for all parameter combinations
MBN_WT   = zeros(n,1000000); % initialization of a vecotor storing WT model outputs for all parameter combinations

MBMP_Chet   = zeros(n,1000000);
MChd_Chet   = zeros(n,1000000);
MNog_Chet   = zeros(n,1000000);
MSzd_Chet   = zeros(n,1000000);
MBC_Chet    = zeros(n,1000000);
MBN_Chet    = zeros(n,1000000);

MBMP_CLF   = zeros(n,1000000);
MChd_CLF   = zeros(n,1000000);
MNog_CLF   = zeros(n,1000000);
MSzd_CLF   = zeros(n,1000000);
MBC_CLF    = zeros(n,1000000);
MBN_CLF    = zeros(n,1000000);

MBMP_NLF  = zeros(n,1000000);
MChd_NLF  = zeros(n,1000000);
MNog_NLF  = zeros(n,1000000);
MSzd_NLF  = zeros(n,1000000);
MBC_NLF   = zeros(n,1000000);
MBN_NLF   = zeros(n,1000000);

MBMP_Ahet  = zeros(n,1000000);
MChd_Ahet  = zeros(n,1000000);
MNog_Ahet  = zeros(n,1000000);
MSzd_Ahet  = zeros(n,1000000);
MBC_Ahet   = zeros(n,1000000);
MBN_Ahet   = zeros(n,1000000);

MBMP_ALF   = zeros(n,1000000);
MChd_ALF   = zeros(n,1000000);
MNog_ALF   = zeros(n,1000000);
MSzd_ALF   = zeros(n,1000000);
MBC_ALF    = zeros(n,1000000);
MBN_ALF    = zeros(n,1000000);

MBMP_TLF   = zeros(n,1000000);
MChd_TLF   = zeros(n,1000000);
MNog_TLF   = zeros(n,1000000);
MSzd_TLF   = zeros(n,1000000);
MBC_TLF    = zeros(n,1000000);
MBN_TLF    = zeros(n,1000000);

MBMP_TALF   = zeros(n,1000000);
MChd_TALF   = zeros(n,1000000);
MNog_TALF   = zeros(n,1000000);
MSzd_TALF   = zeros(n,1000000);
MBC_TALF    = zeros(n,1000000);
MBN_TALF    = zeros(n,1000000);

MBMP_SLF   = zeros(n,1000000);
MChd_SLF   = zeros(n,1000000);
MNog_SLF   = zeros(n,1000000);
MSzd_SLF   = zeros(n,1000000);
MBC_SLF    = zeros(n,1000000);
MBN_SLF    = zeros(n,1000000);
SSE = zeros(10,1000000);  % initialization of a vecotor storing final SSE's for all parameter combinations
j = 1;
for i=1:50000
disp(['i = ', int2str(i)]);
%%%%%%%%%% parameters to be screened
parameters.Vs = 100;
parameters.DN = para_grid(i,1);   % diffusion rate of Noggin          (microns^2*s^-1)*60s*m^-1
parameters.DBC = para_grid(i,2);  % diffusion rate of [BMP,Chd]       (microns^2*s^-1)*60s*m^-1
parameters.DBN = para_grid(i,3);   % diffusion rate of [BMP,Nog]       (microns^2*s^-1)*60s*m^-1
parameters.DC =para_grid(i,4);       % diffusion rate of Chd    (microns^2*s^-1)*60s*m^-1
parameters.decN = para_grid(i,5);   % decay rate of Nog             nM*m^-1
parameters.decS = para_grid(i,6); % decay rate of Sizzled            nM*m^-1
parameters.decBC = para_grid(i,7); % decay rate of BC
parameters.decBN = para_grid(i,8); % decay rate of BN
parameters.j3 = para_grid(i,9);  % production rate of Noggin       nM*m^-1
parameters.k1   = para_grid(i,10);   % binding rates for BMP ligand and Chordin          nM^-1*m^-1
parameters.k_1  = para_grid(i,10);   % unbinding rates for BMP ligand and Chordin        m^-1
parameters.k2   = para_grid(i,11);  % binding rates for BMP ligand and Noggin           nM^-1*m^-1
parameters.k_2  = 0.1*para_grid(i,11);   % unbinding rates for BMP ligand and Noggin         m^-1
parameters.kmt = para_grid(i,12);       % Michaelis constant of the proteinase Tolloid           nM
parameters.kma = para_grid(i,13);       % Michaelis constant of the proteinase bmp1a           nM
parameters.lambda_tld_C = para_grid(i,14);   % tld processing rate of Chd  nM^-1*m^-1
parameters.lambda_tld_BC = para_grid(i,15);    % tld processing rate of BC   nM^-1*m^-1
parameters.lambda_bmp1a_C = para_grid(i,16);    % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC = para_grid(i,17);  % bmp1a processing rate of BC   nM^-1*m^-1
parameters.j1 = para_grid_jBC(i,1);      % production rate of BMP          nM*m^-1 
parameters.j2 = para_grid_jBC(i,2);   % production rate of Chordin      nM*m^-1
parameters.kit = para_grid_ki(i,1);        % the inhibitor constant(dissociation constants) of the protianse tolloid  nM
parameters.kia= para_grid_ki(i,2);       % the inhibitor constant(dissociation constants) of the protianse bmp1a  nM
parameters.k= k(1,i);         % parameter for the hill function
%%%%%%%%% ODE solver
% solve for WT

[B1, C1, N1, S1,BC1, BN1] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B1);
B_WT = B1(m,:);                                     
C_WT = C1(m,:);
N_WT = N1(m,:);
S_WT = S1(m,:);
BC_WT = BC1(m,:);
BN_WT = BN1(m,:);

% %solve for chd hetes

parameters.j2 = 0.5 *para_grid_jBC(i,2);   % production rate of Chordin      nM*m^-1
[B2, C2, N2, S2, BC2, BN2] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B2);
B_Chet = B2(m,:);                                     
C_Chet = C2(m,:);
N_Chet = N2(m,:);
S_Chet = S2(m,:);
BC_Chet = BC2(m,:);
BN_Chet = BN2(m,:);

%solve for chd homo mutant

parameters.j2 = 0;   % production rate of Chordin      nM*m^-1
[B3, C3, N3, S3, BC3, BN3] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B3);
B_CLF = B3(m,:);                                     
C_CLF = C3(m,:);
N_CLF = N3(m,:);
S_CLF = S3(m,:);
BC_CLF = BC3(m,:);
BN_CLF = BN3(m,:);

% %solve for Nog MO

parameters.j2 = para_grid_jBC(i,2);  % production rate of Chordin      nM*m^-1
parameters.j3 = 0;   % production rate of Noggin       nM*m^-1
[B4, C4, N4, S4, BC4, BN4] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B4);
B_NLF = B4(m,:);                                     
C_NLF = C4(m,:);
N_NLF = N4(m,:);
S_NLF = S4(m,:);
BC_NLF = BC4(m,:);
BN_NLF = BN4(m,:);

% %solve for bmp1a hetes
parameters.j3 = para_grid(i,9);   % production rate of Noggin       nM*m^-1
parameters.lambda_bmp1a_C = 0.5* para_grid(i,16);    % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC =0.5* para_grid(i,17);  % bmp1a processing rate of BC   nM^-1*m^-1

[B5, C5, N5, S5, BC5, BN5] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B5);
B_Ahet= B5(m,:);                                     
C_Ahet = C5(m,:);
N_Ahet = N5(m,:);
S_Ahet = S5(m,:);
BC_Ahet = BC5(m,:);
BN_Ahet = BN5(m,:);
% 
% %solve for bmp1a homo mutant

parameters.lambda_bmp1a_C = 0;    % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC =0;  % bmp1a processing rate of BC   nM^-1*m^-1

[B6, C6, N6, S6, BC6, BN6] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B6);
B_ALF= B6(m,:);                                     
C_ALF = C6(m,:);
N_ALF = N6(m,:);
S_ALF = S6(m,:);
BC_ALF = BC6(m,:);
BN_ALF = BN6(m,:);
% 
% %solve for TLd homo mutant
% 
parameters.lambda_bmp1a_C = para_grid(i,16);    % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC = para_grid(i,17);  % bmp1a processing rate of BC   nM^-1*m^-1
parameters.lambda_tld_C = 0;   % tld processing rate of Chd  nM^-1*m^-1
parameters.lambda_tld_BC = 0;    % tld processing rate of BC   nM^-1*m^-1

[B7, C7, N7, S7, BC7, BN7] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B7);
B_TLF= B7(m,:);                                     
C_TLF = C7(m,:);
N_TLF = N7(m,:);
S_TLF = S7(m,:);
BC_TLF = BC7(m,:);
BN_TLF = BN7(m,:);
% 
% %solve for bmp1a homo mutant and Tld MO
% 
parameters.lambda_bmp1a_C = 0;    % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC = 0;  % bmp1a processing rate of BC   nM^-1*m^-1

[B8, C8, N8, S8, BC8, BN8] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B8);
B_TALF= B8(m,:);                                     
C_TALF = C8(m,:);
N_TALF = N8(m,:);
S_TALF = S8(m,:);
BC_TALF = BC8(m,:);
BN_TALF = BN8(m,:);
% 
% %solve for sizzled homo mutant
% 
parameters.lambda_tld_C = para_grid(i,14);    % tld processing rate of Chd  nM^-1*m^-1
parameters.lambda_tld_BC = para_grid(i,15);    % tld processing rate of BC   nM^-1*m^-1
parameters.lambda_bmp1a_C = para_grid(i,16);     % bmp1a processing rate of Chd  nM^-1*m^-1
parameters.lambda_bmp1a_BC = para_grid(i,17);   % bmp1a processing rate of BC   nM^-1*m^-1
parameters.Vs = 0;
[B9, C9, N9, S9, BC9, BN9] = FiniDiffModSzd_lambdaspace(n, tRange, parameters);
[m,~] = size(B9);
B_SLF= B9(m,:);                                     
C_SLF = C9(m,:);
N_SLF = N9(m,:);
S_SLF = S9(m,:);
BC_SLF = BC9(m,:);
BN_SLF = BN9(m,:);

%%%%%%%%%% model results scaling

MBMP_WT(:,j)  = B_WT(1:n);
MChd_WT(:,j)   = C_WT(1:n);
MNog_WT(:,j)   = N_WT(1:n);
MSzd_WT(:,j)   = S_WT(1:n);
MBC_WT(:,j)    = BC_WT(1:n);
MBN_WT(:,j)    = BN_WT(1:n);

MBMP_Chet(:,j)  = B_Chet(1:n);
MChd_Chet(:,j)   = C_Chet(1:n);
MNog_Chet(:,j)   = N_Chet(1:n);
MSzd_Chet(:,j)   = S_Chet(1:n);
MBC_Chet(:,j)    = BC_Chet(1:n);
MBN_Chet(:,j)    = BN_Chet(1:n);

MBMP_CLF(:,j)  = B_CLF(1:n);
MChd_CLF(:,j)   = C_CLF(1:n);
MNog_CLF(:,j)   = N_CLF(1:n);
MSzd_CLF(:,j)   = S_CLF(1:n);
MBC_CLF(:,j)    = BC_CLF(1:n);
MBN_CLF(:,j)    = BN_CLF(1:n);

MBMP_NLF(:,j)  = B_NLF(1:n);
MChd_NLF(:,j)   = C_NLF(1:n);
MNog_NLF(:,j)   = N_NLF(1:n);
MSzd_NLF(:,j)   = S_NLF(1:n);
MBC_NLF(:,j)    = BC_NLF(1:n);
MBN_NLF(:,j)    = BN_NLF(1:n);

MBMP_Ahet(:,j)  = B_Ahet(1:n);
MChd_Ahet(:,j)   = C_Ahet(1:n);
MNog_Ahet(:,j)   = N_Ahet(1:n);
MSzd_Ahet(:,j)   = S_Ahet(1:n);
MBC_Ahet(:,j)    = BC_Ahet(1:n);
MBN_Ahet(:,j)    = BN_Ahet(1:n);

MBMP_ALF(:,j)  = B_ALF(1:n);
MChd_ALF(:,j)   = C_ALF(1:n);
MNog_ALF(:,j)   = N_ALF(1:n);
MSzd_ALF(:,j)   = S_ALF(1:n);
MBC_ALF(:,j)    = BC_ALF(1:n);
MBN_ALF(:,j)    = BN_ALF(1:n);

MBMP_TLF(:,j)  = B_TLF(1:n);
MChd_TLF(:,j)   = C_TLF(1:n);
MNog_TLF(:,j)   = N_TLF(1:n);
MSzd_TLF(:,j)   = S_TLF(1:n);
MBC_TLF(:,j)    = BC_TLF(1:n);
MBN_TLF(:,j)    = BN_TLF(1:n);

MBMP_TALF(:,j)  = B_TALF(1:n);
MChd_TALF(:,j)   = C_TALF(1:n);
MNog_TALF(:,j)   = N_TALF(1:n);
MSzd_TALF(:,j)   = S_TALF(1:n);
MBC_TALF(:,j)    = BC_TALF(1:n);
MBN_TALF(:,j)    = BN_TALF(1:n);

MBMP_SLF(:,j)  = B_SLF(1:n); 
MChd_SLF(:,j)   = C_SLF(1:n);
MNog_SLF(:,j)   = N_SLF(1:n);
MSzd_SLF(:,j)   = S_SLF(1:n);
MBC_SLF(:,j)    = BC_SLF(1:n);
MBN_SLF(:,j)    = BN_SLF(1:n);
j=j+1;
end

%-------------- save data ----------------

save('modeldata_1to5w.mat','MBMP_WT','MChd_WT','MNog_WT','MSzd_WT','MBC_WT','MBN_WT',...
     'MBMP_Chet','MChd_Chet','MNog_Chet','MSzd_Chet','MBC_Chet','MBN_Chet',...
     'MBMP_CLF','MChd_CLF','MNog_CLF','MSzd_CLF','MBC_CLF','MBN_CLF',...
     'MBMP_NLF','MChd_NLF','MNog_NLF','MSzd_NLF','MBC_NLF','MBN_NLF',...
     'MBMP_Ahet','MChd_Ahet','MNog_Ahet','MSzd_Ahet','MBC_Ahet','MBN_Ahet',...
     'MBMP_ALF','MChd_ALF','MNog_ALF','MSzd_ALF','MBC_ALF','MBN_ALF',...
     'MBMP_TLF','MChd_TLF','MNog_TLF','MSzd_TLF','MBC_TLF','MBN_TLF',...
     'MBMP_TALF','MChd_TALF','MNog_TALF','MSzd_TALF','MBC_TALF','MBN_TALF',...
     'MBMP_SLF','MChd_SLF','MNog_SLF','MSzd_SLF','MBC_SLF','MBN_SLF');
