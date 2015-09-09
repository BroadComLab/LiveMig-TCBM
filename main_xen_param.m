%
% MAIN-XEN-LIVEMIG 
% Code for Live Migration Virtual Machines
% Questo codice calcola i parametri delle macchine XEN
%
% author: Danilo Amendola (danilo.amendola@uniroma1.it)
% date: 26/01/2015      at Sapienza University of Rome
% 
% last update: 26/01/2015
% version: 0.1 beta
% -------------------------------------------------------------------
% -------------------------------------------------------------------

close all
clear all
clc

% ----- COSTANTS


% -------------------------------------------------------------------
% ----- VARIABLES --------------------------------------------------
% Final round index
Imax = 29; % 

% Iteration index []
N0 = 50;
% dirty bit rate Mbit/sec
vWavg = 200:0.5:450;
Wavg = vWavg(1);
% max utilizzation of the overall available basndwidth  Cmax  [0,1]
Cmax = 1;%0^(-2); 

% migration trnsmission rate (bandwidth) bit/s
R_i = zeros(N0+1, 1);
% Max transmission bandwidth allowed 
Rmax = 300;%10^3;%*10^6; % Mb/s
% overall bandwidth available at migrating VM bit/sec
Rtot = 2*Rmax;

Rcap=min(Rmax, Cmax*Rtot);

% Power wasted for Ri
K0 = 1.8*10^(-3); % Watt*[s/Mbit]^alfa communication pw
alfa_Pw = 1.45; % exponent power-rate
Pw_i = K0.*R_i.^alfa_Pw; % ATT sistemare!!!! 

% size of memory [Mega bit]
M0 = 2064; %%*10^6;
% binary var for pre-copy > teta = 1 or stop-and-copy > teta = 0
teta = 1;
% Time required by i-th migration round [t]
T_i = ones(Imax, 1)*1;
% energy wasted by connection setim phase Joule
Esetup = 3*10^(-4); % per il wifi
% Number of bits to be midgtrated at i-ty round bit
V_i = zeros(Imax, 1); 
% max number of bit to be migratedper round
Vmax = 2000; 
%communication energy wasted during  the i-th round
En_i = zeros(Imax, 1);
% maximum tolerated total migration time sec
%TimeDT = 0; %TimeDT = M0*Wavg_i(1)^(Imax);

%DeltaTM = 0.17*10^(1); % s 
% maximum tolerated down-time
%DeltaDT = 1.95*10^(-1); % s 1.96
% Data Compression ratio  0< ComprRatio <=1
ComprRatio = 0.5;
% Data protection factor  red >= 1
red = 1;

% Ratio V_i / V_i+1 of data migrated at rounds i and i+1 Beta_i>=1
%beta_i = ones(Imax, 1);
beta = 1.1; %1.2; %beta_i(1);

%beta_ave = 0;

gamma = 1000; % buono tra 100 e 1000
amax = 5*10^(-7);% buono tra 5*10(-7) e 10^(-6)

% ------
%Etot = 0;

DeltaR = (Rmax - Wavg) / (Imax + 1);

TimeDT = zeros(length(Wavg), 1);

% ------- END VAR --------------------------------------------------------


% ------ MAIN ------------------------------------------------------------
for w=1:length(vWavg)
    Wavg = vWavg(w);
    DeltaR = (Rmax - Wavg) / (Imax + 1);

    prodTimeDT = 1;
    for i=0:Imax+1
        prodTimeDT = prodTimeDT * ( 1 / (Wavg + i * DeltaR)); 
    end
    TimeDT(w) = M0 * ((Wavg)^(Imax+1)) * (prodTimeDT);

    sumTimeTM = 0;
    for i=0:(Imax+1)
        prodTimeTM = 1;
        for j=0:i
            prodTimeTM = prodTimeTM * (1 / (Wavg + (j * DeltaR))); 
        end
        sumTimeTM = sumTimeTM + (Wavg^i) * prodTimeTM;
    end
    TimeTM = M0*(sumTimeTM);


    beta_ave = 1 + (1/2) * (Imax / (Imax+1)) * ((Rmax/Wavg)-1);


    sumEtot = 0;
    for i=1:(Imax+1)
        prodEtot = 1;
        for j=0:(i-1)
            prodEtot = prodEtot * (1 / (Wavg + j * DeltaR)); 
        end
        sumEtot = sumEtot + K0*M0* (Wavg^i) * (Wavg + i * DeltaR)^(alfa_Pw-1) * prodEtot;
    end

    Etot = K0 * M0 * (Wavg^(alfa_Pw-1)) + sumEtot + Esetup;
end
% ------ Visualizza variabili risultato

TimeDT
TimeTM
Etot
beta_ave

%ImaxTilde = ceil((log(M0/(Rcap*TimeDT))/log(Rcap/Wavg)) - 1)
