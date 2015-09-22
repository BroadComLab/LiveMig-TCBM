% LiveMig Main Dynamic rate Code V4 
% for Live Migration Dynamic Virtual Machines 
% +  Tunable Complexy Bandwidth Manager (TCBM)
%
% author: Danilo Amendola (danilo.amendola@uniroma1.it)
% start date: 08/01/2015      at Sapienza University of Rome
% 
% update: 11/05/2015 sitemati bug
% last update: 08/09/2015 aggingo tou update v3 manoscritto pag.144
%
% version: 2.0e 
% -------------------------------------------------------------------
% -------------------------------------------------------------------

close all
clear all
clc

% ----- COSTANTS
TOU_UPDATE  =   3;  % { 0 = koushner  --> verificato
                    %   1 = our vers 2
                    %   2 = fast update k/n 
                    %   3 = update v3 da multipath-tpc }
PLOT_FIG    =   1;

ENE_SELECTNEW = 1;  % new energy = 1

% -------------------------------------------------------------------
% ----- VARIABLES --------------------------------------------------
Imax        =   3;              % Final round index
N0          =   5*10^3;         % Iteration index [] 
Wavg        =   3;              % dirty bit rate Mbit/sec
Cmax        =   1;              % max utilizzation of the overall available basndwidth  Cmax  [0,1]
Rmax        =   30;             % Mb/s    % Max transmission bandwidth allowed 
K0          =   1.2*10^-3;  % 1.8*10^(-3) Watt*[s/Mbit]^alfa communication pw % Power wasted for Ri 
alfa_Pw     =   2;              % exponent power-rate
M0          =   256;            % size of memory [Mega bit] 
Esetup      =   3*10^-4;%; per il wifi %% energy wasted by connection setim phase Joule 
gamma       =   1000;           % buono tra 100 e 1000
amax        =   0.05;      %2*10^(-4);%5*10^(-7);% buono tra 5*10(-7) e 10^(-6)
passoTou    =   100;
DeltaTM     =   117;%2*10^200;%141;%2*10^300;%94.6;%21.33;%2*10^600;2*10^300;       %0.17*10^(1); % maximum tolerated Time Migration
DeltaDT     =   0.06;%10^200;%10^300;%0.1; %1.17*10^-3;%19.31;%10^600;         %11.63636364; % maximum tolerated down-time
beta        =   0;            %1.05; %1.2; %beta_i(1);

% ris r0 = 11 lamda 1= 0 lamda2 1689.6 

% binary var for pre-copy > teta = 1 or stop-and-copy > teta = 0
if(Imax == -1)
    teta = 0; 
    Q=1; S=0;
elseif(Imax==0)
    teta = 1; 
    Q=1; S=0;
else
    teta = 1;
    % ----- DATI TC-BM ######
    Q = Imax;   %number of round head 
    S = Imax/Q;
    % -------------  ########
    if(Q>Imax || Q<1 || rem(Q,1)~=0 || (rem(S,1)~=0 && Imax>=1))
       error(['Errore S deve essere intero.']); 
    end
end

maxExpVal   =   10^50;         % Limita gli esponenziali nei gradienti dei lamda.
C           =    ones(Imax+2,1); % max utilizzation of the overall available basndwidth  Cmax  [0,1]
for m=2:Imax+2
    C(m)=Wavg^(m-1);
end

gammaV = ones(1, Imax+2);
gammaV(1)=gamma;
for m=2:Imax+1
    gammaV(m) = gamma*1; %% inibito
end
gammaV(Imax+2)=gamma;

amaxV = ones(1, Imax+2);
amaxV(1)=amax;
for m=2:Imax+2
    amaxV(m) = amax*1;%
end
%amaxV(Imax+2)=amax;

% overall bandwidth available at migrating VM bit/sec
Rtot    =   2*Rmax;
Rcap    =   min(Rmax, Cmax*Rtot);
R_i     =   zeros(N0, Imax+2); % migration trnsmission rate (bandwidth) bit/s

%communication energy wasted during  the i-th round
En_i = zeros(N0, 1); % Joule
feasibility = 1; % init ipothesis to feasible

% ------- END VAR --------------------------------------------------------
% ------------------------------------------------------------------------

% --- INIT Simulation ---------------------------------------------------
% varibili aggiunte
gradRtilde  = zeros(N0, Imax+2);
gradL1      = zeros(N0, 1);
gradL2      = zeros(N0, 1);
gradL3      = zeros(N0, Imax+1);

% Tou
omega_i     = zeros(N0,Imax+2);
tou_1       = zeros(N0,1);
tou_2       = zeros(N0,1);
Psi_3i      = zeros(N0,Imax+1);

%lamda = zeros(3, N0+1);
lamda_1 = zeros(N0+1, 1); 
lamda_2 = zeros(N0+1, 1);
if(Imax==-1)% INIZIALIZZO FITTIZIAMENTE A 2 
    lamda_3 = zeros(N0+1, 2);
elseif(Imax==0)
    lamda_3 = zeros(N0+1, 2);
else
    lamda_3 = zeros(N0+1, Imax+1);
end

%%%%%% ATT. DA RIMUOVERE
% lamda_3(1,1)= 76.8;
% lamda_3(1,2) = 768;
% %

Rtilde = zeros(N0, Imax+2);

if Imax==-1 || Imax==0
    Rtilde = zeros(N0, 2);
    R_i    = zeros(N0, 2);
end

Rtilde(1,:) = log(Rmax);%(Rmax+Wavg)/2);    %%%%%%%% START POINT %%%%%%%%%
%Rtilde(1,3) = log(10^-44);

% for i=0:Imax+1
%     Rtilde(1,i+1)=log(Wavg+i*((Rmax-Wavg)/(Imax+1)));
% end

% NEW INIT FOR TCBM
% Rtilde(1,1+0)=log(Wavg);
% for i=1:Q
%     Rtilde(1,((i-1)*S+2):((i)*S+1))=log(Wavg+i*((Rmax-Wavg)/(Q+1)));
% end
% Rtilde(1,1+Imax+1)=log(Rmax);

%%%%%%%% TOGLI
%     Rtilde(1,1) = log(5);
%     Rtilde(1,2) = log(5);
%     Rtilde(1,3) = -1000;
% %%%%



R_i(1, :)=exp(Rtilde(1,:)); % Rmax/1000000
    
if TOU_UPDATE == 3 % fast Tou
    omega_i(1,:) = amax;
    tou_1(1) = amax;
    tou_2(1) = amax;
    Psi_3i(1,:) = amax;
elseif TOU_UPDATE == 2 % fast Tou
    omega_i(1,:) = passoTou;
    tou_1(1) = passoTou;
    tou_2(1) = passoTou;
    Psi_3i(1,:) = passoTou;
elseif TOU_UPDATE == 1 % our versione 2
    omega_i(1,:) = amax;
    tou_1(1) = amax;
    tou_2(1) = amax;
    Psi_3i(1,:) = amax;
elseif TOU_UPDATE == 0 % Koushner
    omega_i(1,:) = amaxV;
    tou_1(1) = amax;
    tou_2(1) = amax;
    Psi_3i(1,:) = amax;
end

%% per Update Tou versione 2
D_amax = zeros(N0,1);
D_amax(1,:) = amax;

mV_i = ones(N0, Imax+2);    mV_i(1, :) = 0;
mA_1 = ones(N0, 1);         mA_1(1, :) = 0;
mA_2 = ones(N0, 1);         mA_2(1, :) = 0;
mB_i = ones(N0, Imax+1);    mB_i(1, :) = 0;

exitCondValue = zeros(N0+1, 1);

Etot = 0;

% -------------------------------------------------------------------
% ----- MAIN BODY

% ------- PROTOCOL ****************************************************
% ------- PHASE 1) Collect the input data (did it before)
% Wavg (b/s); Rmax (b/s); Cmax {0 or 1}; Rtot (b/s); Imax >= 1; M0 (bit); 
% teta {0 or 1}; Esetup (J); beta >= 1; K0 (W/(b/s)^alfa_Pw ); alfa_Pw>=2 ;
% DeltaTM (s) ; DeltaDT (s) ;  amax, gamma;

% ------- PHASE 2) Check input data and feasibility condition
%input_error = false;
%if (Cmax < 0 && Cmax > 1)|| Imax < 1 || (teta ~= 0 && teta ~=1) || beta < 1 || alfa_Pw < 2 
    %input_error = true;
    %error('Error: input value not allowed!');
%end

% ----- FEASIBILITY CHECH IN -------------------------------------
feasib1 = teta*( +...
    M0/DeltaTM*( +...
    1/Rcap*(Imax+2)*(delta_f(Wavg/Rcap-1))+ +...
    ((1-(Wavg/Rcap)^(Imax+2))/(Rcap-Wavg))* +...
    (1-delta_f(Wavg/Rcap-1))));
feasib2 = (M0/DeltaDT)*(1/Rcap)*((Wavg/Rcap)^(Imax+1));
% beta >= 1
feasib3 = teta*(beta*(Wavg/Rcap));

if (~(feasib1 < 1 && feasib2 < 1 && feasib3 < 1))
     error(['Error: feasibility check FAILED! ' 'feasib1: ' feasib1 ' feasib2: ' feasib2 ' feasib3: ' feasib3 ]);
     error = 'Error: feasibility check FAILED! '
     feasibility = 0;
end

% DeltaTM (s) ; DeltaDT (s) ;  amax, gamma;

% ------- PHASE 3) Implementation ordered set of iterates in n

% ------ BEGIN ROUND CYCLE ---------------------------------------------
for n = 1:N0
    % SE NON-FEASIBLE esci subito!
    if feasibility==0 
        break;
    end

% ---- #### NUOVA TCBM #### -----------------------
    bV = 1;
    % -------- TCBM Calcolo sommarorie e Ttm Tdt 
    sumTdt1 = 0;
    for k = 1:Q-1
        sumTdt1 = sumTdt1 + Rtilde(n, bV+k*S+1);
    end
    sumTmt1 = 0;
    for k = 0:Q-1
        % -- 
        for l = (k*S+1):((k+1)*S)
            sumTmt2 = 0;
            for p = 1:k-1
                sumTmt2 = sumTmt2 + Rtilde(n, bV+p*S+1);
            end
            sumTmt1 = sumTmt1+ (Wavg^l) * (delta_f(k) * min(maxExpVal ,exp(-l*Rtilde(n, bV+1))) + (1-delta_f(k))*min(maxExpVal ,exp((k*S-l)*Rtilde(n, bV+k*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(k-1))*S*(sumTmt2))));
        end
    end
    
    Tdt =  M0*min(maxExpVal ,exp(-Rtilde(n, bV+0)))* ...+
        ( delta_f(1+Imax)+(1-delta_f(1+Imax))*(Wavg^(1+Imax))*min(maxExpVal ,exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*(sumTdt1))) );% da init zeros(1,Imax+2);
    Ttm =  M0*min(maxExpVal ,exp(-Rtilde(n, bV+0)))*(1+(Wavg^(1+Imax)) * min(maxExpVal ,exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*(sumTdt1))) + ...+
        (1-delta_f(Imax)) * (sumTmt1));% da init zeros(1,Imax+2);
    
    % --- Fine TCBM Codice TM DT
    
    % ------ TCBM calcolo dei LAMDA e gradLambda
    
    gradL1temp = teta*((Ttm / DeltaTM)-1); 
    gradL2temp = ((Tdt/DeltaDT)-1);
    gradL3temp = zeros(1,Imax+1);
    
    gradL3temp(1,bV+0) = teta*(beta*Wavg*min(exp(-Rtilde(n,bV+0)), maxExpVal)-1);
    gradL3temp(1,bV+(1:S)) = teta*(beta*Wavg*min(exp(-Rtilde(n,bV+1)), maxExpVal)-1);
    
    for m=1:Q-1
        gradL3temp(1,bV+((m*S+1):((m+1)*S))) = teta*((beta*Wavg*min(exp(-Rtilde(n,bV+m*S+1)), maxExpVal))-1);%*round(1000*((beta*Wavg*exp(-Rtilde(n))-1)))/1000);
    end
        
    % LIMITIAMO i gradienti
    limiteGradienti = 10^5;%Inf;
    fatRound = 500;
%     gradRtilde(n,:) = max(min( round(gradRtildetemp.*fatRound)./fatRound , ones(1,Imax+2)*limiteGradienti), -ones(1,Imax+2)*limiteGradienti);
    gradL1(n) = max(min( round(gradL1temp.*fatRound)./fatRound , limiteGradienti), -limiteGradienti);
    gradL2(n) = max(min( round(gradL2temp.*fatRound)./fatRound , limiteGradienti), -limiteGradienti);
    gradL3(n,:) = max(min( round(gradL3temp.*fatRound)./fatRound , ones(1,Imax+1).*limiteGradienti), -ones(1,Imax+1).*limiteGradienti);

    % ------ CALCOLO DEI LAMDA
    lamda_1(n+1) = max((lamda_1(n) + tou_1(n).*gradL1(n)), 0);
    lamda_2(n+1) = max((lamda_2(n) + tou_2(n).*gradL2(n)), 0);
    if(teta~=0)
        lamda_3(n+1,:) = max((lamda_3(n,:) + Psi_3i(n,:).*gradL3(n,:)), zeros(1,Imax+1));
    end
    
    % -----------------------------------------------------------   
    % ----- GRAD Rtilde TCBM ----------- Tunable 
    % ------- CALCOLO LE SOMMATORIE DEI GRADIENTI per ith Rtilde
    % implementazione ad if innestati su richiesta del Prof.
    
    gradRtildetemp = zeros(1,Imax+2);
    derGrad = zeros(1,Imax+2);
    
    % -- R0  
    sumDerGrad1 = 0;
    for k = 1:Q-1
        sumDerGrad1 = sumDerGrad1 + Rtilde(n, bV+k*S+1);
    end
    sumDerGrad2 = 0;
    for m=0:Q-1
        sumInDerGrad1 = 0;
        for p=1:m-1
            sumInDerGrad1 = sumInDerGrad1+ Rtilde(n, bV+p*S+1);
        end
        for l=(m*S+1):((m+1)*S)
           sumDerGrad2 = sumDerGrad2 + Wavg^l*(delta_f(m)*exp((alfa_Pw-l)*Rtilde(n, bV+1))+(1-delta_f(m))*exp((alfa_Pw+m*S-l)*Rtilde(n, bV+m*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(m-1))*(sumInDerGrad1)) ); 
        end
    end
    derGrad_R0 = (alfa_Pw-1)*K0*M0*exp((alfa_Pw-1)*Rtilde(n, bV+0))- teta*K0*M0*exp(-Rtilde(n, bV+0)) * ...+
        (Wavg^(1+Imax)*exp((alfa_Pw-1)*Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)- (1-delta_f(Q-1))*S*(sumDerGrad1)) + ...+
        (1-delta_f(Imax))*( sumDerGrad2 )); 
    
    % update aggingiamo limitazioni del gradiente 
    limGrad = log(Rmax);
    gradRtildetemp(bV+0) = max(-limGrad, min(limGrad, derGrad_R0 - teta*lamda_1(n)*(1/DeltaTM)*Ttm - lamda_2(n)*(1/DeltaDT)*Tdt - teta*lamda_3(n, bV+ 0)*beta*Wavg*exp(-Rtilde(n, bV+0))));
    
if(Imax>=0)  %~Imax<0)
 if(Imax~=0)
    % ---------------- R1
    sumDerGradR1a = 0;
    sumDerGradR1b = 0;
    for  k=0:Q-1
        for l=(k*S+1):((k+1)*S)
            sumDerGradR1a = sumDerGradR1a + Wavg^l*(delta_f(k)*l*exp( -l*Rtilde(n, bV+1) ) + ...+
                (1-delta_f(k))*S*exp( (k*S-l)*Rtilde(n, bV+k*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(k-1))*S * sum_f(Rtilde(n,:), 1, k-1, S)) );
            sumDerGradR1b = sumDerGradR1b + Wavg^l*(((l-alfa_Pw)*delta_f(k))*exp((alfa_Pw-l)*Rtilde(n, bV+1) ) + ...+
                (1-delta_f(k))*S*exp( (alfa_Pw+k*S-l)*Rtilde(n, bV+k*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(k-1)) * sum_f(Rtilde(n, :), 1, k-1, S) ));
        end
    end
    
    derGrad_R1_Ttm = -M0*exp(-Rtilde(n, bV+0))*(S*(Wavg^(1+Imax))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n, :), 1, Q-1, S)) + ...+
        (1-delta_f(Imax))*( sumDerGradR1a ));    
    derGrad_R1_Tdt = -S*M0*exp(-Rtilde(n, bV+0))*Wavg^(1+Imax)*(1-delta_f(1+Imax))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n, :), 1, Q-1, S));
    derGrad_R1_en = -teta*K0*M0*exp(-Rtilde(n, bV+0))*(S*Wavg^(1+Imax)*exp((alfa_Pw-1)*Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*(sum_f(Rtilde(n,:) ,1,Q-1,S))) + ...+
        (1-delta_f(Imax))*(sumDerGradR1b));
    
    % AGGIORNAMENTO limitazione dei tre valori dellederivate precedenti
    %derGrad_R1_Ttm = derGrad_R1_Ttm
    % Update aggiungo limite gradienti
    gradRtildetemp(bV+(1:S)) =  max(-limGrad, min(limGrad,  derGrad_R1_en + lamda_1(n)*teta*(1/DeltaTM)*derGrad_R1_Ttm + lamda_2(n)*(1/DeltaDT)*derGrad_R1_Tdt - lamda_3(n, bV+ 1)*teta*beta*Wavg*exp(-Rtilde(n, bV+1))));
    
    % ------- Ri
    for j=1:Q-1 % ???
        sumDerGradRi_a = 0;
        sumDerGradRi_b = 0;
        for  k=j:Q-1
            for l=(k*S+1):((k+1)*S)
                sumDerGradRi_a = sumDerGradRi_a + Wavg^l*(delta_f(k-j)*(j*S-l)*exp( (j*S-l)*Rtilde(n, bV+j*S+1) -S*Rtilde(n, bV+1)-(1-delta_f(k-1))*S*(sum_f(Rtilde(n,:),1,j-1,S))) - ...+
                    (1-delta_f(k-j))*S*(1-delta_f(k-1))*exp( (k*S-l)*Rtilde(n, bV+k*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(k-1))*S * sum_f(Rtilde(n, :), 1, k-1, S) ) );

                sumDerGradRi_b = sumDerGradRi_b + Wavg^l*(1-delta_f(k))*(delta_f(k-j)*(alfa_Pw+j*S-l)*exp((alfa_Pw+j*S-l)*Rtilde(n, bV+j*S+1) -S*Rtilde(n, bV+1)-(1-delta_f(k-1))*sum_f(Rtilde(n, :),1,j-1,S)) - ...+
                    (1-delta_f(k-j))*(1-delta_f(k-1))*exp( (alfa_Pw+k*S-l)*Rtilde(n, bV+k*S+1)-S*Rtilde(n, bV+1)-(1-delta_f(k-1)) * sum_f(Rtilde(n, :), 1, k-1, S) ));
            end
        end
        
        derGrad_Ri_Tdt = -S*(1-delta_f(Q-1))*M0*exp(-Rtilde(n, bV+0))*Wavg^(1+Imax)*(1-delta_f(1+Imax))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n, :), 1, Q-1, S));
        derGrad_Ri_Ttm = M0*exp(-Rtilde(n, bV+0))*((-Wavg^(1+Imax))*S*(1-delta_f(Q-1))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n, :), 1, Q-1, S)) + ...+
            (1-delta_f(Imax))*( sumDerGradRi_a ));
        derGrad_Ri_en = teta*K0*M0*exp(-Rtilde(n, bV+0))*(-(1-delta_f(Q-1))*(Wavg^(1+Imax))*S*exp((alfa_Pw-1)*Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*(sum_f(Rtilde(n, :),1,Q-1,S))) + ...+
            (1-delta_f(Imax))*(sumDerGradRi_b));
        % update aggiunto limite gradiente 
        %GRADIENTE j
        gradRtildetemp(bV+((j*S+1):((j+1)*S))) =  max(-limGrad, min(limGrad,  derGrad_Ri_en + lamda_1(n)*(1/DeltaTM) *teta*derGrad_Ri_Ttm + lamda_2(n)*(1/DeltaDT) * derGrad_Ri_Tdt - lamda_3(n, bV+j*S+1)*teta*beta*Wavg*exp(-Rtilde(n, bV+j*S+1)) ));
    end
 end% Imax~=0
 
    % -------- R Imax+1
    derGrad_RImax_Ttm = -M0*(Wavg^(1+Imax))*exp(-Rtilde(n, bV+0))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n,:), 1, Q-1, S));
    derGrad_RImax_Tdt = -M0*exp(-Rtilde(n, bV+0))*Wavg^(1+Imax)*(1-delta_f(1+Imax))*exp(-Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*sum_f(Rtilde(n,:), 1, Q-1, S));
    derGrad_RImax_en = teta*K0*M0*Wavg^(1+Imax)*exp(-Rtilde(n, bV+0))*(alfa_Pw-1)*exp((alfa_Pw-1)*Rtilde(n, bV+Imax+1)-S*Rtilde(n, bV+1)-(1-delta_f(Q-1))*S*(sum_f(Rtilde(n,:),1,Q-1,S)));
     % update aggiunto limite gradiente
     % Imax+1
    gradRtildetemp(bV+Imax+1) =  max(-limGrad, min(limGrad, derGrad_RImax_en + lamda_1(n)*teta*(1/DeltaTM)*derGrad_RImax_Ttm + lamda_2(n)*(1/DeltaDT) * derGrad_RImax_Tdt));
end % Imax>0
    % ----------------------------------------
    % ------------- END NUOVI TCBM GRAD------
    
    

    % LIMITIAMO i gradienti
%     limiteGradienti = 10^5;%Inf;
%     fatRound = 100;
    gradRtilde(n,:) = max(min( round(gradRtildetemp.*fatRound)./fatRound , ones(1,Imax+2)*limiteGradienti), -ones(1,Imax+2)*limiteGradienti);

    % --- STEP 4 TCBM
    % CALCOLO valore Rtilde i+1th     %%% ATTENZIONE introdotto max(... ,
    % -10^2*ones...) per problema di exp(-...R...) che va infinito
    Rtilde(n+1,:) = max(min((Rtilde(n,:) - omega_i(n,:).*gradRtilde(n,:)),(ones(1,Imax+2).*log(Rcap))), (-10^2)*ones(1,Imax+2));
    
    
    
    %(K0, M0, alfa_Pw, Rtilde, teta, Wavg_i(1), lamda(1,i), lamda(2,i), lamda(3,i), DeltaTM, DeltaDT, Imax, beta))^(log(Rcap));     
    
    % TOU UPDATE, possibili tre differenti aggiornamenti del motiplicatore 
    if TOU_UPDATE == 3 % nuovo aggiornamento dei moltiplicatori pag 144 sett-2015
               
        % ------ AGIORNAMENTO DEI TOU sett-2015 (V3)
        amaxMin = amax/10;
        
        omega_i(n+1, :) = max( ones(1,Imax+2).*amaxMin, min(ones(1,Imax+2).*amax, (1/2).*(Rtilde(n+1,:)).^2)); % 
        tou_1(n+1, :) = max( amaxMin, min(amax, (1/2)*(lamda_1(n+1))^2)); % 
        tou_2(n+1, :) = max( amaxMin, min(amax, (1/2)*(lamda_2(n+1))^2)); % 
        Psi_3i(n+1, :) = max( ones(1,Imax+1).*amaxMin, min(ones(1,Imax+1).*amax, (1/2).*(lamda_3(n+1, :)).^2)); % 
        
        
        
    elseif TOU_UPDATE == 2 % 1) FAST UPDATE
        omega_i(n+1,:) = passoTou/n;
        tou_1(n+1,:) = passoTou/n;
        tou_2(n+1,:) = passoTou/n;
        Psi_3i(n+1,:) = passoTou/n;
    elseif TOU_UPDATE == 1 % 2) OUR v2 version  
        %for i=1:Imax+1
        D_amax(n+1) = amax * abs(sqrt(((Rtilde(n+1,:) - Rtilde(n,:))*(Rtilde(n+1,:) - Rtilde(n,:))') + ((lamda_1(n+1) - lamda_1(n))^2)+(lamda_2(n+1) - lamda_2(n))^2 +((lamda_3(n+1,:)-lamda_3(n,:))*(lamda_3(n+1,:)-lamda_3(n,:))'))); %%% E PER LAMDA3?? 
        %end
        
         % ------ AGIORNAMENTO DEI TOU AMAX SEMPLIFICATO (V2)
        amaxMin = amax/10;
        omega_i(n+1, :) = max( amaxMin, min(amax, D_amax(n+1))); %  passoTou/n; %
        tou_1(n+1, :) = max( amaxMin, min(amax, D_amax(n+1))); %  passoTou/n; %
        tou_2(n+1, :) = max( amaxMin, min(amax, D_amax(n+1))); %  passoTou/n; %
        Psi_3i(n+1, :) = max( amaxMin, min(amax, D_amax(n+1))); %  passoTou/n; %
    
    elseif TOU_UPDATE == 0 % 3) Koushner VERSION
        moltAmax = 1;%7*10^1;
        minValue = amax/20;

        % ------ AGIORNAMENTO DEI TOU
        % AGGIUNTA MIN-VALUE amax/10 
        omega_i(n+1,:) = max(ones(1,Imax+2).*minValue, min(ones(1,Imax+2).*amaxV, (omega_i(n, :) - gammaV.*mV_i(n,:).*gradRtilde(n,:)))); %  passoTou/n; %
        tou_1(n+1) = max(1*minValue, min(amax*moltAmax, tou_1(n) + gamma*mA_1(n).*gradL1(n)));
        tou_2(n+1) = max(1*minValue, min(amax*moltAmax, tou_2(n) + gamma*mA_2(n).*gradL2(n)));
        Psi_3i(n+1,:) = max(ones(1,Imax+1).*minValue, min(ones(1,Imax+1).*amax, (Psi_3i(n, :) + gamma.*mB_i(n,:).*gradL3(n,:)))); % tolto 20*
        % ------- Aggiornamento dei bit trasmessi V
        
        mV_i(n+1, :) = (1 - omega_i(n, :)).* mV_i(n,:) - gradRtilde(n,:);
        mA_1(n+1) = (1 - tou_1(n)).* mA_1(n) + gradL1(n);
        mA_2(n+1) = (1 - tou_2(n)).* mA_2(n) + gradL2(n);
        mB_i(n+1, :) = (1 - Psi_3i(n, :)).* mB_i(n,:) + gradL3(n,:);
        
    end % end TOU_UPDATE

    % ------ STEP 5 TCBM  compute Rate
    R_i(n+1,:) = exp(Rtilde(n+1,:));  
    
    % ------------------------------------------
    % ------ CALCOLO ENERGIA all'iesimo passo
if ENE_SELECTNEW == 0
    sumEi = 0;
    for i=2:Imax+2
        %sumEi = sumEi + K0*M0*C(i)*(R_i(n,i)^(alfa_Pw-1))*prod((R_i(n, 1:i-1).^(-1)));% ATTENZIONE CAMBIATO min(R_i(n, 1:i-1).^(-1), 10^-10))
        sumEi = sumEi + K0*M0*C(i)*(R_i(n,i)^(alfa_Pw-1))*prod(R_i(n, 1:i-1).^(-1));%max(R_i(n, 1:i-1).^(-1), 10^-20));%
    end
    En_i(n) = (K0*M0*(R_i(n,1))^(alfa_Pw-1) + teta * (sumEi) + Esetup);
end

% ----- ENRGIA NUOVA FORMULA TCBM -----------------------
if ENE_SELECTNEW == 1 
    sumEneTC2 = 0;
    prodEneTC1 = 1;
    for m=0:Q-1
        prodEneTC1 = prodEneTC1 * (R_i(n, bV+m*S+1)^(-S));
        prod2EneTC2 = 1;
        for p=0:m-1
           prod2EneTC2 = prod2EneTC2 * (R_i(n, bV+(p*S+1))^(-S));
        end
        for l=(m*S+1):((m+1)*S)
            
            sumEneTC2 = sumEneTC2 + (Wavg^(l)) * (delta_f(m)*(R_i(n, bV+1)^(alfa_Pw-l)) + (1-delta_f(m))*(prod2EneTC2)*(R_i(n, bV+m*S+1)^(alfa_Pw+m*S-l))) ;
        end
    end
    En_i(n) = Esetup + K0*M0*(R_i(n, bV+0)^(alfa_Pw-1)) + teta * K0*M0*(R_i(n, bV+0)^(-1)) * ...+
        ((Wavg^(1+Imax))*(R_i(n,bV+Imax+1)^(alfa_Pw-1))*(prodEneTC1) + (1-delta_f(Imax))*sumEneTC2 );
end
 % ----- FINE ENRGIA NUOVA FORMULA TCBM -----------------------
 
    % boolean exitCondition value true = exit
%    exitCondValue(n+1) = ((R_i(n+1)-R_i(n))^2)/(R_i(n+1)^2);
%    exitCondition = exitCondValue(n+1) <= 10^(-2);
    
    %---------------- CONTROLLO SUI GRADIENTI 
%    if n > N0/10 && abs(gradRtilde(n))<10^(-5) && abs(gradL1(n))<10^(-5) && abs(gradL2(n))<10^(-5) && abs(gradL3(n))<10^(-5)
%        break;
%    elseif n > (N0-10) && ~ (abs(gradRtilde(n))<10^(-5) && abs(gradL1(n))<10^(-5) && abs(gradL2(n))<10^(-5) && abs(gradL3(n))<10^(-5))
%        N0 = 2*N0;
%    end


% % ------- CHECH USCITA ANTICIPATA DAL CICLO.
%     if n > 2*10^3 % cicli minimi da relizzare
%         % -----------------------------------------------------------------
%         % ----- FEASIBILITY OUT CHECH -------------------------------------
%         % IMPLEMENTARE VINCOLI DI CONTROLLO SULL Ropt per capire se sono
%         % soddisfatti. con il valore finale.
%         Ropt = R_i(n+1);
%         feasibOUT1 = teta*( +...
%             M0/DeltaTM*( +...
%             1/Ropt*(Imax+2)*(delta_f(Wavg/Ropt-1))+ +...
%             ((1-(Wavg/Ropt)^(Imax+2))/(Ropt-Wavg))* +...
%             (1-delta_f(Wavg/Ropt-1))));
%         feasibOUT2 = (M0/DeltaDT)*(1/Ropt)*(Wavg/Ropt)^(Imax+1);
%         % beta >= 1
%         feasibOUT3 = teta*(beta*(Wavg/Ropt));
% 
%         % -------- SE abbiamo feasibility in OUT ed IN, ed abbiamo che:
%         % -------- GRADIENTI converge alla K cifra decimale ESCO
%         %BoolGradStop = true;
%         offsetGR= 1000;
% 
%         DeltaGrad = [gradRtilde(n-offsetGR)-gradRtilde(n), gradL1(n-offsetGR)-gradL1(n), gradL2(n-offsetGR)-gradL2(n), gradL3(n-offsetGR)-gradL3(n)];
%         BoolGradStop = false;
% 
%         % Se la varioazione dei gradienti rispetto a OFFSET rount sono
%         % minori di valore
%         if DeltaGrad <= 10^(-4) 
%             % E se o vaòpro dei gradienti sono molto piccoli
%             if abs(gradRtilde(n))<10^(-5) && abs(gradL1(n))<10^(-5) && abs(gradL2(n))<10^(-5) && abs(gradL3(n))<10^(-5)
%             %if gradRtilde(n)==0 && gradL1(n)==0 && gradL2(n)==0 && gradL3(n)==0
%                 BoolGradStop = true;
%             end
%         end
% 
%         if  BoolGradStop && (feasibOUT1 < 1 && feasibOUT2 < 1 && feasibOUT3 < 1) && feasibility==1
%             break;
%         end
%     end
% 
%     % ----------------------------  

    % ------- STEP 5) Replace opt value
%     if n==N0 % exitCondition
%         % break;
%     end

end %for 1-Imax

% -------- CALCOLO ENERGIA ALL'USCITA del ciclo N0
Ropt = R_i(N0,:);
En_opt = En_i(N0);

% -------- CALCOLO TEMPI E DATI OPERATIVI Con Ropt --------

% Time required by i-th migration round [t]
Tm_i = ones(Imax+1, 1);
% Number of bits to be midgtrated at i-ty round bit
Vm_i = zeros(Imax+1, 1); 
% Ratio V_i / V_i+1 of data migrated at rounds i and i+1 Beta_i>=1
Vm_i(1) = M0;
Tm_i(1) = M0/Ropt(1);
for m=2:Imax+1
    Vm_i(m) = Wavg*Tm_i(m-1);
    Tm_i(m) = Vm_i(m)/Ropt(m);
end
thDeltaDT = ones(Imax+1, 1)*DeltaDT;

figure(6)
%plot(V_i)
hold on
plot(Tm_i)
plot(thDeltaDT)
title('TEMPI DI ROUND')

% -------- Calcolo di I tilde ---------------------------
ImaxTilde = ceil((log(M0/(Rcap*DeltaDT))/log(Rcap/Wavg)) - 1)


% -----------------------------------------------------------------
% ----- FEASIBILITY OUT CHECH -------------------------------------
% IMPLEMENTARE VINCOLI DI CONTROLLO SULL Ropt per capire se sono
% soddisfatti. con il valore finale.
% if feasibility==1    
%     Ropt = R_i(n);
%     feasibOUT1 = teta*( +...
%         M0/DeltaTM*( +...
%         1/Ropt(1)*(Imax+2)*(delta_f(Wavg/Ropt(1)-1))+ +...
%         ((1-(Wavg/Ropt(1))^(Imax+2))/(Ropt(1)-Wavg))* +...
%         (1-delta_f(Wavg/Ropt(1)-1))));
% 
%     feasibOUT2 = (M0/DeltaDT)*(1/Ropt(Imax+1))*(Wavg/Ropt(Imax+1))^(Imax+1);
% 
%     beta >= 1
%     feasibOUT3 = teta*(beta*(Wavg/Ropt(Imax+1)));
% 
%     if (~(feasibOUT1 < 1 && feasibOUT2 < 1 && feasibOUT3 < 1))
%          error('Error: feasibility OUT check FAILED! Aumenta N0!')
%          error = 'Error: feasibility OUT check FAILED! Aumenta N0!'
% 
%          feasibility = -1;
%     end
% end


Etot

% vG = [-1, -1, -1, -1];
% if feasibility==1
%     vG = [gradRtilde(n-2)-gradRtilde(n-1), gradL1(n-2)-gradL1(n-1), gradL2(n-2)-gradL2(n-1), gradL3(n-2)-gradL3(n-1)];
% end


% ------------ PLOT FIGURE ----------------------------------
if PLOT_FIG == 1
    %Figura di controllo per la convergenza
    figure(1)
    plot(R_i)
    %hold on
    title('RATE')
    legend('rate R_i')
    
    figure(2)
    plot(lamda_1)
    hold on 
    plot(lamda_2)
    plot(lamda_3)
    legend('lamda1', 'lamda2', 'lamda3');
    title('Plot dei lamda');
    
    figure(4)
    plot(gradL1)
    hold on
    plot(gradL2)
    plot(gradL3)
    legend('Grad1', 'Grad2', 'Grad3');
    title('Gradienti')
    
    figure(5)
    plot(En_i)
    title('ENERGIA')
    legend('Energia')
end
% ------------ FINE PLOT FIGURE ----------------------------------

%% FEASIBILITY CHECK OUT 

[feasibility_res, TM, TD] = feasibility_check(R_i(N0,:), M0, Wavg, DeltaTM, DeltaDT)

if feasibility_res == 1
    error('Error: downtime check FAILED!');
elseif feasibility_res == 2
    error('Error: Total migration time exceeded');
end
