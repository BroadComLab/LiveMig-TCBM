% Funzione di controllo feasibility (in uscita)
% 
% Date 18/09/2015 
% Last update  current vers. 0.01

function [notfeasible, TTM, TTD] = feasibility_check(R, M0, Wavg, TM, DT)
% function [err] = controllo(R,M0,Wavg,TM,DT) return 1 if not feasible
% R is the vector of rates per round
% M0 is the vm memory value
% Wavg is the avg dirty rate
% TM and TD are the time contstraints

    notfeasible = 0;
    %R       = ones(Imax+2,1)*2;     % Vettore dei rate
    Imax = length(R)-2;
    Volume  = zeros(Imax+2,1);      % Vettore che contiene la dimensione dei  bit migrati ad ogni round
    Time    = zeros(Imax+2,1);      % Vettore che contiene la durata di ogni round 

    Volume(1)=M0;
    Time(1) = Volume(1)/R(1);

    for i=2:length(R)
        Volume(i) = Wavg * Time(i-1);
        Time(i) = Volume(i) / R(i);
    end

    TTD = Time(Imax+2);
    if TTD > DT
        notfeasible = 1; 
        %error('Error: downtime check FAILED!');
    end
    
    TTM = sum(Time);
    if (TTM > TM)
        notfeasible = 2;
        %error('Error: Total migration time exceeded');
    end
    
    % figure(1)
    % bar(Volume)
    % title('Grafico del Volume')
    % grid on
    % 
    % figure(2)
    % bar(Time)
    % title('Grafico del Tempo')
    % grid on
