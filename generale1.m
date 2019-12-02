%LAB GEONUMERICA: esercizio 1
clear;
clc

% NOTA: qui non è implementato il "taglio" per la soluzione tre. 
% usare il programma esercizio1_3


% Risoluzione di un problema di Cauchy 

% parametri 
h = 0.1;  % intervallo di integrazione
ti = 0; % tempo iniziale 
tf = 10; % tempo finale 
% t = [ti:h:tf]; %vettore dei nodi su cui calcolare la soluzione temporale (dominio)
u0 = 1;   % condizione iniziale 


funzione = 1; % può essere 1,2 o 3 a seconda della funzione da integrare 
% 1: u_t = u
% 2: u_t=f=ucos(t)
% 3: u_t=f=u+u^2

% calcolo della soluzione con metodi diversi

if funzione == 1
    [uEU,t] =EU(@f1, u0, ti, tf, h);
    [uAB2,t] =AB2(@f1, u0, ti, tf, h);
    [uAB3,t] =AB3(@f1, u0, ti, tf, h);
    [uAB4,t] =AB4(@f1, u0, ti, tf, h);
elseif funzione == 2
    [uEU,t] =EU(@f2, u0, ti, tf, h);
    [uAB2,t] =AB2(@f2, u0, ti, tf, h);
    [uAB3,t] =AB3(@f2, u0, ti, tf, h);
    [uAB4,t] =AB4(@f2, u0, ti, tf, h);
else 
    [uEU,t] =EU(@f3, u0, ti, tf, h);
    [uAB2,t] =AB2(@f3, u0, ti, tf, h);
    [uAB3,t] =AB3(@f3, u0, ti, tf, h);
    [uAB4,t] =AB4(@f3, u0, ti, tf, h);
end 

% Soluzione analitica
if funzione == 1
    u = u0 * exp(t); 
elseif funzione == 2
    u = u0*exp(sin(t)-sin(0));
else 
    if u0 ~= 0
        u = (1)./((1+(1/u0))*exp(-(t-ti))-1);
    else
        u = 0;
    end 
end 

% figure 

% u(t)
figure(1)
subplot(1,3,1)
plot(t,u,t,uEU,t,uAB2,t,uAB3,t,uAB4,'LineWidth',1.5);
grid on;
title('Equazione')
xlabel('t'), ylabel('u(t)')
legend('Esatta','Eulero','AB2','AB3','AB4')


subplot(1,3,2)
semilogy(t,abs(u-uEU'),t,abs(u-uAB2'),t,abs(u-uAB3'),t,abs(u-uAB4'),'LineWidth',1.5);
grid on;
title('Errori assoluti')
xlabel('t'), ylabel('Errore assoluto')
legend('Eulero','AB2','AB3','AB4')

subplot(1,3,3)
semilogy(t,abs(u-uEU')./u,t,abs(u-uAB2')./u,t,abs(u-uAB3')./u,t,abs(u-uAB4')./u,'LineWidth',1.5);
grid on;
title('Errori relativi')
xlabel('t'), ylabel('Errore relativo')
legend('Eulero','AB2','AB3','AB4')


% funzioni integratrici

function [uE,t] = EU(f, u0, ti, tf, h)
%f(u, t) E' la fuinzione da integrare
%u_i E' la condizione iniziale
%ti E' il tempo iniziale
%tf E' il tempo finale
%h E' l'intervallo di integrazione

    t = [ti:h:tf]; %vettore di tempi
    uE = zeros(length(t),1);
    uE(1) = u0;  % condizione iniziale per la soluzione numerica vE 
    for i = 2 : length(t)                      % ciclo for per il calcolo della soluzione numerica vE a ogni step i
        uE(i) = uE(i-1) + h * f(uE(i-1),t(i-1));     
    end
    
end

function [uAB2,t] = AB2(f, u0, ti, tf, h)
    t = [ti:h:tf]; %vettore di tempi
    [uAB2,tb] = EU(f, u0, ti, ti+h, h);
    for i = 3:length(t)
        uAB2(i) = uAB2(i-1) + h*(3/2*f(uAB2(i-1), t(i-1)) -1/2*f(uAB2(i-2), t(i-2)));
    end
end


function [uAB3, t] = AB3(f, u0, ti, tf, h)
    t = [ti:h:tf]; %vettore di tempi
    [uAB3,tb] = AB2(f, u0, ti, ti+2*h, h);
    for i = 4:length(t)
        uAB3(i) = uAB3(i-1) + h*(23/12*f(uAB3(i-1), t(i-1)) -4/3*f(uAB3(i-2), t(i-2)) +5/12*f(uAB3(i-3), t(i-3)));
    end
end

function [uAB4, t] = AB4(f, u0, ti, tf, h)
    t = [ti:h:tf]; %vettore di tempi
    [uAB4,tb] = AB3(f, u0, ti, ti+3*h, h);
    for i = 5:length(t)
        uAB4(i) = uAB4(i-1) + h*(55/24*f(uAB4(i-1), t(i-1)) -59/24*f(uAB4(i-2), t(i-2)) +37/24*f(uAB4(i-3), t(i-3)) - 9/24*f(uAB4(i-4), t(i-4)));
    end
end



% funzioni da integrare

% Soluzione analitica

function d = f1(u, t)
    d = u;
end

function d = f2(u, t)
    d = u*cos(t);
end

function d = f3(u, t)
    d = u + u^2;
end
























