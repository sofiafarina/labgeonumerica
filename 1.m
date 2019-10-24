%LAB GEONUMERICA: esercizio 1
clear;
clc

% Risoluzione di un problema di Cauchy u_t = u


h = 0.1;  % intervallo di integrazione
t = [0:h:10] %vettore dei nodi su cui calcolare la soluzione temporale (dominio)
u0 = 1;   % condizione iniziale 
Nt = length(t)  % numero di elementi nel dominio temporale

% Metodo di Eulero esplicito
vE(1) = u0;  % condizione iniziale per la soluzione numerica vE 
for i = 2 : Nt                        % ciclo for per il calcolo della soluzione numerica vE a ogni step i
    vE(i) = vE(i-1) + h * vE(i-1);     
end

% Soluzione analitica
u = u0 * exp(t);   % calcolo della soluzione analitica

% Metodo di Eulero implicito
vvE(1) = u0;
for i = 2 : Nt                        % ciclo for per il calcolo della soluzione numerica vE a ogni step i
     vvE(i) = vvE(i-1)/(1 - h);     
end

% Plot della soluzione numerica e della soluzione analitica
figure(1)
subplot(1,3,1)
plot (t,u,t,vE,t,vvE,'LineWidth',1.5)
legend('Soluzione analitica' , 'Eulero esplicito', 'Eulero implicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-vE),'LineWidth',1.5)
legend('Errore assoluto Eulero esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-vE)./u,'LineWidth',1.5)
legend('Errore relativo Eulero esplicito','LineWidth',1.5)

%%

% Metodi espliciti di Adams-Bashforth secondo ordine 

wE(1) = u0;   % condizione iniziale per la soluzione numerica wE
wE(2) = vE(2);
for i = 3 : Nt                        % ciclo for per il calcolo della soluzione numerica wE a ogni step i
    wE(i) = wE(i-1) + h * ((3/2)*wE(i-1)-(1/2)*wE(i-2));     
end

% Plot della soluzione numerica e della soluzione analitica
figure(2)
subplot(1,3,1)
plot (t,u,t,wE,'LineWidth',1.5)
legend('Soluzione analitica' , 'Adams-Bashforth 2 esplicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-wE),'LineWidth',1.5)
legend('Errore AB2 esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-wE)./u,'LineWidth',1.5)
legend('Errore relativo AB2 esplicito','LineWidth',1.5)


% Metodi espliciti di Adams-Bashforth terzo ordine
uE(1) = u0;
uE(2:3) = wE(2:3);
for i = 4 : Nt
    uE(i) = uE(i-1) + h*((23/12)*uE(i-1)-(4/3)*uE(i-2)+(5/12)*uE(i-3));
end

% Plot della soluzione numerica e della soluzione analitica
figure(3)
subplot(1,3,1)
plot (t,u,t,uE,'LineWidth',1.5)
legend('Soluzione analitica' , 'Adams-Bashforth 3 esplicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-uE),'LineWidth',1.5)
legend('Errore AB3 esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-uE)./u,'LineWidth',1.5)
legend('Errore relativo AB3 esplicito','LineWidth',1.5)

% Metodi espliciti di Adams-Bashforth quarto ordine
zE(1) = u0; 
zE(2:4) = uE(2:4);
for i = 5 : Nt
    zE(i) = zE(i-1) + h*((55/24)*zE(i-1)-(59/24)*zE(i-2)+(37/24)*zE(i-3)-(9/24)*zE(i-4));
end

% Plot della soluzione numerica e della soluzione analitica
figure(4)
subplot(1,3,1)
plot (t,u,t,zE,'LineWidth',1.5)
legend('Soluzione analitica' , 'Adams-Bashforth 4 esplicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-zE),'LineWidth',1.5)
legend('Errore AB3 esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-zE)./u,'LineWidth',1.5)
legend('Errore relativo AB4 esplicito','LineWidth',1.5)

% Plot comparativi
figure(5)
subplot(1,2,1)
plot(t,u,t,vE,t,vvE,t,wE,t,uE,t,zE,'LineWidth',1.5)
legend('Soluzione analitica', 'Eulero esplicito', 'Eulero implicito','AM 2 ordine', 'AM 3 ordine', 'AM 4 ordine','LineWidth',1.5)
title('Confronto tra le soluzioni')
subplot(1,2,2)
semilogy(t,abs(u-vE)./u,t,abs(u-wE)./u,t,abs(u-uE)./u,t,abs(u-zE)./u,'LineWidth',1.5,'LineWidth',1.5)
legend('Eulero esplicito', 'AM 2 ordine', 'AM 3 ordine', 'AM 4 ordine','LineWidth',1.5)
title('Confronto tra gli errori relativi')









