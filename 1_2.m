%LAB GEONUMERICA: esercizio 1
clear;
clc

% Risoluzione del problema u_t=f=ucos(t)

h = 0.1;  % intervallo di integrazione
t = [0:h:10] %vettore dei nodi su cui calcolare la soluzione temporale (dominio)
u0 = 1;   % condizione iniziale 
Nt = length(t)  % numero di elementi nel dominio temporale

% metodo di Eulero esplicito

vE(1) = u0;
for i = 2 : Nt                        % ciclo for per il calcolo della soluzione numerica vE a ogni step i
    vE(i) = vE(i-1) + h * vE(i-1)*cos(t(i-1));     
end

% soluzione analitica
u = u0*exp(sin(t)-sin(0));


% metodo di Eulero implicito

vvE(1) = u0;
for i = 2 : Nt
    vvE(i) = vvE(i-1)/(1-h*cos(t(i)))
end 


figure(1)
subplot(1,3,1)
plot (t,u,t,vE,t,vvE,'LineWidth',1.5)
grid on;
legend('Soluzione analitica' , 'Eulero esplicito', 'Eulero implicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-vE),'LineWidth',1.5)
grid on;
legend('Errore assoluto Eulero esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-vE)./u,'LineWidth',1.5)
grid on;
legend('Errore relativo Eulero esplicito','LineWidth',1.5)


% metodo di AB secondo ordine 
AB1(1) = u0;
AB1(2) = vE(2);
for i = 3 : Nt 
    AB1(i) = AB1(i-1)+h*((3/2)*AB1(i-1)*cos(t(i-1))-1/2*AB1(i-2)*cos(t(i-2)));
end

figure(2)
subplot(1,3,1)
plot (t,u,t,AB1,'LineWidth',1.5)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 2','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-AB1),'LineWidth',1.5)
grid on;
legend('Errore assoluto AB2','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-AB1)./u,'LineWidth',1.5)
grid on;
legend('Errore relativo AB2','LineWidth',1.5)


% metodo AB di terzo ordine
AB2(1) = u0;
AB2(2) = AB1(2);
AB2(3) = AB1(3);
for i = 4 : Nt
    AB2(i) = AB2(i-1) + h*((23/12)*AB2(i-1)*cos(t(i-1))-(4/3)*AB2(i-2)*cos(t(i-2))+(5/12)*AB2(i-3)*cos(t(i-3)));
end

% Plot della soluzione numerica e della soluzione analitica
figure(3)
subplot(1,3,1)
plot (t,u,t,AB2,'LineWidth',1.5)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 3 esplicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-AB2),'LineWidth',1.5)
grid on;
legend('Errore AB3 esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-AB2)./u,'LineWidth',1.5)
grid on;
legend('Errore relativo AB3 esplicito','LineWidth',1.5)

% Metodo AB quarto ordine
AB3(1) = u0; 
AB3(2:4) = AB2(2:4);
for i = 5 : Nt
    AB3(i) = AB3(i-1) + h*((55/24)*AB3(i-1)*cos(t(i-1))-(59/24)*AB3(i-2)*cos(t(i-2))+(37/24)*AB3(i-3)*cos(t(i-3))-(9/24)*AB3(i-4)*cos(t(i-4)));
end

% Plot della soluzione numerica e della soluzione analitica
figure(4)
subplot(1,3,1)
plot (t,u,t,AB3)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 4 esplicito','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(u-AB3),'LineWidth',1.5)
grid on;
legend('Errore AB3 esplicito','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(u-AB3)./u,'LineWidth',1.5)
grid on;
legend('Errore relativo AB4 esplicito','LineWidth',1.5)

% Plot comparativi
figure(5)
subplot(1,2,1)
plot(t,u,t,vE,t,vvE,t,AB1,t,AB2,t,AB3,'LineWidth',1.5)
grid on;
legend('Soluzione analitica', 'Eulero esplicito','Eulero implicito', 'AB 2 ordine', 'AB3 ordine', 'AB 4 ordine')
title('Confronto tra le soluzioni')
subplot(1,2,2)
semilogy(t,abs(u-vE)./u,t,abs(u-AB1)./u,t,abs(u-AB2)./u,t,abs(u-AB3)./u,'LineWidth',1.5)
grid on; 
legend('Eulero esplicito', 'AM 2 ordine', 'AM 3 ordine', 'AM 4 ordine')
title('Confronto tra gli errori relativi')












