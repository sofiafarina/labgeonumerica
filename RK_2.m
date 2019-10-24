% LAB GEONUMERICA 
% metodi di runge kutta
% esercitazione numero due 
clear;
clc

% COSE DA FARE: 
% provare a fare una tendina per i vari metodi 
% fare i confronti con AB3 o AB4, in termini di errore ad esempio 

% equazione 1: u_t = u*cos(t)

h = 0.01; %intervallo di tempo prova anche 0.1, 0.001, 0.0001
t_i = 0; %tempo iniziale
t_f = 10; %Tempo finale
u_i = 1; %condizione iniziale
t = [t_i:h:t_f]; %vettore di tempi
u(1) = u_i; % inizializzazione condizione iniziale

% soluzione esatta
uex = u_i*exp(sin(t)-sin(0));


%RK 2 ORDINE  

tic
% scelta del metodo 
meth = 0; % può essere 0 oppure 1
if meth == 0  % metodo di heun
    c = [0,1];
    a = [0,0; 1,0];
    b = [0.5,0.5];
else   % metodo del punto medio 
    c = [0,0.5];
    a = [0,0; 1,0];
    b = [0.5,0.5];
end 
toc

s = length(c); % ciclo interno al ciclo n 
tau = zeros(s,1); % valori tempi intermedi
utilde = zeros(s,1); % valori provvisori intermedi
ftilde = zeros(s,1); % valori provvisori della funzione
u2(1)=u_i;

for i = 2:length(t)
    for j = 1:s
        tau(j) = t(i-1)+h*c(j);
        utilde(j) = u2(i-1)+h*a(j,:)*ftilde; 
        ftilde(j) = utilde(j)*cos(tau(j));
    end
    u2(i) = u2(i-1) + h*b*ftilde;
end 

figure(1)
subplot(1,3,1)
plot(t,uex,t,u2,'LineWidth',1.5)
grid on;
legend('Soluzione esplicita', 'RK2','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(uex-u2),'LineWidth',1.5)
grid on;
legend('Errore assoluto','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(uex-u2)./u2,'LineWidth',1.5)
grid on;
legend('Errore relativo','LineWidth',1.5)
sgtitle('Runge Kutta 2° ordine')
toc 

% RK 3 ORDINE 

tic
% scelta del metodo 
method = 0; % può essere 0 oppure 1
if method == 0  % metodo di heun (3)
    c = [0,1/3,2/3];
    a = [0,0,0;1/3,0,0;0,2/3,0];
    b = [1/4,0,3/4];
elseif method == 1  % metodo RK 3
    c = [0,0.5,1];
    a = [0,0,0;0.5,0,0;0,1,0];
    b = [1/3,1/3,1/3];
elseif method == 2  % ODE 23 
    c = [0,0.5,3/4];
    a = [0,0,0;0.5,0,0;0,3/4,0];
    b = [2/9,1/3,4/9];
end 
toc

s = length(c); % ciclo interno al ciclo n 
tau = zeros(s,1); % valori tempi intermedi
utilde = zeros(s,1); % valori provvisori intermedi
ftilde = zeros(s,1); % valori provvisori della funzione
u3(1)=u_i;

for i = 2:length(t)
    for j = 1:s
        tau(j) = t(i-1)+h*c(j);
        utilde(j) = u3(i-1)+h*a(j,:)*ftilde; 
        ftilde(j) = utilde(j)*cos(tau(j));
    end
    u3(i) = u3(i-1) + h*b*ftilde;
end 

figure(2)
subplot(1,3,1)
plot(t,uex,t,u3,'LineWidth',1.5)
grid on;
legend('Soluzione esplicita', 'RK3','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(uex-u3),'LineWidth',1.5)
grid on;
legend('Errore assoluto','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(uex-u3)./u3,'LineWidth',1.5)
grid on;
legend('Errore relativo','LineWidth',1.5)
sgtitle('Runge Kutta 3° ordine')
toc


% RK 4 ORDINE 
tic 

c = [0,0.5,0.5,1];
a = [0,0,0,0;0.5,0,0,0;0,0.5,0,0;0,0,1,0];
b = [1/6,1/3,1/3,1/6];
toc

s = length(c); % ciclo interno al ciclo n 
tau = zeros(s,1); % valori tempi intermedi
utilde = zeros(s,1); % valori provvisori intermedi
ftilde = zeros(s,1); % valori provvisori della funzione
u4(1)=u_i;

for i = 2:length(t)
    for j = 1:s
        tau(j) = t(i-1)+h*c(j);
        utilde(j) = u4(i-1)+h*a(j,:)*ftilde; 
        ftilde(j) = utilde(j)*cos(tau(j));
    end
    u4(i) = u4(i-1) + h*b*ftilde;
end 


figure(3)
subplot(1,3,1)
plot(t,uex,t,u3,'LineWidth',1.5)
grid on;
legend('Soluzione esplicita', 'RK4','LineWidth',1.5)
subplot(1,3,2)
plot(t,abs(uex-u3),'LineWidth',1.5)
grid on;
legend('Errore assoluto','LineWidth',1.5)
subplot(1,3,3)
semilogy(t,abs(uex-u3)./u3,'LineWidth',1.5)
grid on;
legend('Errore relativo','LineWidth',1.5)
sgtitle('Runge Kutta 4° ordine')
toc 




























