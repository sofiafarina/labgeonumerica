%LAB GEONUMERICA: esercizio 1
clear;
clc

% Risoluzione del problema u_t=f=ucos(t)

h = 0.01;  % intervallo di integrazione
t = [0:h:10]; %vettore dei nodi su cui calcolare la soluzione temporale (dominio)
u0 = 1;   % condizione iniziale 
Nt = length(t);  % numero di elementi nel dominio temporale
t0 = 0;


% soluzione analitica
if u0 ~= 0
    u = (1)./((1+(1/u0))*exp(-(t-t0))-1);
else
    u = 0;
end 


% metodo di Eulero esplicito

vE(1) = u0;
i1 = 0;
for i = 2 : Nt 
    vE(i) = vE(i-1) + h * (vE(i-1)+vE(i-1)^2);
    if vE(i) > 100.
        i1 = i; 
        disp(i1);
        break
    end 
end


figure(1)
subplot(1,3,1)
plot (t,u,t(1:i1),vE)
grid on;
legend('Soluzione analitica' , 'Eulero esplicito')
subplot(1,3,2)
plot(t(1:i1),abs(u(1:i1)-vE))
grid on;
legend('Errore assoluto Eulero esplicito')
subplot(1,3,3)
semilogy(t(1:i1),abs(u(1:i1)-vE)./u(1:i1))
grid on;
legend('Errore relativo Eulero esplicito')


% metodo di Eulero implicito
% FARE 

% metodo di AB secondo ordine 
AB1(1) = u0;
AB1(2) = vE(2);
i2 = 0; 
for i = 3 : Nt 
    AB1(i) = AB1(i-1)+h*((3/2)*(AB1(i-1)+AB1(i-1)^2)-1/2*(AB1(i-2)+AB1(i-2)^2));
    if AB1(i) > 100
        i2 = i;
        break 
    end 
end

figure(2)
subplot(1,3,1)
plot (t,u,t(1:i2),AB1)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 2')
subplot(1,3,2)
plot(t(1:i2),abs(u(1:i2)-AB1(1:i2)))
grid on;
legend('Errore assoluto AB2')
subplot(1,3,3)
semilogy(t(1:i2),abs(u(1:i2)-AB1(1:i2))./u(1:i2))
grid on;
legend('Errore relativo AB2')



% metodo AB di terzo ordine
AB2(1) = u0;
AB2(2) = AB1(2);
AB2(3) = AB1(3);
i3 = 0;
for i = 4 : Nt
    AB2(i) = AB2(i-1) + h*((23/12)*(AB2(i-1)+AB2(i-1)^2)-(4/3)*(AB2(i-2)+AB2(i-2)^2)+(5/12)*(AB2(i-3)+AB2(i-3)^2));
    if AB2(i) > 100
        i3 = i;
        break 
    end 
end

% Plot della soluzione numerica e della soluzione analitica
figure(3)
subplot(1,3,1)
plot (t,u,t(1:i3),AB2)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 3 esplicito')
subplot(1,3,2)
plot(t(1:i3),abs(u(1:i3)-AB2(1:i3)))
grid on;
legend('Errore AB3 esplicito')
subplot(1,3,3)
semilogy(t(1:i3),abs(u(1:i3)-AB2(1:i3))./u(1:i3))
grid on;
legend('Errore relativo AB3 esplicito')

% Metodo AB quarto ordine
AB3(1) = u0; 
AB3(2:4) = AB2(2:4);
i4 = 0;

for i = 5 : Nt
    AB3(i) = AB3(i-1) + h*((55/24)*(AB3(i-1)+AB3(i-1)^2)-(59/24)*(AB3(i-2)+AB3(i-2)^2)+(37/24)*(AB3(i-3)+AB3(i-3)^2)-(9/24)*(AB3(i-4)+AB3(i-4)^2));
    if AB3(i) > 100
       i4 = i;
       break 
    end 
end

% Plot della soluzione numerica e della soluzione analitica
figure(4)
subplot(1,3,1)
plot (t,u,t(1:i4),AB3)
grid on;
legend('Soluzione analitica' , 'Adams-Bashforth 4 esplicito')
subplot(1,3,2)
plot(t(1:i4),abs(u(1:i4)-AB3(1:i4)))
grid on;
legend('Errore AB3 esplicito')
subplot(1,3,3)
semilogy(t(1:i4),abs(u(1:i4)-AB3(1:i4))./u(1:i4))
grid on;
legend('Errore relativo AB4 esplicito')












