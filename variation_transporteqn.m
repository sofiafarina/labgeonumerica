% ** EQUAZIONE DEL TRASPORTO **
% esercizio numero tre - laboratorio di geofisica numerica 

% variante: qui c'è la c che varia nel tempo (cresce)

clear;
clc

% dominio spaziale e temporale 
xi = -10;
xf = +10; 
h = 0.1; 
ti = 0; 
tf = 10; 
k = 0.1; 
lambda = k/h; 
x = [xi:h:xf];
t = [ti:k:tf];
L = 12;
U = 5; % ampiezza iniziale dell'oscillazione 
c = linspace(0.45,0.80,length(t)); % velocità di fase

%%

N = length(x);

u0 = zeros(length(x),1); % inizializzazione condizione iniziale
u = zeros(length(t),length(x)); % inizializzazione soluzione analitica 

% pay attention: lambda esattamente uguale a 1 è un caso particolare 
% per il metodo upwind perchè si sta proprio al limite della instabilità
% dunque allontanati da 1 

% perturbazione iniziale 
for i = 1:length(x)
    if abs(x(i))<(L/2)
    u0(i) = U/2*(1+cos(2*pi*x(i)/L));
    elseif abs(x(i))>(L/2)
    u0(i) = 0;
    end 
end 

figure(1)
plot(x, u0,'LineWidth',1.1)

%% 
% soluzione analitica 
% u(x,t) = uo(x+ct)

for i = 2:length(t)
    for j = 1:length(x)
        if abs(x(j)+c(i)*t(i))<(L/2)
        u(i,j) = U/2*(1+cos(2*pi*(x(j)+c(i)*t(i))/L));
        elseif abs(x(j)+c(i)*t(i))>(L/2)
        u(i,j) = 0;
        end 
    end
    figure(2)
    plot(x,u(i,:),'LineWidth',1.1)
    title("soluzione analitica");
end 
    
figure(3)
surf(x,t,u, 'Linestyle', 'None');
title("soluzione analitica");


%%
% soluzione numerica della PDE 
% METODO UPWIND o REVERSE UPWIND 

v = zeros(length(t),length(x));
v(1,:) = u0;

if c > 0
    % condizione al contorno (per c > 0)
    for i = 1:length(t)
        v(i,N)= 0;
    end 

    for i = 2:(length(t)-1)
        v(i,1)=v(i-1,1)+c(i)*lambda*(v(i-1,2)-v(i-1,1)); % condizione al contorno 
        for j = 1:(length(x)-1)
            v(i,j)= v(i-1,j)+c(i)*lambda*(v(i-1,j+1)-v(i-1,j));     
        end 
        figure(4)
        plot(x,u(i,:),x,v(i,:),'LineWidth',1.1);
        title("metodo upwind");
        pause(0.01)
    end 

    figure(5)
    surf(x,t,v, 'Linestyle', 'None');
    title("metodo upwind");
    
    % calcolo dell'errore
    err_UW = errore(u,v);
    
else 
    % condizione al contorno (per c < 0)
    for i = 1:length(t)
        v(i,1)=0;
    end 

    for i = 2:(length(t))
        v(i,N)=v(i-1,N)+c(i)*lambda*(v(i-1,N)-v(i-1,N-1)); % condizione al contorno 
        for j = 2:(length(x))
            v(i,j)=v(i-1,j)+c(i)*lambda*(v(i-1,j)-v(i-1,j-1));
        end 
        figure(4)
        plot(x,u(i,:),x,v(i,:),'LineWidth',1.1);
        title("metodo reverse upwind");
        pause(0.1)
    end 

    figure(5)
    surf(x,t,v, 'Linestyle', 'None');
    title("metodo reverse upwind");
    
    % calcolo dell'errore
    err_RUW = errore(u,v);
end 

 
 

%%
% METODO BACKWARD EULERO
% ricorda che essendo un metodo implicito le condizioni al contorno sono 
% già contenute nelle matrici A e B 

vBEU = zeros(length(t),length(x));
vBEU(1,:) = u0; % inizializzazione con la condizione iniziale 

% matrici per il metodo implicito
A = zeros(length(x),length(x));

for i = 1:length(x)
    A(i,i)=1.;
end 
for i = 2:(length(x)-1)
        A(i,i+1) = -c(1)*lambda*0.5;
end
for i = 1:(length(x)-1)
        A(i+1,i) = c(1)*lambda*0.5;
end
A(length(x),length(x)-1)=0.;

B = zeros(length(x),length(x));
if c(1) > 0
    % caso c > 0 
    for i = 1:(length(x)-1)
        B(i,i)=1.;
    end 
    B(1,1)=1-c(1)*lambda; 
    B(1,2)=c(1)*lambda;
else 
    % caso c < 0 
    for i = 1:(length(x)-1)
        B(i,i)=1.;
    end 
    B(1,1)=0; 
    B(length(x),length(x))=1+c(1)*lambda;
    B(length(x),length(x)-1)=-c(1)*lambda;
end 

for i = 2:(length(t))
    vBEU(i,:) = inv(A)*B*vBEU(i-1,:)';
  
    figure(8)
    plot(x,u(i,:),x,vBEU(i,:),'LineWidth',1.1);  
    title("metodo backward Eulero");
    pause(0.001)
end 

% calcolo dell'errore
err_BE = errore(u,vBEU);


%%
% METODO LEAP FROG

vLF = zeros(length(t),length(x));
vLF(1,:) = u0;
vLF(2,:) = v(2,:);

% condizioni al contorno
if c(1) > 0 
    for i = 1:length(t)
    vLF(i,N)= 0;
    end 
else 
    for i = 1:length(t)
    vLF(i,1)=0;
    end       
end 

for i = 3:(length(t))
    % condizioni al contorno 
    if c(1) > 0
        vLF(i,1)=vLF(i-1,1)+c(i)*lambda*(vLF(i-1,2)-vLF(i-1,1)); 
    else 
        vLF(i,N)=vLF(i-1,N)+c(i)*lambda*(vLF(i-1,N)-vLF(i-1,N-1));
    end 
     
    for j = 2:(length(x)-1)
        vLF(i,j)= vLF(i-2,j)+c(i)*lambda*(vLF(i-1,j+1)-vLF(i-1,j-1));     
    end 
    
    figure(9)
    plot(x,u(i,:),x,vLF(i,:),'LineWidth',1.1);
    title("metodo leap frog");
    pause(0.001)
end 

% calcolo dell'errore
err_LF = errore(u,vLF);

% nota: se scegli un h più piccolo si vede meglio l'oscillazione che fa 
% alla fine della "campana" che è proprio un errore che deriva da come 
% abbiamo scritto il metodo e che viene essenzialmente eliminata in LW 



%%
% METODO LAX WENDROFF 

vLW = zeros(length(t),length(x));
vLW(1,:) = u0;

% condizioni al contorno
if c(1) > 0 
    for i = 1:length(t)
    vLW(i,N)= 0;
    end 
else 
    for i = 1:length(t)
    vLW(i,1)=0;
    end       
end 

for i = 2:(length(t))
    % condizioni al contorno 
    if c(1) > 0
        vLW(i,1)=vLW(i-1,1)+c(i)*lambda*(vLW(i-1,2)-vLW(i-1,1)); 
    else 
        vLW(i,N)=vLW(i-1,N)+c(i)*lambda*(vLW(i-1,N)-vLW(i-1,N-1));
    end 
     
    for j = 2:(length(x)-1)
        vLW(i,j)= vLW(i-1,j)+c(i)*lambda*0.5*(vLW(i-1,j+1)-vLW(i-1,j-1))+c(i)^2*lambda^2*0.5*(vLW(i-1,j+1)-2*vLW(i-1,j)+vLW(i-1,j-1));     
    end 
    
    figure(10)
    plot(x,u(i,:),x,vLW(i,:),'LineWidth',1.1);
    title("metodo Lax Wendroff");
    pause(0.001)
end 

% calcolo dell'errore
err_LW = errore(u,vLW);

%%
% METODO CRANCK NICOLSON 
% qui servono le giuste matrici A e B 
% controlla se le matrici devono essere diverse per c < o > 0 

vCN = zeros(length(t),length(x));
vCN(1,:) = u0; % inizializzazione con la condizione iniziale 

% matrici per il metodo implicito
% CONTROLLARE LE MATRICI 
A = eye(length(x));
for i = 2:length(x)-1
    A(i,i-1) = c(1)*lambda/4;
    A(i,i+1) = -c(1)*lambda/4;
end

B = eye(length(x));
if c(1) > 0
    B(1,1) = 1 - c(1)*lambda;
    B(1,2) = c(1)*lambda;
    B(length(x),length(x)) = 0;
    for i = 2:length(x)-1
        B(i,i-1) = -c(1)*lambda/4;
        B(i,i+1) = c(1)*lambda/4;
    end
else
    B(1,1)=0;
    B(length(x),length(x)-1) = -c(1)*lambda;
    B(length(x),length(x)) = 1+c(1)*lambda;
    for i = 2:length(x)-1
        B(i,i-1) = -c(1)*lambda/4;
        B(i,i+1) = c(1)*lambda/4;
    end
end

for i = 2:(length(t))
    vCN(i,:) = A\B*vCN(i-1,:)';
  
    figure(11)
    plot(x,u(i,:),x,vCN(i,:),'LineWidth',1.1);  
    ylim([-0.5,5.5]);
    title("metodo Cranck Nicholson");
    pause(0.01)
end 

% calcolo dell'errore
err_CN = errore(u,vCN);

%%
% plot comparativo degli errori

if c(1) > 0
    figure(12)
    plot([1:1:length(err_UW)-1],err_UW(1:length(err_UW)-1),[1:1:length(err_BE)],err_BE,[1:1:length(err_LF)],err_LF,[1:1:length(err_CN)],err_CN,[1:1:length(err_LW)],err_LW,'LineWidth',1.1)
    legend('Upwind','Backward Eulero','Leap Frog','Cranck Nicolson','Lax Wendroff')
else 
    figure(12)
    plot([1:1:length(err_RUW)-1],err_RUW(1:length(err_RUW)-1),[1:1:length(err_BE)],err_BE,[1:1:length(err_LF)],err_LF,[1:1:length(err_CN)],err_CN,[1:1:length(err_LW)],err_LW,'LineWidth',1.1)
    legend('Reverse upwind','Backward Eulero','Leap Frog','Cranck Nicolson','Lax Wendroff')
end 


%%

% definizione di funzioni
% funzione per il calcolo dell'errore

function err = errore(u,v) 
    % u = soluzione analitica
    % v = metodo numerico 
    diff = (u - v).^2; 
    num = sum(diff,2); % il 2 è per sommare lungo le righe 
    den = sum(u.^2,2);
    err = sqrt(num./den);
end 









