% EQUAZIONE DI DIFFUSIONE 

clear;
clc

% definizione di parametri e variabili 
% dominio spaziale e temporale 
xi = -20;
xf = +20; 
h = 0.1; 
ti = 0; 
tf = 100; 
k = 0.1; 

x = [xi:h:xf];
t = [ti:k:tf];
N = length(x);

% parametri
sigma = k/(h^2); 
b = 0.2; 
L = 12;
U = 2.;
W = -U; 

x0 = 0; 
x01 = x0-L;
x02 = x0+L; 


perturbazione = 1; 
% può essere 1,2,3
% 1: box function 
% 2: semi-coseno
% 3: seno

% ricorda la condizione di stabilità 

% inizializzazioni 
u0 = zeros(length(x),1); % inizializzazione condizione iniziale
u = zeros(length(t),length(x)); % inizializzazione soluzione analitica 
u1 = zeros(length(t),length(x)); % inizializzazione soluzione analitica 
u2 = zeros(length(t),length(x)); % inizializzazione soluzione analitica 


% perturbazione iniziale 

if perturbazione == 1
    for i = 1:N
        if abs(x(i)-x0)<(L/2)
        u0(i) = U;
        elseif abs(x(i)-x0)>(L/2)
        u0(i) = 0;
        end     
    end 
    figure(1)
    plot(x,u0,'LineWidth',1.5)
    grid on; 
    ylim([-0.5,2.5])
    title("condizione iniziale")    
elseif perturbazione == 2
    L = xf-xi;
    for i = 1:N
        u0(i)=U*cos(pi/L*x(i));
    end 
    figure(1)
    plot(x,u0,'LineWidth',1.5)
    grid on; 
    ylim([-0.5,2.5])
    title("condizione iniziale")
elseif perturbazione == 3
    L = xf-xi;
    for i = 1:N
        u0(i)=U*sin(2*pi/L*x(i));
    end 
    figure(1)
    plot(x,u0,'LineWidth',1.5)
    grid on; 
    ylim([-2.5,2.5])
    title("condizione iniziale")
end 



% soluzione analitica
% CONDIZIONI AL CONTORNO PER L'ANALITICA CON LA DOPPIA BOX 

if perturbazione == 1
    for i = length(x)
        u(1,i) = u0(i); 
        u1(1,i) = u0(i); 
        u2(1,i) = u0(i); 
    end

    for i = 2:length(t)
        for j = 1:length(x) 
            aux1 = ((x(j)-x0)+L/2)/(2*b*sqrt(t(i)));
            aux2 = (L/2-(x(j)-x0))/(2*b*sqrt(t(i)));
            if abs(x(j)-x0)<=(L/2)
                u(i,j) = U/2*(erf(aux1)+erf(aux2));
            elseif abs(x(j)-x0)>=(L/2)
                u(i,j) = U/2*(erf(aux1)-erf(-aux2));
            end 
        end
        for j = 1:length(x)
            aux1 = ((x(j)-x01)+L/2)/(2*b*sqrt(t(i)));
            aux2 = (L/2-(x(j)-x01))/(2*b*sqrt(t(i)));
            if (x(j)-x01)<=(L*3/2) && (x(j)-x01)>=(L*1/2)
                u1(i,j) = W/2*(erf(aux1)+erf(aux2));
            else
                u1(i,j) = W/2*(erf(aux1)-erf(-aux2));
            end 
        end
        for j = 1:length(x)
            aux1 = ((x(j)-x02)+L/2)/(2*b*sqrt(t(i)));
            aux2 = (L/2-(x(j)-x02))/(2*b*sqrt(t(i)));
            if (x(j)-x02)<=(-L*1/2) && (x(j)-x02)>=(-L*3/2)
                u2(i,j) = W/2*(erf(aux1)+erf(aux2));
            else 
                u2(i,j) = W/2*(erf(aux1)-erf(-aux2));
            end 
        end
        
        for j = 1:length(x)
            utot(i,j) = u(i,j)+u1(i,j)+u2(i,j); 
        end 
        
        figure(2)
        plot(x,u(i,:),x,u0,'LineWidth',1.5)
        grid on; 
        %ylim([-0.5,2.5])
        title("soluzione analitica");
    end  
    
    
    
    
elseif perturbazione == 2
    for i = 1:length(t)
        for j = 1:length(x)
            u(i,j) = U*exp(-(pi^2)/(L^2)*b^2*t(i))*cos(pi/L*x(j));
        end
        figure(2)
        plot(x,u(i,:),x,u0,'LineWidth',1.5)
        grid on; 
        ylim([-0.5,2.5])
        title("soluzione analitica");
    end 
    
elseif perturbazione == 3
    for i = 1:length(t)
        for j = 1:length(x)
            u(i,j) = U*exp(-(4*pi^2)/(L^2)*b^2*t(i))*sin(2*pi/L*x(j));
        end
        figure(2)
        plot(x,u(i,:),x,u0,'LineWidth',1.5)
        grid on; 
        ylim([-2.5,2.5])
        title("soluzione analitica");
    end 
    
end 




figure(3)
surf(x,t,u,'LineStyle','None')

% risoluzione numerica della PDE 

% AGGIUNGERE LE CONDIZIONI AL CONTORNO PER LA ANALITICA 
% CON IL METODO DELLE DOPPIE BOX

%% Metodo Eulero 

vEU = zeros(length(t),length(x));
vEU(1,:) = u0;

% condizioni al contorno 
for i = 1:length(t)
    vEU(i,1)= 0;
    vEU(i,N)= 0;
end 

for i = 2:(length(t)-1)
    for j = 2:(length(x)-1)
       vEU(i,j)= vEU(i-1,j)+b^2*sigma*(vEU(i-1,j+1)-2*vEU(i-1,j)+vEU(i-1,j-1));     
    end 
    figure(4)
    plot(x,u(i,:),x,vEU(i,:),'LineWidth',1.5);
    title("metodo di Eulero");
    
    err_EU = errore(u,vEU);
end 



%% Metodo Backward Eulero 

vBEU = zeros(length(t),length(x));
vBEU(1,:) = u0;

% definizione delle matrici A e B 
% (perchè è un metodo implicito)

A = zeros(length(x),length(x));
A(N,N) = 1;
A(1,1) = 1;

for i = 2:length(x)-1
    A(i,i)=1+2*b^2*sigma;
end 
for i = 2:(length(x)-1)
        A(i,i+1) = -b^2*sigma;
end
for i = 1:(length(x)-1)
        A(i+1,i) = -b^2*sigma;
end

B = eye(length(x));
B(1,1) = 0;
B(N,N) = 0;

for i = 2:(length(t)-1)
    vBEU(i,:) = A\B*vBEU(i-1,:)';
    figure(4)
    plot(x,u(i,:),x,vBEU(i,:),'LineWidth',1.5);
    title("metodo Backward Eulero");
end 


%% Metodo di Cranck-Nicolson 

vCN = zeros(length(t),length(x));
vCN(1,:) = u0;

% matrici A e B 
A = eye(length(x));
for i = 2:length(x)-1
    A(i,i-1) = -b^2*sigma/2;
    A(i,i) = 1+2*b^2*sigma/2;
    A(i,i+1) = -b^2*sigma/2;
end

B = zeros(length(x));
for i = 2:length(x)-1
    B(i,i-1) = +b^2*sigma/2;
    B(i,i) = 1-2*b^2*sigma/2;
    B(i,i+1) = +b^2*sigma/2;
end

for i = 2:(length(t)-1)
    vCN(i,:) = A\B*vCN(i-1,:)';
    figure(5)
    plot(x,u(i,:),x,vCN(i,:),'LineWidth',1.5);
    title("metodo Cranck Nicolson");
    ylim([-0.5,2.5])
end 


%% Metodo Leap Frog 
% nota: è sempre instabile, è normale che faccia così

vLF = zeros(length(t),length(x));
vLF(1,:) = u0;
vLF(2,:) = vEU(2,:);

% condizioni al contorno 
for i = 1:length(t)
    vLF(i,1)= 0;
    vLF(i,N)= 0;
end 

for i = 3:(length(t)-1)
    for j = 2:(length(x)-1)
       vLF(i,j)= vLF(i-2,j)+2*b^2*sigma*(vLF(i-1,j+1)-2*vLF(i,j)+vLF(i-1,j-1));     
    end 
        figure(6)
        plot(x,u(i,:),x,vLF(i,:),'LineWidth',1.5);
        title("metodo Leap Frog");
end 



%% definizione di funzioni
% funzione per il calcolo dell'errore

function err = errore(u,v) 
    % u = soluzione analitica
    % v = metodo numerico 
    diff = (u - v).^2; 
    num = sum(diff,2); % il 2 è per sommare lungo le righe 
    den = sum(u.^2,2);
    err = sqrt(num./den);
end 



















