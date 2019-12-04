% LAB GEONUMERICA 
% l'equazione delle onde
% esercitazione numero cinque  
clear;
clc

% EQUAZIONE: 
% u_{tt} = c^2u_{xx}
% u(x,0) = u_0(x)
% u_t(x,0) = u_1(x)
% u è uno spostamento, c è la velocità di fase dell'onda

% variabili e parametri 

xi = -20;
xf = +20; 
ti = 0; 
tf = 200; 
k = 0.1; 
h = 0.15; 
lambda = k/h; 
c = 1.; 

x = [xi:h:xf]; % vettore delle posizioni
t = [ti:k:tf]; % vettore dei tempi 
N = length(x);

% parametri
sigma = k/(h^2); 
b = 0.2; 
L = 40; % dominio spaziale
 

lenx = length(x);
lent = length(t);

% selezione del metodo:
% può essere: 
% (1) trasparenza, 
% (2) Dirichlet,
% (3) Neumann

method = 1; 

% selezione delle condizioni al contorno:
% può essere:
% ('a') gaussiana con derivata iniziale nulla 
% ('b') spostamento iniziale nullo e derivata diversa da 0
% ('c') tenendo solo una delle due componenti dell'onda 




% soluzione esatta
% qui ci sono le tre soluzioni analitiche 
% relative alle tre diverse condizioni iniziali
% in termini di u_0 e u_1 

% uex(i,j)=uex(t,x)

switch method
    case 1
        % inizializzo tutto
        uex = zeros(length(t),length(x));
        u_0 = zeros(length(x),1);
        u_1 = zeros(length(x),1);
        
        for i = 1:length(x)
            u_0(i) = exp(-x(i).^2);
        end 
        
        for i = 1:length(t)
            
            for j = 1:length(x)
            uex(i,j) = 0.5*(exp(-(x(j)+c*t(i))^2)+exp(-(x(j)-c*t(i))^2)); 
            end
            if i < 350
                figure(1)
                plot(x,uex(i,:),'LineWidth',1.1);
                ylim([-0.2,1.1])
                title("caso 1");
            end 
        end        

    case 2
        % inizializzo tutto
        uex = zeros(length(t),length(x));
        u_0 = zeros(length(x),1);
        u_1 = zeros(length(x),1);
        
        for i = 1:length(x)
            u_1(i) = -2*c*x(i).*exp(-x(i).^2);
        end 
        
        for i = 1:length(t)
            
            for j = 1:length(x)
            uex(i,j) = 0.5*(exp(-(x(j)+c*t(i))^2)-exp(-(x(j)-c*t(i))^2)); 
            end
            
            figure(1)
            plot(x,uex(i,:),'LineWidth',1.1);
            ylim([-0.8,0.8])
            title("caso 2");
        end 
        
    case 3
        % inizializzo tutto
        uex = zeros(length(t),length(x));
        u_0 = zeros(length(x),1);
        u_1 = zeros(length(x),1);
        
        for i = 1:length(x)
            u_0(i) = exp(-x(i).^2);
            u_1(i) = +2*c*x(i).*exp(-x(i).^2);
        end 
        
        for i = 1:length(t)
            
            for j = 1:length(x)
            uex(i,j) = (exp(-(x(j)-c*t(i))^2)); 
            end
            
            figure(1)
            plot(x,uex(i,:),'LineWidth',1.1);
            ylim([-0.2,1.2])
            title("caso 3");
        end       

    otherwise
        disp('Le opzioni sono 1,2 e 3')
end

figure(2)
plot(x,u_0,x,u_1,'LineWidth',1.5);
title('perturbazioni iniziali');
legend('u_0(x)','u_1(x)');

figure(3)
title('soluzione analitica');
surf(x,t,uex,'LineStyle','None')




%%

boundary = 'b'; 

% Soluzione per integrazione
% inizializzazione 
v = zeros(length(t),length(x));

% condizioni iniziali 
% definizione dei primi due step temporali 

switch method
    case 1
        % t = 1
        for j = 1:length(x)
            v(1,j) = exp(-x(j).^2);
        end 
        % t = 2
        for j = 1:length(x)
            v(2,j) = exp(-x(j).^2) + k*0;
        end 
    case 2
        % t = 1
        % rimane tutto zero         
        % t = 2
        for j = 1:length(x)
            v(2,j) = 0 - k*2*c*x(j).*exp(-x(j).^2);
        end
  
    case 3
        % t = 1
        for j = 1:length(x)
            v(1,j) =  exp(-x(j).^2);
        end 
        % t = 2
        for j = 1:length(x)
            v(2,j) = exp(-x(j).^2) + k*2*c*x(j).*exp(-x(j).^2);
        end 
end 


% integrazione e definizione delle condizioni al contorno (tre opzioni)

switch boundary
    case 'a' % trasparenza
        for i = 3:length(t)
            % condizioni al contorno
            v(i,1) = v(i-1,1)+c*lambda*(v(i-1,2)-v(i-1,1));
            v(i,length(x)) = v(i-1,length(x))-c*lambda*(v(i-1,length(x))-v(i-1,length(x)-1));
            
            for j = 2:length(x)-1
                v(i,j)=2*v(i-1,j)-v(i-2,j)+c^2*lambda^2*(v(i-1,j+1)-2*v(i-1,j)+v(i-1,j-1));
            end
            figure(4)
            plot(x,uex(i,:),x,v(i,:),'LineWidth',1.5);
            grid on 
            title("condizione di trasparenza");
            ylim([-0.2,1.1])
            disp(t(i))
        end 
    case 'b' % dirichlet
        for i = 3:length(t)
            % condizioni al contorno
            v(i,1) = 0;
            v(i,length(x)) = 0;
            
            for j = 2:length(x)-1
                v(i,j)=2*v(i-1,j)-v(i-2,j)+c^2*lambda^2*(v(i-1,j+1)-2*v(i-1,j)+v(i-1,j-1));
            end
            figure(4)
            plot(x,uex(i,:),x,v(i,:),'LineWidth',1.5);
            grid on 
            title("condizione di Dirichlet");
            ylim([-1.1,1.1])
        end         
    case 'c' % von neumann 
        for i = 3:length(t)
            
            for j = 2:length(x)-1
                v(i,j)=2*v(i-1,j)-v(i-2,j)+c^2*lambda^2*(v(i-1,j+1)-2*v(i-1,j)+v(i-1,j-1));
            end
            % condizioni al contorno
            v(i,1) = v(i,2);
            v(i,length(x)) = v(i,length(x)-1);
            
            figure(4)
            plot(x,uex(i,:),x,v(i,:),'LineWidth',1.5);
            grid on 
            title("condizione di Von Neumann");
            ylim([-0.1,1.0])
            
        end        
        
    otherwise
        disp('Le alternative sono: a,b o c')
end

%% calcolo dell'errore

switch boundary
    case 'a'
        %
    case 'b'
        % calcolo dei tempi del passaggio al centro 
        intervallo = (xf-xi)/c;  
        pos = round(intervallo/k); 
        posizioni = zeros(3,1);
        posizioni(1) = 1; 
        for i = 1:2
            posizioni(i+1) = pos*i+1; 
        end 
        
        errori = zeros(3,1);
        for i = 1:3
            errori(i) = errore(uex(1,:),abs(v(posizioni(i),:)));
        end 
        disp(errori)
        
        % per vederli
        figure(5)
        plot(x,uex(1,:),x,abs(v(posizioni(2),:)),x,abs(v(posizioni(3),:)),'LineWidth',1.1)
        grid on; 
        
    case 'c'
        %
end 



%% definizione di funzioni 

function err = errore(u,v) 
    % u = soluzione analitica
    % v = metodo numerico 
    diff = (u - v).^2; 
    num = sum(diff,2); % il 2 è per sommare lungo le righe 
    den = sum(u.^2,2);
    err = sqrt(num./den);
end 
















