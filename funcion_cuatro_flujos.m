function [R, T, A, ocupacion, saltos_intercapa, tiempo, scat] = funcion_cuatro_flujos;
% FUNCION_MONOCROMATICA Cálculo para una longitud de onda.
%   INPUTS:
%
%     NT: número total de fotones entrantes en el sistema.
%
%     PK: array con las probabilidades de absorción de todas las capas,
%     incluido el sustrato.
%
%     PS: array con las probabilidades de scattering hacia atrás de
%     todas las capas, incluido el sustrato.
%
%     PA: array con las probabilidades de scattering hacia delante de
%     todas las capas, incluido el sustrato.
%
%     PR: array con las probabilidades de rebote normal de
%     todas las capas, incluido el sustrato.
%
%   OUTPUTS:
%
%     R: array lógico de dimensión NT con valor 'true' en aquellas
%     posiciones que corresponden a fotones reflejados.
%
%     T: array lógico de dimensión NT con valor 'true' en aquellas
%     posiciones que corresponden a fotones transmitidos.
%
%     A: array de dimensión NT que almacena el número de capa donde fue
%     absorbido cada fotón ('0' en caso de no serlo).
%
%     ocupacion: array de la misma dimensión que PK, PS y PS_atras que
%     registra el paso de los fotones por cada capa.
%
%     saltos_intercapa: array de dimensión NT que registra el número de
%     saltos intercapa de cada fotón.
%     
%     tiempo: array de dimensión NT que cuenta el tiempo que pasa el fotón
%     dentro del material hasta que es absorbido, reflejado o transmitido.
%
%     scat: array de dimensiones NT que cuenta el número de veces que el
%     fotón sufre dispersión por scattering.


NT=100000;
N=1000;

%CAPAS DE TRANSICIÓN
Cap1=40*N/100;
Cap2=60*N/100;


medio2=(Cap1+1:Cap2-1);


% Inicialización
R = false(NT,1);
T = false(NT,1);
A = zeros(NT,1);

ocupacion = zeros(N,1);
saltos_intercapa = -ones(NT,1); % lo inicializamos con '-1' para considerar
                                % que los fotones reflejados en la primera 
                                % capa dan '0' saltos
%tiempo=zeros(NT,1);                              
tiempo=(N*0.1).*randn(NT,1);
tiempo=tiempo-min(tiempo);

scat=zeros(NT,1);
%scat(1:N/2)=1;

PK=linspace(0,0,N)'; 
PS=linspace(0,0,N)';
PA=linspace(0,0,N)';
PR=linspace(0,0,N)';   


for I=1:NT
   
    
        if scat(I)>0 %PROBABILIDADES para FOTONES DIFUNDIDOS
    
            %probabilidades MEDIO 1
            PK(1:N)=0.004/N;
            PS(1:N)=0.04/N;
            PA(1:N)=0.04/N;
            PR(1:N)=0;
            %probabilidades MEDIO 2
            PK(medio2)=0.007/N;
            PS(medio2)=0.07/N;
            PA(medio2)=0.07/N;
            PR(medio2)=0;
            
            
        else %PROBABILIDADES para FOTONES COLIMADOS

            %probabilidades MEDIO 1
            PK(1:N)=0.002/N;
            PS(1:N)=0.02/N;
            PA(1:N)=0.02/N;
            PR(1:N)=0;
            %probabilidades MEDIO 2
            PK(medio2)=0.005/N;
            PS(medio2)=0.05/N;
            PA(medio2)=0.05/N;
            PR(medio2)=0;
        end
    
    
    
    % Cada fotón empieza por la primera capa moviéndose hacia la derecha
    J = 1;
    moviendose = true;
    moviendose_hacia_la_derecha = true;
    
    
    
    
    while moviendose % mientras el fotón siga moviéndose
        ocupacion(J) = ocupacion(J) + 1; % registramos el paso del fotón 
                                         % por la capa J
        saltos_intercapa(I) = saltos_intercapa(I) + 1; % registramos el 
                                                       % salto del fotón a
                                                       % otra capa    
  
  


        
     %PROBABILDIADES DE LAS CAPAS DE TRANSICIÓN SEGÚN EL SENTIDO                                                  
        if  moviendose_hacia_la_derecha 
            if scat(I)>0 %son difusos
                PR(Cap1)=0.42; 
                PR(Cap2)=1;
                
            else         %son colimados
                PR(Cap1)=0.04; 
                PR(Cap2)=0.04;
                
            end
        else
            moviendose_hacia_la_derecha=false; %se mueve hacia la izq
            if scat(I)>0 %son difusos
                PR(Cap1)=1; 
                PR(Cap2)=0.42;
                
            else         %son colimados
                PR(Cap1)=0.04; 
                PR(Cap2)=0.04;
                
            end
        end
                                                
                                                       
                                                       
        rho = rand; % número aleatorio entre 0 y 1
        if rho<PK(J) % hay absorción
            moviendose = false; % el fotón deja de moverse
            A(I) = J; % el fotón 'I' es absorbido en la capa 'J'
            
        elseif rho<(PK(J)+PS(J)) % hay scattering hacia atrás
            
           scat(I)=scat(I)+1;
           
           if scat(I)>0
                tiempo(I)=tiempo(I)+2;
           else
                tiempo(I)=tiempo(I)+1;
           end
            
            if moviendose_hacia_la_derecha                    
                if isequal(J,1) % el fotón es reflejado en la primera capa
                    moviendose = false; % el fotón deja de moverse
                    R(I) = true; % el fotón 'I' es reflejado
                   
                else
                    moviendose_hacia_la_derecha = false; % invertimos el sentido
                    J = J-1; % retrocedemos una capa
                end                                                                
            else
                moviendose_hacia_la_derecha = true; % invertimos el sentido
                J = J+1; % avanzamos una capa
                
            end    
            
        elseif rho<(PK(J)+PS(J)+PA(J)) %HAY SCATTERING HACIA DELANTE
            
            scat(I)=scat(I)+1;
            
            if scat(I)>0
                tiempo(I)=tiempo(I)+2;
            else
                tiempo(I)=tiempo(I)+1;
            end
            
            if moviendose_hacia_la_derecha
                if isequal(J,N) %sufre scattering hacia adelante en la última capa
                    moviendose=false;
                    T(I)=true;
                else
                    J=J+1;
                end
            else 
                if isequal(J,1) % el fotón es reflejado en la primera capa
                    moviendose = false; % el fotón deja de moverse
                    R(I) = true; % el fotón 'I' es reflejado
                else
                    J = J-1;
                end
            end
            
        elseif rho<(PK(J)+PS(J)+PA(J)+PR(J)) %REBOTE NORMAL
            
            if scat(I)>0
                tiempo(I)=tiempo(I)+2;
            else
                tiempo(I)=tiempo(I)+1;
            end
            
            if moviendose_hacia_la_derecha                    
                if isequal(J,1) % el fotón es reflejado en la primera capa
                    moviendose = false; % el fotón deja de moverse
                    R(I) = true; % el fotón 'I' es reflejado
                   
                else
                    moviendose_hacia_la_derecha = false; % invertimos el sentido
                    J = J-1; % retrocedemos una capa
                    
                end                                                                
            else
                moviendose_hacia_la_derecha = true; % invertimos el sentido
                J = J+1; % avanzamos una capa
                
            end    
           
            
        else % fotón atraviesa    
            
            if scat(I)>0
                tiempo(I)=tiempo(I)+2;
            else
                tiempo(I)=tiempo(I)+1;
            end
            
            if moviendose_hacia_la_derecha                
                if isequal(J,N) % el fotón es transmitido en la última capa
                    moviendose = false; % el fotón deja de moverse
                    T(I) = true; % el fotón 'I' es transmitido
                else
                    J = J+1; % avanzamos una capa                    
                end
            else                
                if isequal(J,1) % el fotón es reflejado en la primera capa
                    moviendose = false; % el fotón deja de moverse
                    R(I) = true; % el fotón 'I' es reflejado
                else
                    J = J-1; % retrocedemos una capa
                end
            end
        end        
    end
 
end

%COMANDOS:
%length(R(R==true))
%nnz(A) número de veces que A es distinto de cero
%beerlambert(exponencial con sólo absorción)
%histogram(tiempo(R==true))
%[i]=find(A) posiciones en las que A no es cero
%nonzeros(A) valores cuando A es distinto de cero
%{
CONVERGENCIA
RT=(((1-4)*sinh(sqrt(99))+sqrt(99)*0.4*cosh((sqrt(99))))/((10-0.4)*sinh(sqrt(99))+sqrt(99)*cosh((sqrt(99)))))
RT=0.0501
Inc=abs((RT-RMC))*100
%}
%{
LEY BEER-LAMBERT (solo absorcion)
N=(1:500)';
Nf=NT*exp(-PK.*N);
hold on
plot(N,Nf*100/NT,'bx')
plot(N,ocupacion*100/NT,'rx')
%}
%{
LEY STOKES
Rg=0.4;
x=length(R(R==true))/NT; Sin sustrato x=0.2256
y=length(T(T==true))/NT; Sin sustrato y=0.3929
Rf=x+(((y^2)*Rg)/(10000-x*Rg));
%}
%{
HISTOGRAMAS
[R, T, A, ocupacion, saltos_intercapa, tiempo] = funcion_monocromatica(NT,PK,PS,PA,PR);
subplot(211), histogram(tiempo(R==true))
subplot(212), histogram(tiempo(T==true))
%}