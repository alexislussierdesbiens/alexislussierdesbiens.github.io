clear all; clc; 

% Garder dt à 0.01 sec pour le contrôleur
dt  = 0.01; 
t = 0:dt:16; 

% Choisir cercle (Stitch=0) ou pattern de zigzag (Stitch=1)
Stitch = 1;

% Choisir un controleur (1 à 4) 
SelectDemo = 1; 

LowGainsWithoutFF   = 1; 
LowGainsWithFF      = 2; 
HighGainsWithoutFF  = 3; 
HighGainsWithFF     = 4; 

switch SelectDemo
    % Gains ajustés pour dt = 0.01 sec
    case LowGainsWithoutFF
        kpA = 1000; kdA = 1000; 
        kpB = 1000; kdB = 1000;
        EnableFF = 0; 
    case LowGainsWithFF
        kpA = 1000; kdA = 1000; 
        kpB = 1000; kdB = 1000;
        EnableFF = 1; 
    case HighGainsWithoutFF
        kpA = 10000; kdA = 10000; 
        kpB = 10000; kdB = 10000;
        EnableFF = 0; 
    case HighGainsWithFF
        kpA = 10000; kdA = 10000; 
        kpB = 10000; kdB = 10000;
        EnableFF = 1;
    otherwise
        error('Cas non supporté')
end

% Valeurs initiales 
qA = 0*pi/180;      
qB = 40*pi/180;     
qC = 0.2649472990121713;  % ATTENTION: dépends de qA et qB selon les équations de contraintes
wy = 0;
qAp = 0;
qBp = 0;
VAR = [qA; qB; qC; wy; qAp; qBp]; 

% Trajectoire désirée (cercle ou cercle + zigzag)
% incluant les dérivées (important pour FF)
%   - Idéalement  on calculerait également qC, qC' et qC'' ici... 
qB0 = 40*pi/180;
tF  = 16; 
qADesired = 2*pi*t/tF; 
qBDesired = qB0*ones(size(t)) + Stitch*10*pi/180*sin(2*pi*t/tF*100);
qADesiredp     = [0, diff(qADesired)]/dt;
qADesiredpp    = [0, diff(qADesiredp)]/dt;
qBDesiredp     = [0, diff(qBDesired)]/dt;
qBDesiredpp    = [0, diff(qBDesiredp)]/dt;

% Valeurs initiales de l'erreur du contrôleur 
qAError = 0;
qBError = 0;

% Calculs sur un pas de temps à répéter sur toute la durée de la simulation... 
for i = 1:length(t)
    % Renommer les variables en fonction des états précédents
    qA  = VAR(1); 
    qB  = VAR(2); 
    qC = VAR(3);
    wy  = VAR(4);

    % Sauvegarder l'erreur précédente pour calcul de la dérivée de l'erreur
    qAErrorOld = qAError; 
    qBErrorOld = qBError; 
    
    % Calcul de l'erreur
    qAError = qA - qADesired(i);
    qBError = qB - qBDesired(i);

    % Calcul de la commande selon contrôleur PD
    TApid = -kpA*qAError - kdA*(qAError - qAErrorOld)/dt; 
    TBpid = -kpB*qBError - kdB*(qBError - qBErrorOld)/dt; 
    
    % Calcul de la commande selon contrôleur Feedforward
    [TAff, TBff] = ComputeFFTorque(qADesired(i), qBDesired(i), ...
        qADesiredp(i), qBDesiredp(i), ...
        qADesiredpp(i), qBDesiredpp(i), ...
        wy, qC, dt); 
    
    % Combinaison des commandes PD et FF
    TA = TApid+EnableFF*TAff; 
    TB = TBpid+EnableFF*TBff; 
    
    % Calcul de la dynamique pour un pas de temps
    [tint,VAR,Output] = robotballmodified(TA, TB, VAR, dt);
    SavedOutput(i,:) = Output;
end

%% Affichage des résultats
Enx = SavedOutput(:,2);
Eny = SavedOutput(:,3);
qA = SavedOutput(:,4);
qB = SavedOutput(:,6);

figure(1)

plot(Enx, Eny, 'b-')
axis equal; grid on
xlabel('x (m)'); ylabel('y (m)')

figure(2)

subplot(311)
plot(t, qA', t, qADesired*180/pi)
xlabel('t (sec)'); ylabel('qA (deg)')
legend('Mesured','Desired')
grid on

subplot(312)
plot(t, qB', t, qBDesired*180/pi)
xlabel('t (sec)'); ylabel('qB (deg)')
legend('Mesured','Desired')
grid on


subplot(313)
plot(t, qADesired*180/pi -  qA', ...
    t, qBDesired*180/pi - qB')
xlabel('t (sec)'); ylabel('erreur (deg)')
legend('qA','qB')
grid on
