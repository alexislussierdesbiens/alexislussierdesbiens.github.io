function [TA, TB] = computeFFTorque(qA, qB, qAp, qBp, qApp, qBpp, wy, qC, dt)

% Constantes du modele
h                               =  0.6;                    % m                   Constant
LB                              =  0.8;                    % m                   Constant
LC                              =  1.2;                    % m                   Constant
m                               =  700;                    % kg                  Constant
mE                              =  70;                     % kg                  Constant
r                               =  0.4;                    % m                   Constant
I = 0.4*mE*r^2;

% Les équations proviennent de MG

% Calcul de qCp, qCpp, wx et wz selon la trajectoire désirée de qA et qB... 
qCp     = LB.*cos(qB).*qBp./(LC.*cos(qC));
qCpp    = -(LB.*sin(qB).*qBp.^2-LC.*sin(qC).*qCp.^2-LB.*cos(qB).*qBpp)./(LC*cos(qC));
wx = -(LB*cos(qB)+LC*cos(qC)).*qAp/r;
wz = (LB*sin(qB).*qBp+LC*sin(qC).*qCp)/r;

% Calcul des couples requis selon les équations dynamiques
TA      = I.*wy' + (I.*(1+(LB.*cos(qB)+LC.*cos(qC)).^2./r.^2)+m.*(LB.^2.*cos(qB).^2+LC.^2.*cos(qC).^2+2.*LB.*LC.*cos(qB).*cos(qC))+mE.*(LB.^2.*cos(qB).^2+LC.^2.*cos(qC).^2+2.*LB.*LC.*cos(qB).*cos(qC))).*qApp - qAp.*(I.*(LB.*cos(qB)+LC.*cos(qC)).*(wz+(LB.*sin(qB).*qBp+LC.*sin(qC).*qCp)./r)./r+2.*m.*(LB.*LC.*sin(qB).*cos(qC).*qBp+LB.^2.*sin(qB).*cos(qB).*qBp+LB.*LC.*sin(qC).*cos(qB).*qCp+LC.^2.*sin(qC).*cos(qC).*qCp)+2.*mE.*(LB.*LC.*sin(qB).*cos(qC).*qBp+LB.^2.*sin(qB).*cos(qB).*qBp+LB.*LC.*sin(qC).*cos(qB).*qCp+LC.^2.*sin(qC).*cos(qC).*qCp));
TB      = LB.*(LB.*(m.*(1-cos(qB).*cos(qB+qC)./cos(qC))+mE.*(1-cos(qB).*cos(qB+qC)./cos(qC))+I.*sin(qB).*(sin(qB)+cos(qB).*tan(qC))./r.^2).*qBpp-I.*(sin(qB)+cos(qB).*tan(qC)).*(wx.*qAp-(LB.*cos(qB).*qBp.^2+LC.*cos(qC).*qCp.^2+LC.*sin(qC).*qCpp)./r)./r-m.*(LC.*cos(qB+qC).*qCpp-sin(qB).*(LB.*cos(qB)+LC.*cos(qC)).*qAp.^2-cos(qB).*tan(qC).*(LB.*cos(qB)+LC.*cos(qC)).*qAp.^2-LC.*sin(qB+qC).*qCp.^2-cos(qB).*(LB.*sin(qB+qC).*qBp.^2+LC.*qCpp)./cos(qC))-mE.*(LC.*cos(qB+qC).*qCpp-sin(qB).*(LB.*cos(qB)+LC.*cos(qC)).*qAp.^2-cos(qB).*tan(qC).*(LB.*cos(qB)+LC.*cos(qC)).*qAp.^2-LC.*sin(qB+qC).*qCp.^2-cos(qB).*(LB.*sin(qB+qC).*qBp.^2+LC.*qCpp)./cos(qC))); 
