function TAmodel = modelFourBarLinkage(states,qApp)

% Renames the states from function inputs: 
qA                              =  states(1);     % rad                 
qB                              =  states(2);     % rad                 
qC                              =  states(3);     % rad                 
qAp                             =  states(4);     % rad/sec             
qBp                             =  states(5);     % rad/sec              
qCp                             =  states(6);     % rad/sec     

% Constants: 
b                               =  50;                     % N*m/rad             Constant
g                               =  9.81;                   % m/sec^2             Constant
H                               =  200;                    % N                   Constant
LA                              =  1;                      % m                   Constant
LB                              =  2;                      % m                   Constant
LC                              =  2;                      % m                   Constant
LN                              =  1;                      % m                   Constant
mA                              =  10;                     % kg                  Constant
mB                              =  20;                     % kg                  Constant
mC                              =  20;                     % kg                  Constant

IA = 0.08333333333333333*mA*LA^2;
IB = 0.08333333333333333*mB*LB^2;
IC = 0.08333333333333333*mC*LC^2;

% Feedforward equation: 
TAmodel = mB*g*LA*sin(qA) + 0.5*mA*g*LA*sin(qA) + b*qAp + 0.5*mB*LA*LB*sin(qA-qB)*qBp^2 + 0.25*(mA*LA^2+4*mB*LA^2+4*IA)*qApp + 0.5*mB*LA*cos(qA-qB)*(LC*qCp^2-sin(qC)*(LA*sin(qA)*qAp^2+LB*sin(qB)*qBp^2-2*LA*cos(qA)*qApp)-cos(qC)*(LA*cos(  ...
qA)*qAp^2+LB*cos(qB)*qBp^2+2*LA*sin(qA)*qApp))/sin(qB-qC) + 0.25*LA*sin(qA)*(cos(qB)*(2*LC*(2*H*cos(qC)-mC*g*sin(qC))+(mC*LC^2+4*IC)*(LB*qBp^2+sin(qB)*(LA*sin(qA)*qAp^2-LC*sin(qC)*qCp^2-2*LA*cos(qA)*qApp)-cos(qB)*(LC*cos(qC)*qCp^2-LA*cos(  ...
qA)*qAp^2-LA*sin(qA)*qApp))/(LC*sin(qB-qC)))/(LC*sin(qB-qC))-cos(qC)*(2*mB*LB*(g*sin(qB)-LA*sin(qA-qB)*qAp^2)+(mB*LB^2+4*IB)*(LC*qCp^2-cos(qC)*(LA*cos(qA)*qAp^2+LB*cos(qB)*qBp^2+LA*sin(qA)*qApp)-sin(qC)*(LA*sin(qA)*qAp^2+LB*sin(qB)*qBp^  ...
2-2*LA*cos(qA)*qApp))/(LB*sin(qB-qC)))/(LB*sin(qB-qC))) - 0.25*LA*cos(qA)*(sin(qB)*(2*LC*(2*H*cos(qC)-mC*g*sin(qC))+(mC*LC^2+4*IC)*(LB*qBp^2+cos(qB)*(LA*cos(qA)*qAp^2-LC*cos(qC)*qCp^2)+sin(qB)*(LA*sin(qA)*qAp^2-LC*sin(qC)*qCp^2-LA*cos(  ...
qA)*qApp))/(LC*sin(qB-qC)))/(LC*sin(qB-qC))-sin(qC)*(2*mB*LB*(g*sin(qB)-LA*sin(qA-qB)*qAp^2)+(mB*LB^2+4*IB)*(LC*qCp^2-cos(qC)*(LA*cos(qA)*qAp^2+LB*cos(qB)*qBp^2)-sin(qC)*(LA*sin(qA)*qAp^2+LB*sin(qB)*qBp^2-LA*cos(qA)*qApp))/(LB*sin(qB-qC)))/(  ...
LB*sin(qB-qC)));
