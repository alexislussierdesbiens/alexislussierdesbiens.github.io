function [t,VAR,Output] = MGFourBarKaneAugmented(states, TA, tcontroller)
%===========================================================================
% File: MGFourBarKaneAugmented.m created on Wed Jul 27 2016 by MotionGenesis 5.7.
% Advanced Student Licensee: AlexisLussierDesbiens (until May 2018).
% Portions copyright (c) 2009-2015 Motion Genesis LLC.  Rights reserved.
% Paid-up MotionGenesis Advanced Student licensees are granted the right
% right to distribute this code for legal student-academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
%  Masses et inerties
%  Rotations
%  Translations
%  Constraintes
%  EDM selon Kane
%  Solution à l'équilibre statique
%  Pour la résolution dynamique, on aimerait avoir des équation de contraintes linéaires.
%  Pour ce faire, on différentie les équations de contraintes et on trouve les conditions initiales.
%    - Une stabilisation des contraintes pourrait être requise…
%  Trouver les conditions initiales (qB, qC, qB', qC') pour des valeurs de qA et qA'
% % Solution Dynamique
%  Methode #1
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
% TA=0;
FCx=0; FCy=0; qApp=0; qBpp=0; qCpp=0; Abx=0; Aby=0; Bcx=0; Bcy=0; IA=0; IB=0; IC=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
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

qA                              =  states(1);     % rad                 Initial Value
qB                              =  states(2);     % rad                 Initial Value
qC                              =  states(3);     % rad                 Initial Value
qAp                             =  states(4);     % rad/sec             Initial Value
qBp                             =  states(5);     % rad/sec             Initial Value
qCp                             =  states(6);     % rad/sec             Initial Value

tInitial                        =  0.0;                    % second              Initial Time
tFinal                          =  tcontroller;            % sec                 Final Time
tStep                           =  tcontroller;            % second              Integration Step
printIntScreen                  =  0;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  0;                      % 0 or +integer       0 is NO file   output
absError                        =  1.0E-07;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;

% Evaluate constants
IA = 0.08333333333333333*mA*LA^2;
IB = 0.08333333333333333*mB*LB^2;
IC = 0.08333333333333333*mC*LC^2;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files



%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );

% Quantities to be specified (not assigned in MotionGenesis).
% TA = 0;
COEF = zeros(5,5);
COEF(1,1) = mB*LA^2 + 0.25*mA*LA^2 + IA;
COEF(1,2) = 0.5*mB*LA*LB*cos(qA-qB);
COEF(1,4) = LA*sin(qA);
COEF(1,5) = -LA*cos(qA);
COEF(2,1) = 0.5*mB*LA*LB*cos(qA-qB);
COEF(2,2) = 0.25*mB*LB^2 + IB;
COEF(2,4) = LB*sin(qB);
COEF(2,5) = -LB*cos(qB);
COEF(3,3) = 0.25*mC*LC^2 + IC;
COEF(3,4) = -LC*sin(qC);
COEF(3,5) = LC*cos(qC);
COEF(4,1) = -LA*cos(qA);
COEF(4,2) = -LB*cos(qB);
COEF(4,3) = LC*cos(qC);
COEF(5,1) = LA*sin(qA);
COEF(5,2) = LB*sin(qB);
COEF(5,3) = -LC*sin(qC);
RHS = zeros(1,5);
RHS(1) = TA - mB*g*LA*sin(qA) - 0.5*mA*g*LA*sin(qA) - b*qAp - 0.5*mB*LA*LB*sin(qA-qB)*qBp^2;
RHS(2) = -0.5*mB*LB*(g*sin(qB)-LA*sin(qA-qB)*qAp^2);
RHS(3) = 0.5*LC*(2*H*cos(qC)-mC*g*sin(qC));
RHS(4) = LC*sin(qC)*qCp^2 - LA*sin(qA)*qAp^2 - LB*sin(qB)*qBp^2;
RHS(5) = LC*cos(qC)*qCp^2 - LA*cos(qA)*qAp^2 - LB*cos(qB)*qBp^2;
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
qApp = SolutionToAlgebraicEquations(1);
qBpp = SolutionToAlgebraicEquations(2);
qCpp = SolutionToAlgebraicEquations(3);
FCx = SolutionToAlgebraicEquations(4);
FCy = SolutionToAlgebraicEquations(5);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros(1,6);
VAR(1) = qA;
VAR(2) = qB;
VAR(3) = qC;
VAR(4) = qAp;
VAR(5) = qBp;
VAR(6) = qCp;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qA = VAR(1);
qB = VAR(2);
qC = VAR(3);
qAp = VAR(4);
qBp = VAR(5);
qCp = VAR(6);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros(1,6);
VARp(1) = qAp;
VARp(2) = qBp;
VARp(3) = qCp;
VARp(4) = qApp;
VARp(5) = qBpp;
VARp(6) = qCpp;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
Abx = LA*cos(qA);
Aby = LA*sin(qA);
Bcx = LA*cos(qA) + LB*cos(qB);
Bcy = LA*sin(qA) + LB*sin(qB);

Output = zeros(1,17);
Output(1) = t;
Output(2) = qA;
Output(3) = qB;
Output(4) = qC;
Output(5) = qAp;
Output(6) = qBp;
Output(7) = qCp;
Output(8) = FCx;
Output(9) = FCy;

Output(10) = 0.0;
Output(11) = 0.0;
Output(12) = Abx;
Output(13) = Aby;
Output(14) = Bcx;
Output(15) = Bcy;
Output(16) = 0.0;
Output(17) = LN;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files MGFourBarKaneAugmented.i  (i=1,2)\n' );
      fprintf( 1, '\n Note: Plots are automatically generated by issuing the OutputPlot command in MotionGenesis\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''MGFourBarKaneAugmented.1'' );\n' );
      fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             qA             qB             qC             qA''            qB''            qC''            FCx            FCy\n' );
      fprintf( 1,                '%%     (sec)          (rad)          (rad)          (rad)        (rad/sec)      (rad/sec)      (rad/sec)      (Newtons)      (Newtons)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros(1,2);
      FileIdentifier(1) = fopen('MGFourBarKaneAugmented.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file MGFourBarKaneAugmented.1'); end
      fprintf(FileIdentifier(1), '%% FILE: MGFourBarKaneAugmented.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             qA             qB             qC             qA''            qB''            qC''            FCx            FCy\n' );
      fprintf(FileIdentifier(1), '%%     (sec)          (rad)          (rad)          (rad)        (rad/sec)      (rad/sec)      (rad/sec)      (Newtons)      (Newtons)\n\n' );
      FileIdentifier(2) = fopen('MGFourBarKaneAugmented.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file MGFourBarKaneAugmented.2'); end
      fprintf(FileIdentifier(2), '%% FILE: MGFourBarKaneAugmented.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%      Nox            Noy            Abx            Aby            Bcx            Bcy            Cnx            Cny\n' );
      fprintf(FileIdentifier(2), '%%    (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)        (UNITS)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:9) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(10:17) );  end
end


%===========================================================================
function WriteNumericalData( fileIdentifier, Output )
%===========================================================================
numberOfOutputQuantities = length( Output );
if( numberOfOutputQuantities > 0 ),
   for( i = 1 : numberOfOutputQuantities ),
      fprintf( fileIdentifier, ' %- 14.6E', Output(i) );
   end
   fprintf( fileIdentifier, '\n' );
end
end



%===========================================================================
function [functionsToEvaluateForEvent, eventTerminatesIntegration1Otherwise0ToContinue, eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1] = EventDetection( t, VAR, uSimulink )
%===========================================================================
% Detects when designated functions are zero or cross zero with positive or negative slope.
% Step 1: Uncomment call to mdlDerivatives and mdlOutputs.
% Step 2: Change functionsToEvaluateForEvent,                      e.g., change  []  to  [t - 5.67]  to stop at t = 5.67.
% Step 3: Change eventTerminatesIntegration1Otherwise0ToContinue,  e.g., change  []  to  [1]  to stop integrating.
% Step 4: Change eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1,  e.g., change  []  to  [1].
% Step 5: Possibly modify function EventDetectedByIntegrator (if eventTerminatesIntegration1Otherwise0ToContinue is 0).
%---------------------------------------------------------------------------
% mdlDerivatives( t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
% mdlOutputs(     t, VAR, uSimulink );        % UNCOMMENT FOR EVENT HANDLING
functionsToEvaluateForEvent = [];
eventTerminatesIntegration1Otherwise0ToContinue = [];
eventDirection_AscendingIs1_CrossingIs0_DescendingIsNegative1 = [];
eventDetectedByIntegratorTerminate1OrContinue0 = eventTerminatesIntegration1Otherwise0ToContinue;
end


%===========================================================================
function [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents )
%===========================================================================
isIntegrationFinished = eventDetectedByIntegratorTerminate1OrContinue0( nIndexOfEvents );
if( ~isIntegrationFinished ),
   SetNamedQuantitiesFromMatrix( VAR );
%  Put code here to modify how integration continues.
   VAR = SetMatrixFromNamedQuantities;
end
end



%===========================================================================
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, tStep, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',tStep, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*tStep;                   tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + tStep;
printCounterFile   = 0;       isIntegrationFinished = 0;
mdlDerivatives( t, VAR, 0 );
while 1,
   if( (integrateForward && t >= tFinalMinusEpsilonT) || (~integrateForward && t <= tFinalMinusEpsilonT) ), isIntegrationFinished = 1;  end
   shouldPrintToScreen = printIntScreen && ( isIntegrationFinished || printCounterScreen <= 0.01 );
   shouldPrintToFile   = printIntFile   && ( isIntegrationFinished || printCounterFile   <= 0.01 );
   if( isIntegrationFinished || shouldPrintToScreen || shouldPrintToFile ),
      Output = mdlOutputs( t, VAR, 0 );
      OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile );
      if( isIntegrationFinished ), break;  end
      if( shouldPrintToScreen ), printCounterScreen = printIntScreen;  end
      if( shouldPrintToFile ),   printCounterFile   = printIntFile;    end
   end
   [TimeOdeArray, VarOdeArray, timeEventOccurredInIntegrationStep, nStatesArraysAtEvent, nIndexOfEvents] = ode45( @mdlDerivatives, [t tAtEndOfIntegrationStep], VAR, OdeMatlabOptions, 0 );
   if( isempty(timeEventOccurredInIntegrationStep) ),
      t = TimeOdeArray( length(TimeOdeArray) );
      VAR = VarOdeArray( length(TimeOdeArray), : );
      printCounterScreen = printCounterScreen - 1;
      printCounterFile   = printCounterFile   - 1;
      if( abs(tAtEndOfIntegrationStep - t) >= abs(epsilonT) ), warning('numerical integration failed'); break;  end
      tAtEndOfIntegrationStep = t + tStep;
   else
      t = timeEventOccurredInIntegrationStep;
      VAR = nStatesArraysAtEvent;
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents );
   end
end
end


%===============================================
end    % End of function MGFourBarKaneAugmented
%===============================================
