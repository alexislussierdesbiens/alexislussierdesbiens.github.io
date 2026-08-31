function [t,VAR,Output] = robotballWOctl(TA, TB, states, dt)
%===========================================================================
% File: robotballWOctl.m created on Wed Jul 22 2015 by MotionGenesis 5.5.
% Advanced Research Licensee: AlexisLussierDesbiens (until February 2016).
% Portions copyright (c) 1988-2014 Paul Mitiguy and 2009-2014 Motion Genesis.
% Paid-up MotionGenesis Advanced Research licensees are granted the right
% right to distribute this code for legal academic (non-professional) purposes only,
% provided this copyright notice appears in all copies and distributions.
%===========================================================================
% The software is provided "as is", without warranty of any kind, express or    
% implied, including but not limited to the warranties of merchantability or    
% fitness for a particular purpose. In no event shall the authors, contributors,
% or copyright holders be liable for any claim, damages or other liability,     
% whether in an action of contract, tort, or otherwise, arising from, out of, or
% in connection with the software or the use or other dealings in the software. 
%===========================================================================
%  Rotations
%  Translations
%  Forces/Moments
%  Constraints
% --------------------------------------------------------------------
%        Desired motion for qA and qB
%===========================================================================
eventDetectedByIntegratorTerminate1OrContinue0 = [];
%%%% Commenter cette ligne
% TA=0; TB=0;
wx=0; wz=0; qCp=0; wyp=0; qApp=0; qBpp=0; qCpp=0; CNy=0; ENx=0; ENz=0; I=0;


%-------------------------------+--------------------------+-------------------+-----------------
% Quantity                      | Value                    | Units             | Description
%-------------------------------|--------------------------|-------------------|-----------------
h                               =  0.6;                    % m                   Constant
LB                              =  0.8;                    % m                   Constant
LC                              =  1.2;                    % m                   Constant
m                               =  800;                    % kg                  Constant
mE                              =  50;                     % kg                  Constant
r                               =  0.4;                    % m                   Constant

%%%% Commenter les variables d'état initiales (vient des entrées de la
%%%% fonction)
qA                              =  states(1);                      % rad                 Initial Value
qB                              =  states(2);                      % rad                 Initial Value
qC                              =  states(3);                      % rad                 Initial Value
wy                              =  states(4);                      % rad/sec             Initial Value
qAp                             =  states(5);                      % rad/sec             Initial Value
qBp                             =  states(6);                      % rad/sec             Initial Value

tInitial                        =  0.0;                    % second              Initial Time
%%%% Calculer seulement un pas de temps discret
tFinal                          =  dt;                   % second              Final Time
integStp                        =  0.01;                    % second              Integration Step
%%% Ne pas afficher a l'ecran/fichier
printIntScreen                  =  0;                      % 0 or +integer       0 is NO screen output
printIntFile                    =  0;                      % 0 or +integer       0 is NO file   output

absError                        =  1.0E-05;                %                     Absolute Error
relError                        =  1.0E-08;                %                     Relative Error
%-------------------------------+--------------------------+-------------------+-----------------

% Unit conversions
DEGtoRAD = pi / 180.0;
RADtoDEG = 180.0 / pi;

% Evaluate constants
I = 0.4*mE*r^2;


VAR = SetMatrixFromNamedQuantities;
[t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, integStp, absError, relError, VAR, printIntScreen, printIntFile );
OutputToScreenOrFile( [], 0, 0 );   % Close output files



%===========================================================================
function sys = mdlDerivatives( t, VAR, uSimulink )
%===========================================================================
SetNamedQuantitiesFromMatrix( VAR );
qCp = LB*cos(qB)*qBp/(LC*cos(qC));
wx = -(LB*cos(qB)+LC*cos(qC))*qAp/r;
wz = (LB*sin(qB)*qBp+LC*sin(qC)*qCp)/r;


%%%% Commenter ces lignes aussi
% TA = 0;
% TB = 0;
COEF = zeros(4,4);
COEF(1,1) = I*(1+(LB*cos(qB)+LC*cos(qC))^2/r^2) + m*(LB^2*cos(qB)^2+LC^2*cos(qC)^2+2*LB*LC*cos(qB)*cos(qC)) + mE*(LB^2*cos(qB)^2+LC^2*cos(qC)^2+2*LB*LC*cos(qB)*cos(qC));
COEF(1,3) = I;
COEF(2,2) = LB^2*(m*(1-cos(qB)*cos(qB+qC)/cos(qC))+mE*(1-cos(qB)*cos(qB+qC)/cos(qC))+I*sin(qB)*(sin(qB)+cos(qB)*tan(qC))/r^2);
COEF(2,4) = -LB*LC*(m*(cos(qB+qC)-cos(qB)/cos(qC))+mE*(cos(qB+qC)-cos(qB)/cos(qC))-I*sin(qC)*(sin(qB)+cos(qB)*tan(qC))/r^2);
COEF(3,1) = 1;
COEF(3,3) = 1;
COEF(4,2) = LB*cos(qB);
COEF(4,4) = -LC*cos(qC);
RHS = zeros(1,4);
RHS(1) = TA + qAp*(I*(LB*cos(qB)+LC*cos(qC))*(wz+(LB*sin(qB)*qBp+LC*sin(qC)*qCp)/r)/r+2*m*(LB*LC*sin(qB)*cos(qC)*qBp+LB^2*sin(qB)*cos(qB)*qBp+LB*LC*sin(qC)*cos(qB)*qCp+LC^2*sin(qC)*cos(qC)*qCp)+2*mE*(LB*LC*sin(qB)*cos(qC)*qBp+LB^2*sin(  ...
qB)*cos(qB)*qBp+LB*LC*sin(qC)*cos(qB)*qCp+LC^2*sin(qC)*cos(qC)*qCp));
RHS(2) = TB + LB*(I*(sin(qB)+cos(qB)*tan(qC))*(wx*qAp-(LB*cos(qB)*qBp^2+LC*cos(qC)*qCp^2)/r)/r-m*(sin(qB)*(LB*cos(qB)+LC*cos(qC))*qAp^2+LB*cos(qB)*sin(qB+qC)*qBp^2/cos(qC)+cos(qB)*tan(qC)*(LB*cos(qB)+LC*cos(qC))*qAp^2+LC*sin(qB+qC)*qCp^  ...
2)-mE*(sin(qB)*(LB*cos(qB)+LC*cos(qC))*qAp^2+LB*cos(qB)*sin(qB+qC)*qBp^2/cos(qC)+cos(qB)*tan(qC)*(LB*cos(qB)+LC*cos(qC))*qAp^2+LC*sin(qB+qC)*qCp^2));
RHS(4) = LB*sin(qB)*qBp^2 - LC*sin(qC)*qCp^2;
SolutionToAlgebraicEquations = COEF \ transpose(RHS);

% Update variables after uncoupling equations
qApp = SolutionToAlgebraicEquations(1);
qBpp = SolutionToAlgebraicEquations(2);
wyp = SolutionToAlgebraicEquations(3);
qCpp = SolutionToAlgebraicEquations(4);

sys = transpose( SetMatrixOfDerivativesPriorToIntegrationStep );
end



%===========================================================================
function VAR = SetMatrixFromNamedQuantities
%===========================================================================
VAR = zeros(1,6);
VAR(1) = qA;
VAR(2) = qB;
VAR(3) = qC;
VAR(4) = wy;
VAR(5) = qAp;
VAR(6) = qBp;
end


%===========================================================================
function SetNamedQuantitiesFromMatrix( VAR )
%===========================================================================
qA = VAR(1);
qB = VAR(2);
qC = VAR(3);
wy = VAR(4);
qAp = VAR(5);
qBp = VAR(6);
end


%===========================================================================
function VARp = SetMatrixOfDerivativesPriorToIntegrationStep
%===========================================================================
VARp = zeros(1,6);
VARp(1) = qAp;
VARp(2) = qBp;
VARp(3) = qCp;
VARp(4) = wyp;
VARp(5) = qApp;
VARp(6) = qBpp;
end



%===========================================================================
function Output = mdlOutputs( t, VAR, uSimulink )
%===========================================================================
ENx = cos(qA)*(LB*cos(qB)+LC*cos(qC));
ENz = -sin(qA)*(LB*cos(qB)+LC*cos(qC));
CNy = h + LB*sin(qB) - LC*sin(qC);

Output = zeros(1,9);
Output(1) = t;
Output(2) = ENx;
Output(3) = ENz;

Output(4) = qA*RADtoDEG;
Output(5) = qAp*RADtoDEG;
Output(6) = qB*RADtoDEG;
Output(7) = qB*RADtoDEG;
Output(8) = qC*RADtoDEG;
Output(9) = CNy;
end


%===========================================================================
function OutputToScreenOrFile( Output, shouldPrintToScreen, shouldPrintToFile )
%===========================================================================
persistent FileIdentifier hasHeaderInformationBeenWritten;

if( isempty(Output) ),
   if( ~isempty(FileIdentifier) ),
      for( i = 1 : 2 ),  fclose( FileIdentifier(i) );  end
      clear FileIdentifier;
      fprintf( 1, '\n Output is in the files robotballWOctl.i  (i=1,2)\n' );
      fprintf( 1, '\n Note: Plots are automatically generated by issuing the OutputPlot command in MotionGenesis\n' );
      fprintf( 1, '\n To load and plot columns 1 and 2 with a solid line and columns 1 and 3 with a dashed line, enter:\n' );
      fprintf( 1, '    someName = load( ''robotballWOctl.1'' );\n' );
      fprintf( 1, '    plot( someName(:,1), someName(:,2), ''-'', someName(:,1), someName(:,3), ''--'' )\n\n' );
   end
   clear hasHeaderInformationBeenWritten;
   return;
end

if( isempty(hasHeaderInformationBeenWritten) ),
   if( shouldPrintToScreen ),
      fprintf( 1,                '%%       t             Enx            Enz\n' );
      fprintf( 1,                '%%     (sec)           (m)            (m)\n\n' );
   end
   if( shouldPrintToFile && isempty(FileIdentifier) ),
      FileIdentifier = zeros(1,2);
      FileIdentifier(1) = fopen('robotballWOctl.1', 'wt');   if( FileIdentifier(1) == -1 ), error('Error: unable to open file robotballWOctl.1'); end
      fprintf(FileIdentifier(1), '%% FILE: robotballWOctl.1\n%%\n' );
      fprintf(FileIdentifier(1), '%%       t             Enx            Enz\n' );
      fprintf(FileIdentifier(1), '%%     (sec)           (m)            (m)\n\n' );
      FileIdentifier(2) = fopen('robotballWOctl.2', 'wt');   if( FileIdentifier(2) == -1 ), error('Error: unable to open file robotballWOctl.2'); end
      fprintf(FileIdentifier(2), '%% FILE: robotballWOctl.2\n%%\n' );
      fprintf(FileIdentifier(2), '%%      qA             qA''            qB             qB             qC             CNy\n' );
      fprintf(FileIdentifier(2), '%%     (deg)        (deg/sec)        (deg)        (deg/sec)        (deg)           (m)\n\n' );
   end
   hasHeaderInformationBeenWritten = 1;
end

if( shouldPrintToScreen ), WriteNumericalData( 1,                 Output(1:3) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(1), Output(1:3) );  end
if( shouldPrintToFile ),   WriteNumericalData( FileIdentifier(2), Output(4:9) );  end
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
function [t,VAR,Output] = IntegrateForwardOrBackward( tInitial, tFinal, integStp, absError, relError, VAR, printIntScreen, printIntFile )
%===========================================================================
OdeMatlabOptions = odeset( 'RelTol',relError, 'AbsTol',absError, 'MaxStep',integStp, 'Events',@EventDetection );
t = tInitial;                 epsilonT = 0.001*integStp;                tFinalMinusEpsilonT = tFinal - epsilonT;
printCounterScreen = 0;       integrateForward = tFinal >= tInitial;    tAtEndOfIntegrationStep = t + integStp;
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
      tAtEndOfIntegrationStep = t + integStp;
   else
      t = timeEventOccurredInIntegrationStep;
      VAR = nStatesArraysAtEvent;
      printCounterScreen = 0;
      printCounterFile   = 0;
      [isIntegrationFinished, VAR] = EventDetectedByIntegrator( t, VAR, nIndexOfEvents );
   end
end
end


%==============================================
end   % End of embedded function robotballWOctl
%==============================================
