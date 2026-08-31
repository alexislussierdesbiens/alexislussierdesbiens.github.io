clear all; 
clc; 

tcontroller = 0.05;
tend  = 5; 
t = 0:tcontroller:tend; 

% États initiaux: 
qA_ini                              =  0.5235987755982988;     % rad                 Initial Value
qB_ini                              =  1.299877806337475;      % rad                 Initial Value
qC_ini                              =  0.7945172960540738;     % rad                 Initial Value
qAp_ini                             =  0;                      % rad/sec             Initial Value
qBp_ini                             =  0;                      % rad/sec             Initial Value
qCp_ini                             =  0;                      % rad/sec             Initial Value
states = [qA_ini, qB_ini, qC_ini, qAp_ini, qBp_ini, qCp_ini]; 

% Trajectoire désirée: 
des_qA  = (30*pi/180)*sin(2*pi*1*t); 
des_qAp = (30*pi/180)*(2*pi*1)*cos(2*pi*1*t); 
des_qApp = -(30*pi/180)*(2*pi*1)^2*sin(2*pi*1*t); 

tgain = 1;
limit = 2000;
ActivateFF = 1; 
for i = 1:(length(t)-1)
    % Rename variables
    qA = states(i,1); qAp = states(i,4); 
    
    % Controller
    TA_FF = modelFourBarLinkage(states(i,:),des_qApp(i));
    
    TA(i) = -5000*tgain*(qA - des_qA(i)) - 1500*tgain*(qAp - des_qAp(i)) + ActivateFF*TA_FF;
    
    
    if TA(i)>limit  % Saturation
        TA(i) = limit;
    elseif TA(i)<-limit
        TA(i) = -limit;
    end
    
    % Simulate forward
    [tsim,VARsim,Outputsim] = MGFourBarKaneAugmented(states(i,:), TA(i), tcontroller);
    
    % Save states for display
    states = [states; VARsim]; 
end

%%
figure(2)
subplot(211)
plot(t,states(:,1)*180/pi,t,states(:,2)*180/pi,t,states(:,3)*180/pi,t,des_qA*180/pi)
xlabel('temps (s)')
ylabel('angles (deg)')
legend('qA','qB','qC','qA_des')
grid on
axis([0 tend -100 100])

subplot(212)
plot(t(1:end-1), TA)
xlabel('temps (s)')
ylabel('TA (Nm)')
grid on
axis([0 tend -4000 4000])

