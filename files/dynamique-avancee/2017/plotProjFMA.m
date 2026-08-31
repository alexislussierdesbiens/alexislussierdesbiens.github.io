clear all; 

% Plot data from MG
test = load('MGProjectileMotionFMA.1'); 

% Plot data
figure(2)
plot(test(:,2),test(:,3),'*-')
xlabel('x (m)')
ylabel('y (m)')
grid on


figure(1)
plot(test(:,1),test(:,3),'*-')
xlabel('t (sec)')
ylabel('y (m)')
grid on