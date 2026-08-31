clear all; 

% Run MG with event detection and load data
MGProjectileMotionFMAeventDetection;
test = load('MGProjectileMotionFMAeventDetection.1'); 

% Plot data
figure(1)
plot(test(:,1),test(:,2),'*-')
xlabel('x (m)')
ylabel('y (m)')
grid on