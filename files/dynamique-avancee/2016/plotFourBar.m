function plotFourBar 

data = load('MGFourBarKaneAugmented.1');
pos = load('MGFourBarKaneAugmented.2');

t = data(:,1);
qA = data(:,2); 
qB = data(:,3); 
qC = data(:,4); 

xpos = pos(:,1:2:end);
ypos = pos(:,2:2:end);

% figure(1)
% skipNFrames = 5; 
% plot(ypos(1:skipNFrames:end,:)', xpos(1:skipNFrames:end,:)','-'); 
% ylabel('nx')
% xlabel('ny')
% set(gca,'YDir','reverse');

skipNFrames = 5; 
for i = 1:skipNFrames:length(xpos)
    figure(1)
    plot(ypos(i,:)', xpos(i,:)','-o','LineWidth',3); 
    ylabel('nx')
    xlabel('ny')
    set(gca,'YDir','reverse');
    axis equal
    axis([-2.5 2.5 -2.5 2.5])
    grid on
    pause(0.02)
end

figure(2)
plot(t,qA, t, qB, t, qC)
xlabel('temps (s)')
ylabel('angles (deg)')
legend('qA','qB','qC')
grid on
