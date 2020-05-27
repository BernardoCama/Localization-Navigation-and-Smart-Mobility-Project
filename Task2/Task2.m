clear all, clc, close all

set(0,'DefaultTextFontSize',18)
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesFontSize',16)


%% scenario settings (4000x4000 m)
parameters.xmin = -2000; parameters.ymin = -2000;
parameters.xmax =  2000; parameters.ymax =  2000;

%% load measurements Task 2

load('Task2_trajectory_GR31.mat')

%% parameters
parameters.simulationTime = 200; %s
parameters.samplingTime = 1; %s
parameters.numberTrajectory = 100;

x = zeros (parameters.numberTrajectory,  parameters.simulationTime, 4);

for i = 1:parameters.numberTrajectory
    
    x(i, : ,:) = UEtrajectory{i};
    
end

% extract velocity
v = x(:,:,3:4);

%% plot UE trajectories
figure,hold on

for n=1:parameters.numberTrajectory
    
    hold on
    
    plot( x(n,:,1) , x(n,:,2))
    
    legend('UE')
    
    xlabel('[m]'), ylabel('[m]');
    
    title('Trajectories');
    
    xlim([parameters.xmin parameters.xmax])
    
    ylim([parameters.ymin parameters.ymax])
    
    axis equal
    
    grid on

end


%% plot statistics of first trajectory
figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,1));

plot([0:parameters.samplingTime:parameters.simulationTime-1],v(1,:,2));

title('Velocity');xlabel('time');ylabel('m/s')

acceleration = squeeze(diff(v(:,:,1:2),1,2));

figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,2));

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:,2));

title('Acceleration');xlabel('time');ylabel('m/s2')



sigma_velocity = squeeze(std(v(1,:,:),0,2))
mean_velocity = squeeze(mean(v(1,:,:),2))

sigma_acceleration = squeeze(std(acceleration(1,:,:),0,2))
mean_acceleration = squeeze(mean(acceleration(1,:,:),2))



%% motion model statistics
sigma_tot_velocity = mean(squeeze(std(v(:,:,:),0,2)),1)

sigma_tot_acceleration = mean(squeeze(std(acceleration(:,:,:),0,2)),1)
mean_tot_acceleration = mean(squeeze(mean(acceleration(:,:,:),2)),1)

%{
RIVEDERE
we can clearly see that is Motion Model M2 because of:
- zero sigma_tot_velocity
- constant mean_tot_velocity (about v_x = 1.387, v_y = 5.635)
- zero sigma_tot_acceleration
- zero mean_tot_acceleration

otherwise for M3 we would have obtained:
- zero sigma_tot_velocity
- zero mean_tot_velocity
- zero sigma_tot_acceleration
- zero mean_tot_acceleration

and for M4:
- zero sigma_tot_velocity
- zero mean_tot_velocity
- zero sigma_tot_acceleration
- mean_tot_acceleration != 0
%}

