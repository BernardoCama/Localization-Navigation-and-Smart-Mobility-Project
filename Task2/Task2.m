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

% mean over all trajectories of velocity vector
v_mean = [mean(x(:,:,3),1); mean(x(:,:,4),1)];

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


%% plot statistics
figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-1],v_mean(1,:));

plot([0:parameters.samplingTime:parameters.simulationTime-1],v_mean(2,:));

title('Velocity');xlabel('time');ylabel('m/s')

acceleration = diff(v_mean(1:2,:),1,2);

figure,hold on

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(1,:));

plot([0:parameters.samplingTime:parameters.simulationTime-2],acceleration(2,:));

title('Acceleration');xlabel('time');ylabel('m/s2')


% since, as we can see in the graphs, the values of velocities and
% accelerations are misleading after the 100_th sample (due to constraint
% of delimited space), we calculate the statistics on samples up to 100.

sigma_tot_velocity = std(v_mean(:,1:100),0,2)
mean_tot_velocity = mean(v_mean(:,1:100),2)

sigma_tot_acceleration = std(acceleration(:,1:100),0,2)
mean_tot_acceleration = mean(acceleration(:,1:100),2)

%{
we can clearly see that is Motion Model M2 because of:
- zero sigma_tot_velocity
- constant mean_tot_velocity (about v_x = 0.4, v_y = 0.6)
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




