%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For paper, "A Gradient Descent Method for Finite Horizon Distributed Control of Discrete Time Systems" by S. Heinke and H. Werner.
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under GPLv3. See License.txt in the project root for license information.
% Author: Simon Heinke
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%% Figure 1
load('Data_sf.mat')
figure(1)
plot(100*x_,G2_,'r','LineWidth',1.5)
hold on
grid on
box on
plot(100*x_,G3_,'b','LineWidth',1.5)
plot(100*x_,G4_,'color',[0,0.5,0],'LineWidth',1.5)
plot(100*x_,G5_,'--','color',[0.3,0.3,0.3],'LineWidth',1.5)
legend('Alg. 1','Alg. 2','LMI','upper bound')
xlabel('Fraction of systems [%]')
ylabel('Centralized/ distribtuted')

%% Figure 2
load('Data_sf_Case2.mat')
figure(2)
plot(100*x_,G_,'b','LineWidth',1.5)
grid on
box on
legend('Alg. 2')
xlabel('Fraction of systems [%]')
ylabel('Centralized/ distribtuted')
ylim([0 1])

%% Figure 3
load('Data_sf_comp.mat')
figure(3)
plot(Ni(1:N_dlqr),t1,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_alg1),t3,'r','LineWidth',1.5)
plot(Ni,t2,'b','LineWidth',1.5)
plot(Ni(1:N_lmi),t4,'color',[0 0.5 0],'LineWidth',1.5)
legend('dlqr','Alg. 1','Alg. 2','LMI')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])

%% Figure 4
load('Data_of.mat')
figure(4)
plot(100*x_,P3_,'r','LineWidth',1.5)
hold on
grid on
box on
plot(100*x_,P1_,'b','LineWidth',1.5)
plot(100*x_,P4_,'r-.','LineWidth',1.5)
plot(100*x_,P2_,'b-.','LineWidth',1.5)
legend('Alg. 1_S','Alg. 2_S','Alg. 1_D','Alg. 2_D')
xlabel('Fraction of systems [%]')

ylabel('Centralized/ distribtuted')
ylim([0 1])

%% Figure 5
load('Data_of_comp.mat')
figure(5)
plot(Ni(1:N_h2syn),t1,'k','LineWidth',1.5)
hold on
grid on
plot(Ni(1:N_alg1_sof),t4,'r','LineWidth',1.5)
plot(Ni,t2,'b','LineWidth',1.5)
plot(Ni(1:N_alg1_dof),t5,'r-.','LineWidth',1.5)
plot(Ni,t3,'b-.','LineWidth',1.5)
legend('h2syn','Alg. 1_S','Alg. 2_S','Alg. 1_D','Alg. 2_D')
xlabel('Number of subsystems')
ylabel('Computation time [s]')
ylim([0 60])

