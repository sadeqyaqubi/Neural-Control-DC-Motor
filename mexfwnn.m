clc; clear all; close all;

k1 = 0.04;
k2 = 0.04;
l1 = .200;
l2 = .200;
r1 = 2;
r2 = 2;
j1 = 12e-5;
j2 = 12e-5;
T = 0.02;
tf = 10;

tic
mex control.cpp;
toc
soluMat = control(k1, k2, j1, j2, l1, l2, r1, r2, T, tf);

figure; plot(soluMat(:,1), soluMat(:,2)); xlabel('Time [sec]'); ylabel('X_1 [rad]');
figure; plot(soluMat(:,1), soluMat(:,3)); xlabel('Time [sec]'); ylabel('X_2 [rad/s]');
figure; plot(soluMat(:,1), soluMat(:,4));  xlabel('Time [sec]'); ylabel('X_3 [A]');
figure; plot(soluMat(:,1), soluMat(:,5));  xlabel('Time [sec]'); ylabel('X_4 [A]');

figure; plot(soluMat(:,1), soluMat(:,6));  xlabel('Time [sec]'); ylabel('U_1 [V]');
figure; plot(soluMat(:,1), soluMat(:,7));  xlabel('Time [sec]'); ylabel('U_2 [V]');


figure; plot(soluMat(:,1), soluMat(:,2)); hold all; plot(soluMat(:,1), soluMat(:,8)); xlabel('Time [sec]'); ylabel('X_1 [rad]'); legend('System Response', 'FWNN Estimation');
figure; plot(soluMat(:,1), soluMat(:,3)); hold all; plot(soluMat(:,1), soluMat(:,9)); xlabel('Time [sec]'); ylabel('X_2 [rad/s]'); legend('System Response', 'FWNN Estimation');
figure; plot(soluMat(:,1), soluMat(:,4)); hold all; plot(soluMat(:,1), soluMat(:,10)); xlabel('Time [sec]'); ylabel('X_3 [rad/s]'); legend('System Response', 'FWNN Estimation');
figure; plot(soluMat(:,1), soluMat(:,5)); hold all; plot(soluMat(:,1), soluMat(:,11)); xlabel('Time [sec]'); ylabel('X_3 [rad/s]'); legend('System Response', 'FWNN Estimation');
