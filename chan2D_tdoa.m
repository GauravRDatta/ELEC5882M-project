%% Gaurav Ramesh Datta 201676662
%  University of Leeds
% Dr James Mclaughlan and Dr Joshua Freeman
% Program to solve 2-D trilateration using Chan's algorithm

%% Description
% This source code solves linear algebra based on matrices ... 
% And helps in location of a transmitter based on Time Difference of
% Arrival based on Chan's algorithm.

%% Reference Paper
% IEEE TRANSACTIONS ON SIGNAL PROCESSING, VOL. 42, NO. 8, AUGUST 1994
% 1905 A Simple and Efficient Estimator for Hyperbolic Location ...
% Y. T. Chm, Senior Member, IEEE, and K. C. Ho, Member, IEEE
%%
clear
clc
close all
%% Declaring constants
c = 3*10^8; % c = speed of light

%% Setting the Receiver anchor's points in cartesian plane
rx = [0;100;100]; % x axis
ry = [0;0;200]; % y axis

%% Setting the true location of the Transmitter in cartesian plane
tx = [57;35];

%% Calculating the distance from the true Tx location to the anchor Rx
true_dist1 = sqrt((tx(1) - rx(1)).^2 + (tx(2) - ry(1)).^2);
true_dist2 = sqrt((tx(1) - rx(2)).^2 + (tx(2) - ry(2)).^2);
true_dist3 = sqrt((tx(1) - rx(3)).^2 + (tx(2) - ry(3)).^2);
true_distance = [true_dist1;true_dist2;true_dist3]; %distance vector

%% Time stamps
time = (true_distance./100)./c; % in seconds

%% Let anchor Rx0 be the reference 0th anchor
T01 = time(2) - time(1);
T02 = time(3) - time(1);
T0i = [T01;T02];

%% Calculating the distance between Rx1 anchor and other anchors using speed of light
% ri,1 = c x Timei,1
r_i1 = (c.*(T0i).*100); %ri1 = [r21;r31]

%% The K matrices for All Receiver anchors
%  K is measured from the origin

K = [rx(1)^2 + ry(1)^2;rx(2)^2 + ry(2)^2;rx(3)^2 + ry(3)^2];

%% Solving the Chan's algorithm
% the positions are given by - Matrices arithmetic involving inverse and
% squares

mat1 = -inv([rx(2) - rx(1), ry(2) - ry(1);rx(3) - rx(1), ry(3) - ry(1)]);
mat2 = [r_i1(1);r_i1(2)];                   %[r21;r31]*r1
mat3 = [r_i1(1)^2 - K(2) + K(1);r_i1(2)^2 - K(3) + K(1)]*0.5;

temp_solution = (mat1)*((mat2) + mat3);

%% solving for the temporary solution values using quadratic equation formula
% K1^2 = x1^2 + y1^2
% temp_solution(1,1) = x1; temp_solution(2,1) = y1

%% This section solves equation 9 of the referred paper.

L = mat1(1,1)*mat2(1,1) + mat1(1,2)*r_i1(2,1);
M = mat1(1,1)*mat3(1,1); 
N = mat1(1,2)*mat3(2,1);

P = mat1(2,1)*mat2(1,1) + mat1(2,2)*r_i1(2,1);
Q = mat1(2,1)*mat3(1,1);
S = mat1(2,2)*mat3(2,1);

%% Fitting the above values to quadratic equation to determine the distance
%  between the reference Rx anchor and transmitter
%  This section solves equation 8 of the referred paper

A = 1 - L^2 - P^2;
B = -2*(L*(M+N) + P*(Q+S) - rx(1,1)*L - ry(1,1)*P);
C =  -1*(K(1,1) + ((M+N)^2) + ((Q+S)^2) - (2*(rx(1,1)*(M+N))) - (2*(ry(1,1)*(Q+S))));

Discriminant = (B^2) - (4*A*C);
r1_1 = (-B + sqrt(Discriminant)) /(2*A);
r1_2 = (-B - sqrt(Discriminant)) /(2*A); % negative in this case

%% Puttin the r1_temp value positive value to the chan's algorithm
Accurate_pos = (mat1)*((mat2)*r1_1 + mat3);
disp(Accurate_pos)

%% PLotting figures
figure(1)
scatter(rx, ry,'filled','red')
hold on
scatter(tx(1),tx(2),'filled','green')
%scatter(rx(2), ry(2),'red','lineWidth',10,'MarkerSize',10)
hold on
scatter(Accurate_pos(1),Accurate_pos(2),'filled','green')
hold on
plot([rx(1) tx(1)],[ry(1) tx(2)],'bl--')
hold on
plot([rx(2) tx(1)],[ry(2) tx(2)],'g--')
hold on
plot([rx(3) tx(1)],[ry(3) tx(2)],'black--')
hold off
title('2-D Time Difference of Arrival','FontSize',18)
xlabel('x - axis distance in cm','FontSize',16)
ylabel('y - axis distance in cm','FontSize',16)

grid on
