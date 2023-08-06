%% Gaurav R Datta
% University of Leeds
% Dr James McLaughlan and Dr Joshua Freeman
% 3-D Trilateration algorithm for TDoA
% Cramer's rule used for 3-D trilateration
%%

clear 
clc
% All units in cm
c = 3*(10^8); % speed of light constant

%% Setting up of Receiver positions [x, y, z] coordinates
rx1 = [0,0,0];
rx2 = [100,0,25];
rx3 = [100,200,0];
rx4 = [0,200,50];

%% Setting up Transmitter position (given)
tx = [50,32,72];

%% Cramer's rule used to solve 3x3 linear system
% Step 1: Determinant is calculated
a1 = 2*(rx2(1) - rx1(1));
a2 = 2*(rx3(1) - rx1(1));
a3 = 2*(rx4(1) - rx1(1));

b1 = 2*(rx2(2) - rx1(2));
b2 = 2*(rx3(2) - rx1(2));
b3 = 2*(rx4(2) - rx1(2));

c1 = 2*(rx2(3) - rx1(3));
c2 = 2*(rx3(3) - rx1(3));
c3 = 2*(rx4(3) - rx1(3));

%% we calculate the distance between each anchor and the TX tag
d1 = sqrt((tx(1) - rx1(1)).^2 + (tx(2) - rx1(2)).^2 + (tx(3) - rx1(3)).^2);
d2 = sqrt((tx(1) - rx2(1)).^2 + (tx(2) - rx2(2)).^2 + (tx(3) - rx2(3)).^2);
d3 = sqrt((tx(1) - rx3(1)).^2 + (tx(2) - rx3(2)).^2 + (tx(3) - rx3(3)).^2);
d4 = sqrt((tx(1) - rx4(1)).^2 + (tx(2) - rx4(2)).^2 + (tx(3) - rx4(3)).^2);

%% We form the equations
%     d1^2 - d2^2 -  (x1^2  - x2^2)  - (y1^2    - y2^2)     -  (z1^2 - z2^2)
D1 = (d1^2 - d2^2) - (rx1(1)^2 - rx2(1)^2) - (rx1(2)^2 - rx2(2)^2) - (rx1(3)^2 - rx2(3)^2);

D2 = (d1^2 - d3^2) - (rx1(1)^2 - rx3(1)^2) - (rx1(2)^2 - rx3(2)^2) - (rx1(3)^2 - rx3(3)^2);

D3 = (d1^2 - d4^2) - (rx1(1)^2 - rx4(1)^2) - (rx1(2)^2 - rx4(2)^2) - (rx1(3)^2 - rx4(3)^2);

Denominator = det([a1 b1 c1; a2 b2 c2; a3 b3 c3]);
del_X = det([D1 b1 c1;D2 b2 c2; D3 b3 c3]);
x = abs(del_X/Denominator); % x - location of transmitter

del_Y = det([D1 a1 c1;D2 a2 c2; D3 a3 c3]);
y = abs(del_Y/Denominator); % y - location of transmitter

del_Z = det([D1 a1 b1; D2 a2 b2; D3 a3 b3]);
z = abs(del_Z/Denominator); % z - location of transmitter
tx_guess = [x;y;z];

figure(1)
scatter3(rx1(1),rx1(2),rx1(3),40,'filled',"blue")
hold on
scatter3(rx2(1),rx2(2),rx2(3),40,'filled',"blue")
hold on
scatter3(rx3(1),rx3(2),rx3(3),40,'filled',"blue")
hold on
scatter3(rx4(1),rx4(2),rx4(3),40,'filled',"blue")
hold on
scatter3(x,y,z,80,'filled',"red")
hold off
title('3-D space localisation Cramers rule','Fontsize',14)
xlabel('X-axis','FontSize',14)
ylabel('Y-axis','FontSize',14)
zlabel('Z-axis','FontSize',14)