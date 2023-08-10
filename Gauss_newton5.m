%% Gaurav R Datta 201676662
%  University of Leeds
%  Dr James Mclaughlan and Dr Joshua Freeman
%  Gauss Newton Program that uses Jacobian Matrices method
%  of linear regression for 3-D localization

%% Reference paper: An Improved Trilateration Positioning Algorithm with Anchor ...
% Node Combination and K-Means Clustering...
% Qinghua Luo 1,2,3 , Kexin Yang 1...
% , Xiaozhen Yan 1,2,*, Jianfeng Li 1,2, Chenxu Wang 1,2 and Zhiquan Zhou
% 1,2 1 School of Information Science and Engineering, Harbin Institute of Technology, Weihai 264209, China...
% 2 Shandong Institute of Shipbuilding Technology, Ltd., Weihai 264209, China...
% 3 Shandong New Beiyang Information Technology Co., Ltd., Weihai 264209, China...
% * Correspondence: yanxiaozhen@hit.edu.cn

%% This code does not use filtering algorithm, nonetheless refers to the...
%%   technique of non linear gauss newton technique

clear
clc
close all
iter = 80; % maximum number of iterations
error = 1e-8; % error value
RSSE = 0;
%% Setting the Receiver anchor's points in cartesian plane
rx = [0;100;100;0]; % x
ry = [0;0;200;200];
rz = [0;0;0;0]; % z axis

%% Setting the true location of the Transmitter in cartesian plane

tx = [18;35;65];

%% Calculating the distance from the true Tx location to the anchor Rx r_i matrix
true_dist1 = sqrt((tx(1) - rx(1)).^2 + (tx(2) - ry(1)).^2 + (tx(3) - rz(1)).^2);
true_dist2 = sqrt((tx(1) - rx(2)).^2 + (tx(2) - ry(2)).^2 + (tx(3) - rz(2)).^2);
true_dist3 = sqrt((tx(1) - rx(3)).^2 + (tx(2) - ry(3)).^2 + (tx(3) - rz(3)).^2);
true_dist4 = sqrt((tx(1) - rx(4)).^2 + (tx(2) - ry(4)).^2 + (tx(3) - rz(4)).^2);
true_distance = [true_dist1;true_dist2;true_dist3;true_dist4]; % r_i

%% Setting the initial guess of the Transmitter location
guess = [21;53;72]; % x;y;z
tic
for i = 1:iter

    %% We calculate the distance between the guestimate position and the Rx anchors
    %% This keeps on continuosly updating every iteration
    guess_dist1 = sqrt((guess(1) - rx(1))^2 + (guess(2) - ry(1))^2 + (guess(3) - rz(1))^2);
    guess_dist2 = sqrt((guess(1) - rx(2))^2 + (guess(2) - ry(2))^2 + (guess(3) - rz(2))^2);
    guess_dist3 = sqrt((guess(1) - rx(3))^2 + (guess(2) - ry(3))^2 + (guess(3) - rz(3))^2);
    guess_dist4 = sqrt((guess(1) - rx(4))^2 + (guess(2) - ry(4))^2 + (guess(3) - rz(4))^2);
    guess_dist = [guess_dist1;guess_dist2;guess_dist3;guess_dist4];

    %% We calculate the Residual sum of errors phi function
    f1 = (true_distance(1) - guess_dist(1));
    f2 = (true_distance(2) - guess_dist(2));
    f3 = (true_distance(3) - guess_dist(3));
    f4 = (true_distance(4) - guess_dist(4));
    f_i = [f1;f2;f3;f4]; % stores all the difference between the true distance and the 
    % guesstimatef_i = true_dist - R_estimate

    %RSSE = f_i.^2; % Residual sum of squares
    
    %% We calculate the jacobian matrix of the f_i function above
    %  We calculate the partial derivates

    for k = 1:4 %rows
        for l = 1:3 %columns
            if l == 1
                jacob(k,l) = -2*(guess(1) - rx(k))/(2* sqrt((guess(1) - rx(k))^2 + (guess(2) - ry(k))^2 + (guess(3) - rz(k))^2));
            end
            if l == 2

                jacob(k,l) = -2*(guess(2) - ry(k))/(2* sqrt((guess(1) - rx(k))^2 + (guess(2) - ry(k))^2 + (guess(3) - rz(k))^2));
            end
            if l == 3

                jacob(k,l) = -2*(guess(3) - rz(k))/(2* sqrt((guess(1) - rx(k))^2 + (guess(2) - ry(k))^2 + (guess(3) - rz(k))^2));
            end
        end
    end
    
 %% Residual sum of errors is calculated iteratively
    for j = 1:4
        %residual(j,1) = true_distance(j) - sqrt((guess(1) - rx(j))^2  +  (guess(2) - ry(j))^2  +  (guess(3) - rz(j))^2); % residual value is calculated by...
        residual(j,1) = true_distance(j) - guess_dist(j);
        % subtracting distance of iterative guesses with the actual distance
        RSSE = RSSE + residual(j,1).^2;
    end

    jacob_trans = transpose(jacob);
    J_product = jacob_trans*jacob;
    J_f = jacob_trans*f_i;

    guess_new = guess - transpose(J_product\J_f);
    tolerance = abs(( guess - guess_new).^2);
    % wait until the tolerance reduces sufficiently below the error
    if (tolerance <= error)
        break
    end
    guess = guess_new;
end
toc
%% Plotting the data points in 3-D space
figure(1)
scatter3(rx,ry,rz,'filled','red')
hold on
scatter3(tx(1),tx(2),tx(3),'filled','blue')
hold on
scatter3(guess(1,1),guess(2,1),guess(3,1),'filled','green')
hold off
title('3-D localisation estimation using Gauss-Newton method','FontSize',14)
xlabel('X - axis in cm', 'FontSize',14)
ylabel('Y - axis in cm', 'FontSize',14)
zlabel('Z - axis in cm', 'FontSize',14)