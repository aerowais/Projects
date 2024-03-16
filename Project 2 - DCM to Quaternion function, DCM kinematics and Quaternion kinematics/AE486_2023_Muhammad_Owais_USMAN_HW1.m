clear all
%% Q2
% Generate three example rotation angle sets

% Asymmetric
yaw1   = -pi + (pi-(-pi)) .* rand(1,1); 
roll1  = -pi + (pi-(-pi)) .* rand(1,1);
pitch1 = -pi/2 + (pi/2-(-pi/2)) .* rand(1,1); 

% Symmetric
yaw2   = -pi + (pi-(-pi)) .* rand(1,1);
roll2  = -pi + (pi-(-pi)) .* rand(1,1);
pitch2 = 0 + (pi-0).* rand(1,1);              

% Asymmetric
yaw3   = -pi + (pi-(-pi)) .* rand(1,1);
roll3  = -pi + (pi-(-pi)) .* rand(1,1);
pitch3 = -pi/2 + (pi/2-(-pi/2)) .* rand(1,1); 


% Calculate the DCM matrices for the above rotation angle sets
dcm1 = angle2dcm( yaw1, pitch1, roll1)
dcm2 = angle2dcm( yaw2, pitch2, roll2) 
dcm3 = angle2dcm( yaw3, pitch3, roll3)

% Calculating Quaternion using the first method:
qua_method1 = horzcat(method1(dcm1),method1(dcm2), method1(dcm3))

% Checking the above quaternion using the inbuilt function "dcm2quat"
qua_correct = horzcat(transpose( dcm2quat(dcm1)), transpose( dcm2quat(dcm2)),transpose( dcm2quat(dcm3)))

% Calculating Quaternion using the second method:
qua_method2 = horzcat(method2(dcm1),method2(dcm2), method2(dcm3))

%% Q3
load("hw1_data.mat")

% Defining constants and initialising
t = 1;

delta_t = 1;

A_bi = eye(3);

quater = method2(A_bi);


q_4(1,1) = quater(1,1);
q_1(1,1) = quater(2,1);
q_2(1,1) = quater(3,1);
q_3(1,1) = quater(4,1);

% Just a variable created to make sure that the resulting magnititude of the quaternions is 1

mag (t, 1) = sqrt( (q_4(t,1)^2) + (q_1(t,1)^2) + (q_2(t,1)^2) + (q_3(t,1)^2) );

% The loop for propagating the DCM.

for t  = 1:25000
    %Classical Runge Kutta 4th order is used and K1, K2, K3, K4 are found
    
    w_t = W_act(:,t);
   
    w_th = W_act(:,t); 
    w_t1 = W_act(:,(t));
    
    k_1 = -cpm(w_t) * A_bi;
    
    A_bi_s = A_bi + k_1*0.5*delta_t;
    
    % the angular velocity is assumed to be discreet, hence angular velocity after half a time step is still equal to the angular velocity at the current time step  
    k_2 = -cpm(w_th) * A_bi_s;
    
    A_bi_ss = A_bi + k_2*0.5*delta_t;
    
    k_3 = -cpm(w_th) * A_bi_ss;
    
    A_bi_sss = A_bi + k_3*delta_t;
    
    k_4 = -cpm(w_t1) * A_bi_sss;
    
    A_bi = A_bi + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    
    % method2 is used to avoid singularities rising from q4 being 0.
    quater = method2 (A_bi);
    
    q_1(t+1,1) = quater(1,1);
    q_2(t+1,1) = quater(2,1);
    q_3(t+1,1) = quater(3,1);
    q_4(t+1,1) = quater(4,1);
 
    
    mag (t+1, 1) = sqrt( (q_4(t+1,1)^2) + (q_1(t+1,1)^2) + (q_2(t+1,1)^2) + (q_3(t+1,1)^2) );
    
    % The value of time is incresed by 1
    t = t + 1;
end
% Time vector is made

time = [0:25000]

% Plotting

figure (1)
hold on
tiledlayout(2,2);

% Tile 1
nexttile
plot(time, q_1, 'b')
title('q_1 vs time')
xlabel('time [s]')
ylabel('q_1')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, q_2, 'b')
title('q_2 vs time')
xlabel('time [s]')
ylabel('q_2')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, q_3, 'b')
title('q_3 vs time')
xlabel('time [s]')
ylabel('q_3')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 4
nexttile
plot(time, q_4, 'b')
title('q_4 vs time')
xlabel('time [s]')
ylabel('q_4')
xlim([0 25000])
ylim([-0.2 1.2])

hold off

figure (2)
hold on


plot(time, mag, 'b')
title('magnitude vs time')
xlabel('time [s]')
ylabel('magnitude')
xlim([0 25000])
ylim([-0.2 1.2])
hold off

%% Q4 -- Propagating using Quaternion Kinematics
load("hw1_data.mat")
% Defining constants and initialising
t = 1

delta_t = 1

A_bi = eye(3)

quater = method1(A_bi)


q_1(1,1) = quater(1,1);
q_2(1,1) = quater(2,1);
q_3(1,1) = quater(3,1);
q_4(1,1) = quater(4,1);

% Just a variable created to make sure that the resulting magnititude of the quaternions is 1

mag (t, 1) = sqrt( (q_4(t,1)^2) + (q_1(t,1)^2) + (q_2(t,1)^2) + (q_3(t,1)^2) )

% The loop for propagating the DCM.

for t  = 1:25000
    %Classical Runge Kutta 4th order is used and K1, K2, K3, K4 are found
    
    w_t = W_act(:,t);
   
    w_th = W_act(:,t); 
    w_t1 = W_act(:,t);
    
    k_1 = 0.5*omega(w_t) * quater;
    
    quater_s = quater + k_1*0.5*delta_t;
    
    % the angular velocity is assumed to be discreet, hence angular velocity after half a time step is still equal to the angular velocity at the current time step  
    k_2 = 0.5*omega(w_th) * quater_s;
    
    quater_ss = quater + k_2*0.5*delta_t;
    
    k_3 = 0.5*omega(w_th) * quater_ss;
    
    quater_sss = quater + k_3*delta_t;
    
    k_4 = 0.5*omega(w_t1) * quater_sss;
    
    quater = quater + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    
    % method2 is used to avoid singularities rising from q4 being 0.
    
    mag (t+1, 1) = sqrt( (quater(1,1)^2) + (quater(2,1)^2) + (quater(3,1)^2) + (quater(4,1)^2) );
    
    q_1(t+1,1) = quater(1,1)/mag (t+1, 1);
    q_2(t+1,1) = quater(2,1)/mag (t+1, 1);
    q_3(t+1,1) = quater(3,1)/mag (t+1, 1);
    q_4(t+1,1) = quater(4,1)/mag (t+1, 1);
    
    % The value of time is incresed by 1
    t = t + 1;
end

% Time vector is made

time = [0:25000];

% Plotting

figure (3)
hold on
tiledlayout(2,2);

% Tile 1
nexttile
plot(time, q_1, 'r')
title('q_1 vs time')
xlabel('time [s]')
ylabel('q_1')
xlim([0 25000])
ylim([-1.2 1.2])


% Tile 2
nexttile
plot(time, q_2, 'r')
title('q_2 vs time')
xlabel('time [s]')
ylabel('q_2')
xlim([0 25000])
ylim([-1.2 1.2])


% Tile 3
nexttile
plot(time, q_3, 'r')
title('q_3 vs time')
xlabel('time [s]')
ylabel('q_3')
xlim([0 25000])
ylim([-1.2 1.2])


% Tile 4
nexttile
plot(time, q_4, 'r')
title('q_4 vs time')
xlabel('time [s]')
ylabel('q_4')
xlim([0 25000])
ylim([-1.2 1.2])



hold off

figure (4)
hold on


plot(time, mag, 'r')
title('magnitude vs time')
xlabel('time [s]')
ylabel('magnitude')
xlim([0 25000])
ylim([-1.2 1.2])
hold off
%% Q4 -- Applying constraint on DCM

% Defining constants and initialising
t = 1;

delta_t = 1;

A_bi = eye(3);

quater = method2(A_bi);


q_4(1,1) = quater(1,1);
q_1(1,1) = quater(2,1);
q_2(1,1) = quater(3,1);
q_3(1,1) = quater(4,1);

% Just a variable created to make sure that the resulting magnititude of the quaternions is 1

mag (t, 1) = sqrt( (q_4(t,1)^2) + (q_1(t,1)^2) + (q_2(t,1)^2) + (q_3(t,1)^2) );

% The loop for propagating the DCM.

for t  = 1:25000
    %Classical Runge Kutta 4th order is used and K1, K2, K3, K4 are found
    
    w_t = W_act(:,t);
   
    w_th = W_act(:,t); 
    w_t1 = W_act(:,(t));
    
    k_1 = -cpm(w_t) * A_bi;
    
    A_bi_s = A_bi + k_1*0.5*delta_t;
    
    % the angular velocity is assumed to be discreet, hence angular velocity after half a time step is still equal to the angular velocity at the current time step  
    k_2 = -cpm(w_th) * A_bi_s;
    
    A_bi_ss = A_bi + k_2*0.5*delta_t;
    
    k_3 = -cpm(w_th) * A_bi_ss;
    
    A_bi_sss = A_bi + k_3*delta_t;
    
    k_4 = -cpm(w_t1) * A_bi_sss;
    
    A_bi = A_bi + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4);
    
    % Orthonormalizing the DCM matrix using the Gramschmidt method
    A_bi = gramschmidt(A_bi);
    
    % method2 is used to avoid singularities rising from q4 being 0.
    quater = method2 (A_bi);
    
    q_1(t+1,1) = quater(1,1);
    q_2(t+1,1) = quater(2,1);
    q_3(t+1,1) = quater(3,1);
    q_4(t+1,1) = quater(4,1);
 
    
    mag (t+1, 1) = sqrt( (q_4(t+1,1)^2) + (q_1(t+1,1)^2) + (q_2(t+1,1)^2) + (q_3(t+1,1)^2) );
    
    % The value of time is incresed by 1
    t = t + 1;
end
% Time vector is made

time = [0:25000];

% Plotting

figure (5)
hold on
tiledlayout(2,2);

% Tile 1
nexttile
plot(time, q_1, 'g')
title('q_1 vs time')
xlabel('time [s]')
ylabel('q_1')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, q_2, 'g')
title('q_2 vs time')
xlabel('time [s]')
ylabel('q_2')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, q_3, 'g')
title('q_3 vs time')
xlabel('time [s]')
ylabel('q_3')
xlim([0 25000])
ylim([-1.2 1.2])

% Tile 4
nexttile
plot(time, q_4, 'g')
title('q_4 vs time')
xlabel('time [s]')
ylabel('q_4')
xlim([0 25000])
ylim([-0.2 1.2])

hold off

figure (6)
hold on


plot(time, mag, 'g')
title('magnitude vs time')
xlabel('time [s]')
ylabel('magnitude')
xlim([0 25000])
ylim([-1.2 1.2])

hold off












%% Functions
% First method to convert DCM to Quaternion
function qua = method1(dcm)
    % Calculating the 4th term of the Quaternion from the DCM matrix
    q4 = 0.5* sqrt(1 + trace(dcm));

    % Calculating rest of the terms of the Quaternion from the 4th term of the Quaternion

    q1 = (1/(4*q4))*( dcm(2,3) - dcm(3,2));
    q2 = (1/(4*q4))*( dcm(3,1) - dcm(1,3));
    q3 = (1/(4*q4))*( dcm(1,2) - dcm(2,1));

    % Finding the quaternion vector in the format of the function dcm2quat
    % arranges the elements in, that is the 4th element becomes the first one
    % instead of the 4th one

    qua = [ q1; q2; q3; q4];

    % Checking the solution code using the matlab function dcm2quat
end
% ------------

% Second method to convert DCM to Quaternion
function qua = method2 (dcm)

    % First we calculate the Euler axis/angle from the DCM
    theta = acos( ( dcm(1,1)+dcm(2,2)+dcm(3,3)-1 )/2);
    
    e1 = ( dcm(2,3)-dcm(3,2) )/ (2* sin(theta));
    e2 = ( dcm(3,1)-dcm(1,3) )/ (2* sin(theta));
    e3 = ( dcm(1,2)-dcm(2,1) )/ (2* sin(theta));

    % Calculating the Quaternion from the Euler axis/angle
    
    q1 = e1*sin(theta/2);
    q2 = e2*sin(theta/2);
    q3 = e3*sin(theta/2);
    q4 = cos(theta/2);
    
    % Finding the quaternion vector in the format of the function dcm2quat
    % arranges the elements in, that is the 4th element becomes the first one
    % instead of the 4th one

    qua = [ q1; q2; q3; q4];
    
end


% Function to find the cross product matrix

function cr = cpm (x)

    cr(1,1) = 0;
    cr(2,2) = 0;
    cr(3,3) = 0;
    
    cr(1,2) = -x(3,1);
    cr(1,3) = x(2,1);
    
    cr(2,1) = x(3,1);
    cr(2,3) = -x(1,1);
    
    cr(3,1) = -x(2,1);
    cr(3,2) = x(1,1);


end
% Function to find the omega matrix

function cr = omega (x)

    cr(1,1) = 0;
    cr(2,2) = 0;
    cr(3,3) = 0;
    cr(4,4) = 0;
    
    cr(1,2) = x(3,1);
    cr(1,3) = -x(2,1);
    cr(1,4) = x(1,1);
    
    cr(2,1) = -x(3,1);
    cr(2,3) = x(1,1);
    cr(2,4) = x(2,1);
    
    cr(3,1) = x(2,1);
    cr(3,2) = -x(1,1);
    cr(3,4) = x(3,1);
    
    cr(4,1) = -x(1,1);
    cr(4,2) = -x(2,1);
    cr(4,3) = -x(3,1);


end

% Function for orthonormalization of the DCM

function U = gramschmidt(V)
    U = zeros(3,3);
    U(:,1) = V(:,1) / norm(V(:,1));
    for i = 2:3
        U(:,i) = V(:,i);
        for j = 1:i-1
            U(:,i) = U(:,i) - (U(:,j)'*U(:,i)) * U(:,j);
        end
        U(:,i) = U(:,i) / norm(U(:,i));
    end
end
