clear all

load ("hw2_data");

% Q1

z_b = [ 0; 0; 1 ] ; % initial vector in the body frame which is the unit z vector

z_i = sun_eci(:,1);


a = z_i(1,1);
b = z_i(2,1);
c = z_i(3,1);

psi = atan2( a,(-b) );

phi = atan2( ( a*sin(psi) - b*cos(psi) ), (c) ); 

yaw1 = psi;
roll = phi;
yaw2 = 0;

A_bi = angle2dcm( yaw1, roll, yaw2, 'ZXZ' );

% check = A_bi*z_i
% -------------------------------------------------------------------------
% Part 2
omega_bi = [ 0 ; 0 ; (2*pi)/(15) ];

J = [0.018 0 0;
     0 0.018 0;
     0 0 0.065];

check = cross(omega_bi, (J*omega_bi) ); % Cross product is zero vector hence omega_bi is constant


prompt = "Enter the Question Number: ";  % Enter the number of the desired question:
                                          % Q1                                   == 1
                                          % Q2, with gravity gradient            == 21
                                          % Q2, with Magnetic disturbance torque ==22
                                          
                                        



condition = input(prompt);
% Initializing

t = 1;

delta_t = 0.1;

mu = 3.986004418e+14;



quater = dcm2quat(A_bi);

quater = transpose(quater);

q_1(1,1) = quater(2,1);
q_2(1,1) = quater(3,1);
q_3(1,1) = quater(4,1);
q_4(1,1) = quater(1,1);

quater(1,1) = q_1(1,1);
quater(2,1) = q_2(1,1);
quater(3,1) = q_3(1,1);
quater(4,1) = q_4(1,1);



r_ix = pos_eci(1,:);
r_iy = pos_eci(2,:);
r_iz = pos_eci(3,:);

iqi = zeros(1, 4);

z_eci_x(1,1) = z_i(1,1);
z_eci_y(1,1) = z_i(2,1);
z_eci_z(1,1) = z_i(3,1);

Time = [0:24999];% Just a variable created to make sure that the resulting magnititude of the quaternions is 1

M2 = [-0.09; 0.01; 0.11];

% mag (t, 1) = sqrt( (q_4(t,1)^2) + (q_1(t,1)^2) + (q_2(t,1)^2) + (q_3(t,1)^2) )

% The loop for propagating the DCM.
mi = mag_eci(:,1:25000);
ri = pos_eci(:,1:25000);

r_ii(:,1) = ri(:,1);
m_ii(:,1) = mi(:,1);

Ti = 0.1;
j = 1;

ti    = linspace(0,25000,25000/delta_t);

% for interpolation
for i = 2:250000

    if ( (i-j)/j == 10)

        j = j + 1;

        r_ii(:,i) = ri(:,j);
        m_ii(:,i) = mi(:,j);

    else

        r_ii(:,i) = ri(:,j) + (Ti-ti(1,j))*(ri(:,j+1)-ri(:,j));
        m_ii(:,i) = mi(:,j) + (Ti-ti(1,j))*(mi(:,j+1)-mi(:,j));

    end

    Ti = Ti + delta_t;

end


for t  = 0:249999
    %Classical Runge Kutta 4th order is used and K1, K2, K3, K4 are found
    
    if condition == 1                                   % for no external torques
            N = 0;
    elseif condition == 21                              % for gravity gradient torque
%         r_ix = interp1(Time,pos_eci(1,:),t*0.1);
%         r_iy = interp1(Time,pos_eci(2,:),t*0.1);
%         r_iz = interp1(Time,pos_eci(3,:),t*0.1);
%         
%         r_i = [r_ix ; r_iy; r_iz];
        
        nadir = A_bi * ( (-r_ii(:,t+1))/( norm( r_ii(:,t+1) ) ) );
        n_gg = ( (3 * mu)/(  (norm( r_ii(:,t+1)  ))^3 )) * cross(nadir , J*nadir);
               
        N = n_gg;  
        
        % RK-4 Method to find the value of omega_bi at this point
    
        % omega_bi_dot = inv(J) * ( N - cross( omega_bi, J * omega_bi) );
    
        ok1 = omega_bi_dot(J, N, omega_bi);
     
        omega_bi_s = omega_bi + 0.5*ok1*delta_t;
    
        ok2 = omega_bi_dot(J, N, omega_bi_s);
    
        omega_bi_ss = omega_bi + 0.5*ok2*delta_t;
    
        ok3 = omega_bi_dot(J, N, omega_bi_ss);
    
        omega_bi_sss = omega_bi + ok3*delta_t;
    
        ok4 = omega_bi_dot(J, N, omega_bi_sss);
    
        omega_bi = omega_bi + (1/6) * ( ok1 + 2*ok2 + 2*ok3 + ok4 )*delta_t;
    
    elseif condition == 22                                  % for magnetic disturbance torque
%         b_ix = interp1(Time,mag_eci(1,:),t*0.1);
%         b_iy = interp1(Time,mag_eci(2,:),t*0.1);
%         b_iz = interp1(Time,mag_eci(3,:),t*0.1);

        
        %b_i = [b_ix ; b_iy; b_iz];
  
        b_b = A_bi * m_ii(:,t+1);
        n_md = cross( M2, b_b);
               
        N = n_md * (10^-9);  
        
        % RK-4 Method to find the value of omega_bi at this point
    
        % omega_bi_dot = inv(J) * ( N - cross( omega_bi, J * omega_bi) );
    
        ok1 = omega_bi_dot(J, N, omega_bi);
     
        omega_bi_s = omega_bi + 0.5*ok1*delta_t;
    
        ok2 = omega_bi_dot(J, N, omega_bi_s);
    
        omega_bi_ss = omega_bi + 0.5*ok2*delta_t;
    
        ok3 = omega_bi_dot(J, N, omega_bi_ss);
    
        omega_bi_sss = omega_bi + ok3*delta_t;
    
        ok4 = omega_bi_dot(J, N, omega_bi_sss);
    
        omega_bi = omega_bi + (1/6) * ( ok1 + 2*ok2 + 2*ok3 + ok4 )*delta_t; 
    end
    
    
    
    % RK4 for calculation of the quaternion
    w_t = omega_bi;
   
    w_th = omega_bi; 
    w_t1 = omega_bi;
    
    k_1 = 0.5*omega(w_t) * quater;
    
    quater_s = quater + k_1*0.5*delta_t;
    
    % the angular velocity is assumed to be discreet, hence angular velocity after half a time step is still equal to the angular velocity at the current time step  
    k_2 = 0.5*omega(w_th) * quater_s;
    
    quater_ss = quater + k_2*0.5*delta_t;
    
    k_3 = 0.5*omega(w_th) * quater_ss;
    
    quater_sss = quater + k_3*delta_t;
    
    k_4 = 0.5*omega(w_t1) * quater_sss;
    
    quater = quater + (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)*delta_t;
    
    quater = quater/(norm(quater));
    
    %mag (t+1, 1) = sqrt( (quater(1,1)^2) + (quater(2,1)^2) + (quater(3,1)^2) + (quater(4,1)^2) );
    
    q_1(t+2,1) = quater(4,1);
    q_2(t+2,1) = quater(1,1);
    q_3(t+2,1) = quater(2,1);
    q_4(t+2,1) = quater(3,1);
    
    
    
    iq = [quater(4,1), quater(1,1), quater(2,1), quater(3,1)];
    
    A_bi = quat2dcm(iq);
    
    z_eci = transpose(A_bi) * z_b;
    
    z_eci_x(t+1,1) = z_eci(1,1);
    z_eci_y(t+1,1) = z_eci(2,1);
    z_eci_z(t+1,1) = z_eci(3,1);
    
    % The counter is incresed by 1
%     t = t + 1;
end

% Time vector is made
time = [0:0.1:(250000*0.1-0.1)];
% time = [0:0.1:250000*0.1];

% Plotting

figure (1)
hold on
tiledlayout('flow')

% Tile 1
nexttile
plot(time, z_eci_x, 'r')
title('z_eci_x vs time')
xlabel('time [s]')
ylabel('z_eci_x')
%  xlim([0 50])
%  ylim([-1.2 1.2])


% Tile 2
nexttile
plot(time, z_eci_y, 'r')
title('z_eci_y vs time')
xlabel('time [s]')
ylabel('z_eci_y')
%  xlim([0 50])
%  ylim([-1.2 1.2])


% Tile 3
nexttile
plot(time, z_eci_z, 'r')
title('z_eci_z vs time')
xlabel('time [s]')
ylabel('z_eci_z')
%  xlim([0 50])
%  ylim([-1.2 1.2])

% % Tile 4
% nexttile
% plot(time, q_4, 'r')
% title('q_4 vs time')
% xlabel('time [s]')
% ylabel('q_4')
%  xlim([0 50])
%  ylim([-1.2 1.2])



hold off

% figure ()
% hold on


% plot(time, mag, 'r')
% title('magnitude vs time')
% xlabel('time [s]')
% ylabel('magnitude')
% % xlim([0 25000])
% % ylim([-1.2 1.2])
% hold off
% 




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



function k = omega_bi_dot(J, N, omega_bi)

    k = inv(J) * ( N - cross( omega_bi, J * omega_bi) );




end
