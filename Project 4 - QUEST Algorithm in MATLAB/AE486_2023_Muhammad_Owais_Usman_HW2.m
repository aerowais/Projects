clear all

load ("hw2_data");

%% Q1

% Defining constants and initialising
t = 1;

delta_t = 1;


for t = 1:25000
    
    r_i = pos_eci (:,t);
    v_i = vel_eci (:,t);
    
    rot = c_oi(r_i, v_i);
    
    
    m0_1(:,t) =  rot *( mag_eci(:,t));
    
    
    
    
    t = t+1;

end

m11 = transpose(m0_1(1,:));
m12 = transpose(m0_1(2,:));
m13 = transpose(m0_1(3,:));







% Using dipole model

% initializing and defining constants
t = 1;

m_e = 7.71e15;
w_e = 7.29e-5;
epsilon = deg2rad(9.3);
inc = deg2rad(98.1245);
r_0 = 7056198.5;
mu = 3.986004418e14;
w_0 = sqrt((mu)/(r_0^3));

for t = 1:25000
    
    b_1(1,t) = ((m_e)/(r_0^3)) * ( cos(w_0*t)*( cos(epsilon)*sin(inc) - sin(epsilon)*cos(inc)*cos(w_e*t) ) - sin(w_0*t)*sin(epsilon)*sin(w_e*t) );
    
    b_2(1,t) = -((m_e)/(r_0^3)) * ( cos(epsilon)*cos(inc) + sin(epsilon)*sin(inc)*cos(w_e*t) );
    
    b_3(1,t) = 2 * ((m_e)/(r_0^3)) * ( sin(w_0*t)*( cos(epsilon)*sin(inc) - sin(epsilon)*cos(inc)*cos(w_e*t) ) - 2*sin(w_0*t)*sin(epsilon)*sin(w_e*t) );
    
    m0_2 = (10^(9))*vertcat( b_1, b_2, b_3);
    
    
    t = t+1;

    
end

m21 = transpose(m0_2(1,:));
m22 = transpose(m0_2(2,:));
m23 = transpose(m0_2(3,:));
    

% Plotting to compare

time = [1:25000];

figure (1)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, m11, 'b')
title('m0_1_x vs time')
xlabel('time [s]')
ylabel('m0_1_x [nT]')


% Tile 2
nexttile
plot(time, m12, 'b')
title('m0_1_y vs time')
xlabel('time [s]')
ylabel('m0_1_y [nT]')

% Tile 3
nexttile
plot(time, m13, 'b')
title('m0_1_z vs time')
xlabel('time [s]')
ylabel('m0_1_z [nT]')



% Tile 1
nexttile
plot(time, m21, 'r')
title('m0_2_x vs time')
xlabel('time [s]')
ylabel('m0_2_x [nT]')


% Tile 2
nexttile
plot(time, m22, 'r')
title('m0_2_y vs time')
xlabel('time [s]')
ylabel('m0_2_y [nT]')

% Tile 3
nexttile
plot(time, m23, 'r')
title('m0_2_z vs time')
xlabel('time [s]')
ylabel('m0_2_z [nT]')

hold off

%-----------------------------------------------------------------------------------------------------------------------------

%% Q2


t = 1;

delta_t = 1;


for t = 1:25000
    
    r_i = pos_eci (:,t);
    v_i = vel_eci (:,t);
    
    rot = c_oi(r_i, v_i);
    
    
    s0_1(:,t) =  rot * sun_eci(:,t);
   

end

s11 = transpose(s0_1(1,:));
s12 = transpose(s0_1(2,:));
s13 = transpose(s0_1(3,:));

time = [1:25000];

% Plotting
figure (2)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, s11, 'r')
title('s0_1_x vs time')
xlabel('time [s]')
ylabel('s0_1_x')


% Tile 2
nexttile
plot(time, s12, 'r')
title('s0_1_y vs time')
xlabel('time [s]')
ylabel('s0_1_y')

% Tile 3
nexttile
plot(time, s13, 'r')
title('s0_1_z vs time')
xlabel('time [s]')
ylabel('s0_1_z')

hold off

%%  Q3

% All the angles except roll will be kept at zero as the distrubanc is only in roll

psi   = zeros(25000,1); %yaw
theta = zeros(25000,1); %pitch
phi   = zeros(25000,1); %roll
q1 = zeros(25000,1);
q2 = zeros(25000,1);
q3 = zeros(25000,1);
q4 = zeros(25000,1);

q4(1,1) = 1;

t = 1;



max_amp = deg2rad(10); 

omega = 2*pi/2000;

att  = zeros(3,3);

for t= 1:25000

    phi(t,1) = max_amp*sin(omega*t);
    
    att(:,:,t) = angle2dcm( psi(t,1), theta(t,1), phi(t,1));
    
    quaternion = dcm2quat( att(:,:,t) );
    quaternion = transpose(quaternion);
    
    q4(t,1) = quaternion(1,1);
    q1(t,1) = quaternion(2,1);
    q2(t,1) = quaternion(3,1);
    q3(t,1) = quaternion(4,1);
    

end  

% Time vector is made

time = [1:25000];

% Plotting

figure (4)
hold on
tiledlayout(2,2);

% Tile 1
nexttile
plot(time, q1, 'b')
title('q_1 vs time')
xlabel('time [s]')
ylabel('q_1')


% Tile 2
nexttile
plot(time, q2, 'b')
title('q_2 vs time')
xlabel('time [s]')
ylabel('q_2')

% Tile 3
nexttile
plot(time, q3, 'b')
title('q_3 vs time')
xlabel('time [s]')
ylabel('q_3')


% Tile 4
nexttile
plot(time, q4, 'b')
title('q_4 vs time')
xlabel('time [s]')
ylabel('q_4')


hold off

% Plotting
figure (5)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time,psi, 'b')
title('psi vs time')
xlabel('time [s]')
ylabel('psi [rad]')


% Tile 2
nexttile
plot(time,theta, 'g')
title('theta vs time')
xlabel('time [s]')
ylabel('theta [rad]')

% Tile 3
nexttile
plot(time,phi, 'r')
title('phi vs time')
xlabel('time [s]')
ylabel('phi [rad]')

hold off

%% Q4


t = 1;

for t = 1:25000

    mb(:,t) = att(:,:,t) * m0_1(:,t) + randn*[300;300;300];
    
    sb(:,t) = att(:,:,t) * s0_1(:,t) + randn*[0.002;0.002;0.002];    
    
    
end
    
mbx = transpose(mb(1,:));
mby = transpose(mb(2,:));
mbz = transpose(mb(3,:));

sbx = transpose(sb(1,:));
sby = transpose(sb(2,:));
sbz = transpose(sb(3,:));

    

% Plotting to compare

time = [1:25000];

figure (6)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, mbx, 'b')
title('m_b_x vs time')
xlabel('time [s]')
ylabel('m_b_x [nT]')


% Tile 2
nexttile
plot(time, mby, 'b')
title('m_b_y vs time')
xlabel('time [s]')
ylabel('m_b_y [nT]')

% Tile 3
nexttile
plot(time, mbz, 'b')
title('m_b_z vs time')
xlabel('time [s]')
ylabel('m_b_z [nT]')



% Tile 1
nexttile
plot(time, sbx, 'r')
title('s_b_x vs time')
xlabel('time [s]')
ylabel('s_b_x')


% Tile 2
nexttile
plot(time, sby, 'r')
title('s_b_y vs time')
xlabel('time [s]')
ylabel('s_b_y')

% Tile 3
nexttile
plot(time, sbz, 'r')
title('s_b_z vs time')
xlabel('time [s]')
ylabel('s_b_z')

hold off

%% Q5

t = 1
for t= 1:25000

    quat2(:,1,t) = quest( sb(:,t) , s0_1(:,t) , mb(:,t), m0_1(:,t) );
    
    
    
    
    q12(t,1) = quat2(1,1,t);
    q22(t,1) = quat2(2,1,t);
    q32(t,1) = quat2(3,1,t);
    q42(t,1) = quat2(4,1,t);

    %quat2(1,1) = quat2(4,1);
    %quat2(2,1) = quat2(1,1);
    %quat2(3,1) = quat2(2,1);
    %quat2(4,1) = quat2(3,1);
    
    
    %att2 = quat2dcm(quat2);
    euler = quat2eul( transpose(quat2(:,1,t)), "ZYX" );
    
    psi2(t,1)  = euler(1,1);
    theta2(t,1) = euler(1,2);
    phi2(t,1)   = euler(1,3); 
    
    
    
    
end  

% Time vector is made

time = [1:25000];

% Plotting

figure (7)
hold on
tiledlayout(2,2);

% Tile 1
nexttile
plot(time, q12, 'b')
title('q_1_2 vs time')
xlabel('time [s]')
ylabel('q_1_2')


% Tile 2
nexttile
plot(time, q22, 'b')
title('q_2_2 vs time')
xlabel('time [s]')
ylabel('q_2_2')

% Tile 3
nexttile
plot(time, q32, 'b')
title('q_3_2 vs time')
xlabel('time [s]')
ylabel('q_3_2')


% Tile 4
nexttile
plot(time, q42, 'b')
title('q_4_2 vs time')
xlabel('time [s]')
ylabel('q_4_2')


hold off

% Plotting
figure (8)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time,psi2, 'b')
title('psi_2 vs time')
xlabel('time [s]')
ylabel('psi_2 [rad]')


% Tile 2
nexttile
plot(time,theta2, 'g')
title('theta_2 vs time')
xlabel('time [s]')
ylabel('theta_2 [rad]')

% Tile 3
nexttile
plot(time,phi2, 'r')
title('phi_2 vs time')
xlabel('time [s]')
ylabel('phi_2 [rad]')

hold off












%% ---------------------------------------------FUNCTIONS------------------------------------------------------------------------

% function to get orbital from ECI

function c_oi_matrix = c_oi(r,v)

    o3 = -r*(1/sqrt( (r(1,1)^2) + (r(2,1)^2) + (r(3,1)^2) ));
    cr = cross(r, v);
    o2 = -cr*(1/sqrt( (cr(1,1)^2) + (cr(2,1)^2) + (cr(3,1)^2) ));
    o1 = cross(o2,o3);
    
    c_io = horzcat( o1 ,o2 , o3 );
    
   c_oi_matrix = transpose(c_io);
end

% function for QUEST algorithm

function quater = quest( sb , s0 , mb , m0 )
    
    mb = mb/norm(mb);
    m0 = m0/norm(m0);

    B = (5/6) * sb * transpose(s0) + (1/6) * mb * transpose(m0);

    S = B + transpose(B);
    s = trace(B);
    z = [ B(2,3)-B(3,2) ; B(3,1)-B(1,3); B(1,2)-B(2,1)];
    
    a = (s^2) - trace( adjoint(S) );
    b = (s^2) + transpose(z) * z;
    c = 8* det(B);
    d = transpose(z) * (S^2) * z;
    
    lm = 1;
    j = 1;
    
    for j = 1:15
        
        lm = lm - ( ( (lm^4) - (a + b)*(lm^2) - c*lm + (a*b + c*s - d) )/( 4*(lm^3) - 2*(a+b)*lm - c )   );                          
        
    end

    alpha = (lm^2) - a;
    gamma = alpha*(lm + s) - det(S);
    x     = (  alpha*eye(3) + (lm - s)*S + (S^2) ) * z;
    
    
    quater = 1/sqrt( (gamma^2) + (norm(x)^2) ) * [ x(1,1); x(2,1); x(3,1); gamma];
    
end
