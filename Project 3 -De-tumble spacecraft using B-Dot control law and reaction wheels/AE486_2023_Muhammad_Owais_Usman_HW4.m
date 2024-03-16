clear all

load ("hw2_data");
%% ============================================== Q1 ============================================== %%
% Initializations and constants
w_bi_1(:,1) = [0.2; 0.1;1.5];

J = [6.9 0 0;
     0 7.5 0;
     0 0 8.4];
J_diag = [J(1,1),J(2,2),J(3,3)];

J_min = min(J_diag);


q1 = (0.5*sqrt(2)) * [1;0;0;1];

qua1= [q1(4,1,1), q1(1,1,1), q1(2,1,1), q1(3,1,1)];
               
               
             
           
A_bi_1(:,:,1) = quat2dcm( qua1 );         



epsilon = deg2rad(98.1245);

r_0 = 7056198.5;
G   = 6.67430e-11;
M_e = 5.972e+24;

T_orb = (2*pi*(r_0^1.5))/sqrt(G*M_e);

k = ( (4*pi)/T_orb ) *( 1 + sin(epsilon) ) * J_min;



j = 1;

delta = 1;

for j = 2:25000
    
        col_1 = A_bi_1(:,1,j-1);
        col_2 = A_bi_1(:,2,j-1);
        col_3 = A_bi_1(:,3,j-1);

        vec_1 = col_1;
        vec_2 = col_2 - ( (dot(vec_1,col_2))/(norm(vec_1)^2) )*vec_1;
        vec_3 = col_3 - ( (dot(vec_1,col_3))/(norm(vec_1)^2) )*vec_1;
        vec_3 = vec_3 - ( (dot(vec_2,col_3))/(norm(vec_2)^2) )*vec_2;

        e_1 = vec_1/norm(vec_1);
        e_2 = vec_2/norm(vec_2);
        e_3 = vec_3/norm(vec_3);

        A_bi_1(:,1,j-1) = e_1;
        A_bi_1(:,2,j-1) = e_2;
        A_bi_1(:,3,j-1) = e_3;
    
    
    
    B = A_bi_1(:,:,j-1) * mag_eci(:,j-1);
    B =B/norm(B);
%======= Calculating the control torque =====%
    m(:,j-1) = ( (-k) ) * cross(B,w_bi_1(:,j-1));
    
    L_1(:,j-1) = cross(m(:,j-1),B);
    
%======= Finding new values for the w_bi, and A_bi =====%

    % Calculating w_bi using classical RK-4
        k1w = w_dot( J, L_1(:,j-1), w_bi_1(:,j-1));
        
        k2w = w_dot( J, L_1(:,j-1), w_bi_1(:,j-1)+0.5*delta*k1w );
        
        k3w = w_dot( J, L_1(:,j-1), w_bi_1(:,j-1)+0.5*delta*k2w );
        
        k4w = w_dot( J, L_1(:,j-1), w_bi_1(:,j-1)+delta*k3w );
        
        w_bi_1(:,j) = w_bi_1(:,j-1) + (1/6) * ( k1w + 2*k2w + 2*k3w + k4w ) * delta;
        
    % Calculating A_bi using classical RK-4
    
        k1a = A_bi_dot( w_bi_1(:,j-1), A_bi_1(:,:,j-1) );
        
        k2a = A_bi_dot( w_bi_1(:,j-1), A_bi_1(:,:,j-1)+0.5*delta*k1a );
       
        k3a = A_bi_dot( w_bi_1(:,j-1), A_bi_1(:,:,j-1)+0.5*delta*k2a );
        
        k4a = A_bi_dot( w_bi_1(:,j-1), A_bi_1(:,:,j-1)+delta*k3a );
        
        A_bi_1(:,:,j) = A_bi_1(:,:,j-1) + (1/6) * ( k1a + 2*k2a + 2*k3a + k4a ) * delta;
    
end

time = [1:25000];

%== Ploting ==%

figure (1)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, w_bi_1(1,:), 'b')
title('w_bi_x vs time')
xlabel('time [s]')
ylabel('w_bi_x [rad/s]')
% xlim([0 10000])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, w_bi_1(2,:), 'b')
title('w_bi_y vs time')
xlabel('time [s]')
ylabel('w_bi_y [rad/s]')
% xlim([0 10000])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, w_bi_1(3,:), 'b')
title('w_bi_z vs time')
xlabel('time [s]')
ylabel('w_bi_z [rad/s]')
% xlim([0 10000])
% ylim([-1.2 1.2])

% % Tile 4
% nexttile
% plot(time, q_4, 'b')
% title('q_4 vs time')
% xlabel('time [s]')
% ylabel('q_4')
% xlim([0 25000])
% ylim([-0.2 1.2])

hold off

%------------------------%
time_m = [1:24999];

m = m*(10^9);
figure (2)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time_m, m(1,:), 'b')
title('m_x vs time')
xlabel('time [s]')
ylabel('m_x [A m^2]')
xlim([0 1000])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time_m, m(2,:), 'b')
title('m_y vs time')
xlabel('time [s]')
ylabel('m_y [A m^2]')
xlim([0 1000])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time_m, m(3,:), 'b')
title('m_z vs time')
xlabel('time [s]')
ylabel('m_z [A m^2]')
xlim([0 1000])
% ylim([-1.2 1.2])

hold off

%-----------------%

figure (3)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time_m, L_1(1,:), 'b')
title('L_x vs time')
xlabel('time [s]')
ylabel('L_x [Nm]')
xlim([0 1000])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time_m, L_1(2,:), 'b')
title('L_y vs time')
xlabel('time [s]')
ylabel('L_y [Nm]')
xlim([0 1000])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time_m, L_1(3,:), 'b')
title('L_z vs time')
xlabel('time [s]')
ylabel('L_z [Nm]')
xlim([0 1000])
% ylim([-1.2 1.2])

hold off

%% ============================================== Q2 ============================================== %%

% Initializations
a = 1/sqrt(2);
b = 1/sqrt(2);
c = 1/sqrt(2);
d = 1/sqrt(2);

W_inverse = [1/a, b/((b^2) + (c^2)), 0;
      -1/a, b/((b^2) + (c^2)),0;
      0, c/((b^2) + (c^2)), 1/d;
      0, c/((b^2) + (c^2)), -1/d];
             
q(:,1) = [0.6853;
          0.6953;
          0.1531;
          0.1531];
      
q(:,1) = q(:,1)/norm(q(:,1));

z_b         = [0; 0; 1];
z_i         = sun_eci(:,1);
phi         = 1.524; % [rad]
tht         = 1.591; % [rad]
psi         = 0;     % [rad]
A_bi_c = angle2dcm(phi,tht,psi,"ZXZ");
q_c1    = angle2quat(phi,tht,psi,"ZXZ");

q_c = [q_c1(1,2);
       q_c1(1,3);
       q_c1(1,4);
       q_c1(1,1)]; 

q_c = q_c/norm(q_c);

w_bi(:,1) = [0.53;
             0.53;
             0.053];
         
w_bi(:,1) = (pi/180)*w_bi(:,1);
k_p = 0.1;
k_d = 0.15;
         
j=0;
h(:,1) = [0;0;0];
hw(:,1) = W_inverse * h(:,1);
for j=2:25000

    [ delta_q_1_3 , delta_q_4]  = cross_p( q(:,j-1),q_c );
    
    delta_q0(:,j-1) = [ delta_q_1_3 ; delta_q_4];
    
    delta_q_1_3 = delta_q_1_3/norm(delta_q0(:,j-1));
    delta_q_4   = delta_q_4/norm(delta_q0(:,j-1));
    
    L(:,j-1) = -k_p * sign(delta_q_4) * delta_q_1_3 - k_d * w_bi(:,j-1);
    
    h_dot2(:,j-1) = h_dot( L(:,j-1), w_bi(:,j-1), h(:,j-1) );
    
    %======= Finding new values for the w_bi, and q =====%

    % Calculating w_bi using classical RK-4
  
        k1w2 = w_dot2( J,h(:,j-1),h_dot2(:,j-1),w_bi(:,j-1));
        
        k2w2 = w_dot2( J,h(:,j-1),h_dot2(:,j-1),w_bi(:,j-1)+0.5*delta*k1w2 );
        
        k3w2 = w_dot2( J,h(:,j-1),h_dot2(:,j-1), w_bi(:,j-1)+0.5*delta*k2w2 );
        
        k4w2 = w_dot2( J,h(:,j-1),h_dot2(:,j-1), w_bi(:,j-1)+delta*k3w2 );
        
        w_bi(:,j) = w_bi(:,j-1) + (1/6) * ( k1w2 + 2*k2w2 + 2*k3w2 + k4w2 ) * delta;
    
    % Calculating q using classical RK-4
        k1q = q_dot( w_bi(:,j-1), q(:,j-1) );
        
        k2q = q_dot( w_bi(:,j-1), q(:,j-1)+0.5*delta*k1q );
        
        k3q = q_dot( w_bi(:,j-1), q(:,j-1)+0.5*delta*k2q );
        
        k4q = q_dot( w_bi(:,j-1), q(:,j-1)+delta*k3q );
        
        q(:,j) = q(:,j-1) + (1/6) * ( k1q + 2*k2q + 2*k3q + k4q ) * delta;    

        q(:,j) = q(:,j)/norm(q(:,j));
        
  % Calculating h using classical RK-4
        k1h = h_dot( L(:,j-1), w_bi(:,j-1), h(:,j-1) );
        
        k2h = h_dot( L(:,j-1), w_bi(:,j-1), h(:,j-1)+0.5*delta*k1h );
        
        k3h = h_dot( L(:,j-1), w_bi(:,j-1), h(:,j-1)+0.5*delta*k2h );
        
        k4h = h_dot( L(:,j-1), w_bi(:,j-1), h(:,j-1)+delta*k3h );
        
        h(:,j) = h(:,j-1) + (1/6) * ( k1h + 2*k2h + 2*k3h + k4h ) * delta; 
        
   % Finding the wheel torques using the pseudoinverse law
   
   hw(:,j) = W_inverse * h(:,j);
        
end

%== Ploting ==%

figure (4)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, w_bi(1,:), 'b')
title('w_bi_x vs time')
xlabel('time [s]')
ylabel('w_bi_x [rad/s]')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, w_bi(2,:), 'b')
title('w_bi_y vs time')
xlabel('time [s]')
ylabel('w_bi_y [rad/s]')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, w_bi(3,:), 'b')
title('w_bi_z vs time')
xlabel('time [s]')
ylabel('w_bi_z [rad/s]')
xlim([0 1500])
% ylim([-1.2 1.2])


hold off

figure (5)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, q(1,:), 'b')
title('q_1 vs time')
xlabel('time [s]')
ylabel('q_1')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, q(2,:), 'b')
title('q_2 vs time')
xlabel('time [s]')
ylabel('q_2')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, q(3,:), 'b')
title('q_3 vs time')
xlabel('time [s]')
ylabel('q_3')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 4
nexttile
plot(time, q(4,:), 'b')
title('q_4 vs time')
xlabel('time [s]')
ylabel('q_4')
xlim([0 1500])
% ylim([-1.2 1.2])

hold off


figure (6)
hold on
tiledlayout('flow');

% Tile 1
nexttile
plot(time, hw(1,:), 'b')
title('hw_1 vs time')
xlabel('time [s]')
ylabel('hw_1 [Nms]')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 2
nexttile
plot(time, hw(2,:), 'b')
title('hw_2 vs time')
xlabel('time [s]')
ylabel('hw_2 [Nms]')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 3
nexttile
plot(time, hw(3,:), 'b')
title('hw_3 vs time')
xlabel('time [s]')
ylabel('hw_3 [Nms]')
xlim([0 1500])
% ylim([-1.2 1.2])

% Tile 4
nexttile
plot(time, hw(4,:), 'b')
title('hw_4 vs time')
xlabel('time [s]')
ylabel('hw_4 [Nms]')
xlim([0 1500])
% ylim([-1.2 1.2])

hold off











%============================ Functions ============================%

%=== Function for w_bi using magnetic torquer====
function w_bi_dot = w_dot( J , L , w_bi )
    
    w_bi_dot = inv(J) * ( L - cross(w_bi, J*w_bi) );
    
    
end




%==== Function for A_bi_dot======

% Function to calculate the omega matrix for 
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

% Function to calculate A_bi_dot 
function Abi_dot = A_bi_dot( w_bi, A_bi )
    
    Abi_dot = (-cpm(w_bi) ) * A_bi;
    
end

%========= Function for calculating q_dot====

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

function q_dot_1 = q_dot(w_bi,q_a)

    q_dot_1 = 0.5 * omega(w_bi) * q_a;
    
end

%===== Function for quaternion cross product =====

function [p_cross13,p_cross4] = cross_p(q_1,q_2)

    hat = [q_2(4,1), -q_2(3,1), q_2(2,1);
           q_2(3,1), q_2(4,1), -q_2(1,1);
           -q_2(2,1), q_2(1,1), q_2(4,1);
           -q_2(1,1), -q_2(2,1), -q_2(3,1)];
    
    p_cross13 = transpose(hat)*q_1;
    
    p_cross4 = transpose(q_1)*q_2;

end

%=== Function for w_bi using reaction wheel====
function w_bi_dot_2 = w_dot2( J, h1, h1_dot, w_bi)
    
    w_bi_dot_2 = inv(J) * ( - cross(w_bi, J*w_bi + h1) - h1_dot );
    
end

%==== Function for h_dot =====%

function h_dot_1 = h_dot(L,w_bi,h)

    h_dot_1 = -cross(w_bi,h) - L;
    
end
