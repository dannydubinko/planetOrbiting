%% Planet and Spacecraft Orbiting
clear all;
close all;
clc;

%% System Properties and Initial codition 

%Central Force Constant (analogous to gravitational constant G)
G = (6.674e-11);

%Radius for particle 1
rE = 6.4e6;
%here are some random particles
p1 = [0 0 0];
p2 = [-10*rE+3*rE*rand() 0 0];
p3 = [-10*rE+3*rE*rand() 0 0];
p4 = [-10*rE+3*rE*rand() 0 0];

%random assignment of mass to the particles
m1 = 6e24;  %mass of the earth
m2 = 1.5e4; %mass of spacecrafts (p2, p3, p4)
m3 = 1.5e4;
m4 = 1.5e4;

%initial momentum
mv1 = [0 0 0];
mv2 = [0 3e3*m2 0];
mv3 = [0 3e3*m3 0];
mv4 = [0 3e3*m4 0];

%Time interval and step
tmax = 100000;
dt = 10;

%% Centre of Mass 
%write an expression to compute the center of mass (do not use a loop)

%Center of Mass x direction 
com = (m1*p1 + m2*p2 + m3*p3 + m3*p4) / (m1 + m2 + m3 + m4);
%% Initial Conditions Plot


figure
hold on;
plot3(p1(1,1), p1(1,2), p1(1,3),'.r', 'MarkerSize', 20);
plot3(p2(1,1), p2(1,2), p2(1,3),'.g', 'MarkerSize', 20);
plot3(p3(1,1), p3(1,2), p3(1,3),'.b', 'MarkerSize', 20);
plot3(p4(1,1), p4(1,2), p4(1,3),'.m', 'MarkerSize', 20);
plot3(com(1,1), com(1,2), com(1,3),'*k', 'MarkerSize', 25);
view(45,45);

xlabel("X-Position[m]");
ylabel("Y-Position[m]");
zlabel("Z-Position[m]");
title('Initial Positon of Particles');
lgd = legend('Earth', 'Spacecraft 2', 'Spacecraft 3', 'Spacecraft 4', 'Center of Mass');
title(lgd,'Particles','FontSize',12);

hold off;

%% Momentum Update

%Setup arrays for saving positions and momentum/velocity
m = [m1;m2;m3;m4];
p = [p1;p2;p3;p4];

for i = 2:floor((tmax/dt)) %run the simulation for 1 second
    
    %compute the individual gravitational forces acting on each particle
    %e.g. force of p2, p3, and p4 acting on p1. 
    F1 = NetForce(m,p,G,1);
    F2 = NetForce(m,p,G,2);
    F3 = NetForce(m,p,G,3);
    F4 = NetForce(m,p,G,4);

    %update momentum of each particle
    mv1 = mv1 + F1*dt;
    v1 = (mv1/m1);
    p1 = p1 + v1 * dt;

    mv2 = mv2 + F2*dt;
    v2 = (mv2/m2);
    p2 = p2 + v2 * dt;

    mv3 = mv3 + F3*dt;
    v3 = (mv3/m3);
    p3 = p3 + v3 * dt;

    mv4 = mv4 + F4*dt;
    v4 = (mv4/m4);
    p4 = p4 + v4 * dt;
   
    
    %compute updated position of each particle (set new position equal to variable name of previous position)
    p = [p1;p2;p3;p4];
    %create an array for each point that records time history (e.g. p1_out(i,:) = p1;)
    p1_out(i,:) = p1;
    p2_out(i,:) = p2;
    p3_out(i,:) = p3;
    p4_out(i,:) = p4;

    %momentum arrays
    mv1_out(i,:) = mv1;
    mv2_out(i,:) = mv2;
    mv3_out(i,:) = mv3;
    mv4_out(i,:) = mv4;
   
    %recompute center of mass from updated positions

    com = (m1*p1 + m2*p2 + m3*p3 + m4*p4) / (m1 + m2 + m3 + m4);
   
    %create array for center of mass (e.g. com_out(i,:) = com;)
    com_out(i,:) = com;

    %time array
    time(i,:) = dt * i;

end 

%trajectory plots


figure
hold on;
plot3(p1_out(2:end,1), p1_out(2:end,2), p1_out(2:end,3),'.k','markersize', 20);
plot3(p2_out(2:end,1), p2_out(2:end,2), p2_out(2:end,3),'g','LineWidth', 2);
plot3(p3_out(2:end,1), p3_out(2:end,2), p3_out(2:end,3),'b','LineWidth', 2);
plot3(p4_out(2:end,1), p4_out(2:end,2), p4_out(2:end,3),'m','LineWidth', 2);
plot3(com_out(2:end,1), com_out(2:end,2), com_out(2:end,3),'.r','markersize', 1);
view(45,45);

xlabel("X-Trajectory[m]");
ylabel("Y-Trajectory[m]");
zlabel("Z-Trajectory[m]");
title('Trajectories of Particles');
%lgd = legend('Earth', 'Spacecraft 1', 'Spacecraft 2', 'Spacecraft 3', 'Center of Mass');
%title(lgd,'Trajectories','FontSize',12);
hold off;

%% Animation - Question 1E

figure
view(45,45);
for i = 2:30:length(p1_out)
    hold on;
    plot3(p1_out(i,1), p1_out(i,2), p1_out(i,3),'.r','markersize', 13);
    plot3(p2_out(i,1), p2_out(i,2), p2_out(i,3),'.g','markersize', 13);
    plot3(p3_out(i,1), p3_out(i,2), p3_out(i,3),'.b','markersize', 13);
    plot3(p4_out(i,1), p4_out(i,2), p4_out(i,3),'.m','markersize', 13);
    plot3(com_out(i,1), com_out(i,2), com_out(i,3),'.k','markersize', 16);
    drawnow
    xlabel("X-Trajectory[m]");
    ylabel("Y-Trajectory[m]");
    zlabel("Z-Trajectory[m]");
    title('Trajectories of Particles Animation');
   
end
hold off;

%% Linear Momentum 
%Linear Momentum of Particle 1
figure
hold on;
plot(time(:,1),mv1_out(:,1),'r','LineWidth', 2);
plot(time(:,1),mv1_out(:,2),'g','LineWidth', 2);
plot(time(:,1),mv1_out(:,3),'b','LineWidth', 2);

ylabel("Momentum [kgm/s]");
xlabel("Time [s]");
title('Linear Momentum - Particle 1 vs. Time');
lgd = legend('X - Direction', 'Y - Direction', 'Z - Direction');
title(lgd,'Momentum','FontSize',10);
hold off;
%Linear Momentum of Particle 2
figure
hold on;
plot(time(:,1),mv2_out(:,1),'r','LineWidth', 2);
plot(time(:,1),mv2_out(:,2),'g','LineWidth', 2);
plot(time(:,1),mv2_out(:,3),'b','LineWidth', 2);

ylabel("Momentum [kgm/s]");
xlabel("Time [s]");
title('Linear Momentum - Particle 2 vs. Time');
lgd = legend('X - Direction', 'Y - Direction', 'Z - Direction');
title(lgd,'Momentum','FontSize',10);
hold off;
%Linear Momentum of Particle 3
figure
hold on;
plot(time(:,1),mv3_out(:,1),'r','LineWidth', 2);
plot(time(:,1),mv3_out(:,2),'g','LineWidth', 2);
plot(time(:,1),mv3_out(:,3),'b','LineWidth', 2);

ylabel("Momentum [kgm/s]");
xlabel("Time [s]");
title('Linear Momentum - Particle 3 vs. Time');
lgd = legend('X - Direction', 'Y - Direction', 'Z - Direction');
title(lgd,'Momentum','FontSize',10);
hold off;
%Linear Momentum of Particle 4
figure
hold on;
plot(time(:,1),mv4_out(:,1),'r','LineWidth', 2);
plot(time(:,1),mv4_out(:,2),'g','LineWidth', 2);
plot(time(:,1),mv4_out(:,3),'b','LineWidth', 2);

ylabel("Momentum [kgm/s]");
xlabel("Time [s]");
title('Linear Momentum - Particle 4 vs. Time');
lgd = legend('X - Direction', 'Y - Direction', 'Z - Direction');
title(lgd,'Momentum','FontSize',10);
hold off;




%angular momentum of a the system in x,y,z directions
mv_system = mv2_out + mv3_out + mv4_out;
p_out = p2_out + p3_out + p4_out; %spacecrafts
ang_momentum_earth = cross(p1_out, mv1_out);

ang_momentum_system = cross(p_out, mv_system);

figure
hold on;
plot(ang_momentum_earth(:,1),ang_momentum_system(:,1),'r','LineWidth', 2);
plot(ang_momentum_earth(:,2),ang_momentum_system(:,2),'g','LineWidth', 2);
plot(ang_momentum_earth(:,3),ang_momentum_system(:,3),'b','LineWidth', 2);
view (45,45)
ylim([-1,1]);
view(45,45);
ylabel("Angular Momentum System [kgm/s]");
xlabel("Time [s]");
title('Angular Momentum System vs. Time');
lgd = legend('X - Direction', 'Y - Direction', 'Z - Direction');
title(lgd,'Angular Momentum System','FontSize',10);
hold off;

%% Net Force Function
%Here is a function that will take in all the masses
%and there positions to compute the net force acting on the mass of
%interest
function force = NetForce(m,p,G,ind)
%m = all masses in order m1,m2,m3,m4 (4x1)
%p = all positions p1,p2,p3,p4 (4x3)
%G = gravitational constant
%ind = index of the particle we want the net force for (values 1,2,3,4)

%Pull out the position and mass of the ind particle
m_ind = m(ind, 1);
p_ind = p(ind, :);

%Replace the mass and position of the ind particle with a NaN (don't want the force of it on itself)
m(ind, 1) = NaN;
p(ind, :) = [NaN, NaN, NaN];

F = NaN(3,3);
count = 1;

%Compute the force vectors from each of the particles acting on one another
%Hint: Use a loop to do this that goes through all the particles
for i = 1:length(m)
    if ~isnan(m(i, 1))

        r = p(i,:) - p_ind;
        F_mag = (G * m(i,1) * m_ind) / (norm(r)^2 );
        F_vec = F_mag * (r/norm(r));

        F(count,:) = F_vec;
        count = count + 1;


    end
end
%Sum the force vectors for the total net force
force = sum(F);
end
