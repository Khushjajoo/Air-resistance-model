% Khushkumar Jajoo
% Project 2
% Phase 3
% Modelling a home run, with air resistance
% Exporting data and analyzing it in Excel

clf
clear

v0mph = 112; % exit velocity in mph
phi0deg = 32; % launch angle in degrees

x0 = 0; % starting x coordinate of ball
y0 = 0; % starting y coordinate of ball

g = 10; % gravitaional constant in N/kg

mph2mps =  5280 * 12 * 2.54 / 100 / 3600; % mph to mps conversion
deg2rad =  pi()/180; % degrees to radians conversion

v0 = v0mph*mph2mps; % conversion of exit velocity to mps
phi0 = phi0deg*deg2rad; % conversion of launch angle to degrees

v0x = v0*cos(phi0); % x-component of v
v0y = v0*sin(phi0); % y-component of v

%---------- computing other expressions --------------
thmax = v0y/g; % time to reach the maximum height
tflight = 2*thmax; % time of flight
maxheight = thmax*v0y/2; % maximum height
range = v0x*tflight; % range of projectile

% ------- setting up a time array, compute x(t),y(t) --------
m2ft = 3.3; % converting from m to ft

tmin = 0;
tmax = tflight; 
N = 2000; % number of intervals
t = linspace(tmin,tmax,N+1);

xt = (x0 + v0x*t)*m2ft;
yt = (y0 + v0y*t - 0.5*g*t.^2)*m2ft;

% ------------- numeric solution ---------------
m = 0.145; % mass of baseball in kg
C = input("Enter value for drag coefficient: ");
A = 0.0042; % area of cross section in m^2
p = 1.225;  % air density in kg/m^3
dt = (tmax - tmin)/N;

x = zeros(1, N+1); % initializing x(t)
y = zeros(1, N+1); % initializing y(t)

%------- setting initial velocities and positions for x and y -------
x(1) = x0;
y(1) = y0;
vx = v0x; 
vy = v0y;

k = 0.5*C*A*p;

for n = 1:N 
    v = sqrt(vx^2 + vy^2); % calculating velocity using its two components
    Fx = 0 - k*v*vx; % x-component of the net force, in N
    Fy = -m*g - k*v*vy; % y-component of the net force, in N
    ax = Fx/m;
    ay = Fy/m;
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;
    vy = vy + ay*dt;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
end

x_ft = x*m2ft;
y_ft = y*m2ft;

if C == 0
    checkSum_x = sum(abs(x_ft-xt))
    checkSum_y = sum(abs(y_ft-yt))
end
%---------------plotting the trajectory---------------
plot(xt, yt, x_ft, y_ft,'LineWidth', 2)
grid on
ax = gca; 
ax.FontSize = 16;
ax.GridAlpha = 0.5;
set(gca,'XMinorGrid','on');
set(gca,'YMinorGrid','on');
ax.MinorGridAlpha = 0.5;
ylim([0 130])
title({'ECE 202, Project 2: Phase 3: Trajectory of a baseball ', ... 
        'without drag vs with drag'}, 'FontSize', 22)
xlabel('x (ft)', 'FontSize', 18)
ylabel('y (ft)', 'FontSize', 18)
str = sprintf("with drag, C = %g",C);
legend({'without drag', str}, ...
        'FontSize', 18)
%-------------- exporting data to excel -----------------
export = [t; x_ft; y_ft].';
labels = ["Time t (s)", "x (ft)","y (ft)"];
export = [labels; export];
writematrix(export, 'homerun.csv')