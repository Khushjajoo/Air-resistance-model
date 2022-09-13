% Khushkumar Jajoo
% Project 2
% Phase 1
% Modelling a home run, with air resistance
% Comparing the analytic solutions to the numeric, without drag

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
dt = (tmax - tmin)/N;

x = zeros(1, N+1); % initializing x(t)
y = zeros(1, N+1); % initializing y(t)

%------- setting initial velocities and positions for x and y -------
x(1) = x0;
y(1) = y0;
vx = v0x; 
vy = v0y;

for n = 1:N 
    Fx = 0;    % x-component of the net force, in N
    Fy = -m*g; % y-component of the net force, in N
    ax = Fx/m;
    ay = Fy/m;
    y(n+1) = y(n) + vy*dt + (1/2)*ay*dt^2;
    vy = vy + ay*dt;
    x(n+1) = x(n) + vx*dt + (1/2)*ax*dt^2;
    vx = vx + ax*dt;
end

xft = x*m2ft;
yft = y*m2ft;


checkSum_x = sum(abs(xft-xt))
checkSum_y = sum(abs(yft-yt))

%---------------plotting the trajectory---------------
plot(xt, yt, xft, yft,'LineWidth', 2)
grid on
ax = gca; 
ax.FontSize = 16;
ax.GridAlpha = 0.5;
title({'ECE 202, Project 2: Phase 1: Trajectory of a baseball with ', ... 
       'no drag. Analytic vs Numeric solution'}, 'FontSize', 22)
xlabel('x (ft)', 'FontSize', 18)
ylabel('y (ft)', 'FontSize', 18)
legend({'analytic (behind numeric)', 'numeric'}, ...
        'FontSize', 18)
