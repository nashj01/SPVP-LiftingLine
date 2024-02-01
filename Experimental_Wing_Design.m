clear all
close all

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experimental Wing Design                                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Created by Jordan Nash                                                                  %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% This program is used to evaluate the infinite wing properties of the
% aircrafts wing by considering a user specified profile which is
% influenced by the thickness to chord ratio. It also evaluates the effect
% of flap deflection and determines the aircrafts finite wing properties.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
t_max = 0.12;
x_mc = 0.4;
y_mc = 0.02;
x_f = 0.8;
y_f = 0.005;
theta_deg = 15;
alpha_deg = 12; % Maximum of 14 deg AoA
%
Number_of_stations = 10;
lambda = 0.45;
c_r = 2;
AR = 9;
a_2d = 2*pi;
a_3d = (a_2d/(1+(a_2d/(pi*AR))));
b_span = 14.89;
CD_0 = 0.005;
epsilon_deg = -2;
gamma_add_deg = 2.3; % Additional Sweep
res = 0.01;
v_cruise = 160;
M_mto = 7292.09;
rho_cruise = 0.90926;
a_cruise = 328.6;
M_cruise = v_cruise/a_cruise;
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Initialisation
theta = theta_deg*pi/180;
x_c_0 = 0:res:1;
alpha_geo = alpha_deg*pi/180;
S_wing = b_span^2/AR;
gamma_add = gamma_add_deg*pi/180;
epsilon = epsilon_deg*pi/180;
npts = 1;

% Planform
c_t = 2*(S_wing/b_span)-c_r;
subplot(2,3,2)
plot([0 c_r (b_span/2)*sin(gamma_add)+(c_r-c_t)/2+c_t (b_span/2)*sin(gamma_add)+(c_r-c_t)/2 0], [0 0 (b_span/2)-(b_span-b_span*cos(gamma_add)) (b_span/2)-(b_span-b_span*cos(gamma_add)) 0], 'k')
hold on
plot([0 c_r (b_span/2)*sin(gamma_add)+(c_r-c_t)/2+c_t (b_span/2)*sin(gamma_add)+(c_r-c_t)/2 0], [0 0 (-b_span/2)+(b_span-b_span*cos(gamma_add)) (-b_span/2)+(b_span-b_span*cos(gamma_add)) 0], 'k')
plot([((b_span/2)*sin(gamma_add)+(c_r-c_t)/2)+0.25*c_t 0.25*c_r],[(b_span/2)-(b_span-b_span*cos(gamma_add)) 0],'-.r')
plot([((b_span/2)*sin(gamma_add)+(c_r-c_t)/2)+0.25*c_t 0.25*c_r],[-((b_span/2)-(b_span-b_span*cos(gamma_add))) 0],'-.r')
% (Effective Sweep)
gamma = asin(((b_span/2)*sin(gamma_add)+((c_r-c_t)/2)-0.25*c_r+0.25*c_t)/(b_span/2));
gamma_deg = (180/pi)*gamma;
axis equal
axis([-c_r 3*c_r (-b_span/2)-1 (b_span/2)+1])
c_mac = 2*((((c_t/2)-(c_r/2))/(b_span/2))*(b_span/4) + (c_r/2));
title('Planform')
xlabel('c [m]')
ylabel('b [m]')
% (Leading Edge Sweep)
gamma_le = asin(((b_span/2)*sin(gamma_add)+((c_r-c_t)/2))/(b_span/2));
gamma_le_deg = (180/pi)*gamma_le;

% NACA Profile
for x = 0:res:1
    y_t(1,npts) = (t_max/2)*(2.980*(x^0.5)-1.320*(x)-3.286*(x^2)+2.441*(x^3)-0.815*(x^4));
    if x<x_mc
    y_c_0(1,npts) = y_mc*(2*((x)/(x_mc))-(((x)/(x_mc))^2));
    beta(1,npts) = atan(2*(y_mc/x_mc)*(1-(x/x_mc)));
    else
    y_c_0(1,npts) = y_mc*((2*((1-x)/(1-x_mc))-((1-x)/(1-x_mc))^2));
    beta(1,npts) = atan(-2*(y_mc/(1-x_mc))*(1-((1-x)/(1-x_mc))));
    end
    npts = npts+1;
end
npts = 1;
for x = 0:res:1
n = x;
r(1,npts) = sqrt((((y_c_0(1,npts)-y_f))^2)+((x-x_f)^2));
psi(1,npts) = atan((y_c_0(1,npts)-y_f)/(x-x_f));
    if x<x_f
    x_c(1,npts) = x_c_0(1,npts);
    y_c(1,npts) = y_c_0(1,npts);
    else
    x_c(1,npts) = x_f + r(1,npts)*cos(theta-psi(1,npts));
    y_c(1,npts) = y_f - r(1,npts)*sin(theta-psi(1,npts));
    end
x_u(1,npts) = x_c(1,npts) - y_t(1,npts)*sin(beta(1,npts));
y_u(1,npts) = y_c(1,npts) + y_t(1,npts)*cos(beta(1,npts));
x_l(1,npts) = x_c(1,npts) + y_t(1,npts)*sin(beta(1,npts));
y_l(1,npts) = y_c(1,npts) - y_t(1,npts)*cos(beta(1,npts));
npts = npts+1;
end

% Thin Airfoil Theory
for npts = 1:2
    if npts == 1
        theta_def = 0;
    else
        theta_def = theta;
    end
    syms x phi
    d_y_c1 = simplify(y_mc*((2*((1-x)/(1-x_mc))-((1-x)/(1-x_mc))^2)));
    d_y_c1 = diff(d_y_c1);
    d_y_c1 = subs(d_y_c1,x,(1/2)*(1-cos(phi)));
    d_y_c2 = tan((atan(-2*(y_mc/(1-x_mc))*(1-((1-x)/(1-x_mc)))))-theta_def);
    d_y_c2 = subs(d_y_c2,x,(1/2)*(1-cos(phi)));
    % (A0)
    syms x phi
    a0 = simplify(d_y_c1);
    a0 = double(int(a0, 0, acos(1-2*x_f)));
    a1 = simplify(d_y_c2);
    a1 = double(int(a1, acos(1-2*x_f), pi));
    A0 = alpha_geo - (1/pi)*(a0+a1);
    % (A1)
    syms x phi
    a0 = simplify(d_y_c1*cos(phi));
    a0 = double(int(a0, 0, acos(1-2*x_f)));
    a1 = simplify(d_y_c2*cos(phi));
    a1 = double(int(a1, acos(1-2*x_f), pi));
    A1 = (2/pi)*(a0+a1);
    % (A2)
    syms x phi
    a0 = simplify(d_y_c1*cos(2*phi));
    a0 = double(int(a0, 0, acos(1-2*x_f)));
    a1 = simplify(d_y_c2*cos(2*phi));
    a1 = double(int(a1, acos(1-2*x_f), pi));
    A2 = (2/pi)*(a0+a1);
    % (Cl)
    Cl(1,npts) = pi*(2*A0 + A1);
    % (CL)
    CL(1,npts) = Cl(1,npts)/(1+(Cl(1,npts)/(pi*AR)));
    % (Cm_le)
    Cm_le(1,npts) = -((Cl(1,npts)/4)+(pi/4)*(A1-A2));
    % (Cm_c4)
    Cm_c4(1,npts) = -(pi/4)*(A1-A2);
    % (x_cp)
    x_cp(1,npts) = (1/4)*(1+(pi/Cl(1,npts))*(A1-A2));
    % (Alpha_L0)
    syms x phi
    a0 = simplify(d_y_c1*(cos(phi)-1));
    a0 = double(int(a0, 0, acos(1-2*x_f)));
    a1 = simplify(d_y_c2*(cos(phi)-1));
    a1 = double(int(a1, acos(1-2*x_f), acos(1-2*1)));
    A_L0 = (2/pi)*(a0+a1);
    alpha_L0(1,npts) = (-1/pi)*A_L0;
    alpha_L0_deg(1,npts) = alpha_L0(1,npts)*180/pi;
end

% Airfoil Plot
subplot(2,3,1)
plot(x_c,y_c)
hold on
plot(x_u,y_u)
plot(x_l,y_l)
axis equal 
title('Airfoil Profile - Flap Model')
xlabel('c [%]')
ylabel('$\frac{t}{c}$ [\%]','Interpreter','latex')

% Angle of Incidence
Cl_cruise = 2*M_mto*9.81/(rho_cruise*(v_cruise^2)*S_wing);
alpha_inc_root = Cl_cruise/((2*pi)/(1+(2*pi)/(pi*AR))) + alpha_L0(1,1) - 0.4*epsilon_deg*(pi/180);
alpha_inc_root_deg = alpha_inc_root*(180/pi);
alpha_inc_tip = alpha_inc_root + epsilon_deg*pi/180;
alpha_inc_tip_deg = alpha_inc_tip*180/pi;


% Evaluating Resolution
arr_size = 0:res:1;

% Determining Angle 
for station = 1:Number_of_stations
    upsilon(station, 1) = ((pi/2)/Number_of_stations)*station;
end

% Calculating the trig functions
odd = 1;
for column = 1:Number_of_stations
    for row = 1:Number_of_stations
    sin_functions(row, column) = sin(odd*upsilon(row, 1));
    end
    odd = odd+2;
end

for row = 1:Number_of_stations
    cos_functions(row, 1) = cos(upsilon(row, 1));
end
for row = 1:Number_of_stations
mu(row, 1) = ((a_2d)/(2*AR*(1+lambda))*(1+(lambda-1)*cos_functions(row, 1)));
end

% Evaluating RHS 
odd = 1;
for column = 1:Number_of_stations
    for row = 1:Number_of_stations
        A(row, column) = sin_functions(row, column)*(odd*mu(row, 1)+sin_functions(row, 1));
    end
    odd = odd+2;
end

% Evaluating LHS 
for ver = 1:4
    for row = 1:Number_of_stations
        epsilon(row,1) = (epsilon_deg*pi/180)*(row/Number_of_stations);
        if ver == 1
        B(row, 1) = mu(row, 1)*(alpha_geo-alpha_L0(1,1)-epsilon(row,1)-(((alpha_inc_tip-alpha_inc_root)/(Number_of_stations-1))*row+alpha_inc_root-((alpha_inc_tip-alpha_inc_root)/(Number_of_stations-1))))*sin_functions(row, 1);
        end
        if ver == 2
        B(row, 1) = mu(row, 1)*(alpha_geo-alpha_L0(1,2)-epsilon(row,1)-(((alpha_inc_tip-alpha_inc_root)/(Number_of_stations-1))*row+alpha_inc_root-((alpha_inc_tip-alpha_inc_root)/(Number_of_stations-1))))*sin_functions(row, 1);
        end
        if ver == 3
        B(row, 1) = mu(row, 1)*(alpha_geo-alpha_L0(1,1))*sin_functions(row, 1);
        end
        if ver == 4
        B(row, 1) = mu(row, 1)*(alpha_geo-alpha_L0(1,2))*sin_functions(row, 1);
        end
    end
    
    % Inverting the matrix 
    X = A\B;
    
    % Determining CL 
    CL_lifting_line = X(1,1)*pi*AR;
    
    % Wing Loading 
    column = 1;
    odd = 1;
    for row = 1:Number_of_stations
        for spanwise = -1:res:0
            Cl_normal(row,column) = X(row, 1)*sin(odd*acos(spanwise));
            Chord_percentage(ver, column) = -(lambda-1)*(spanwise)+1;
            column = column + 1;
        end
        column = 1;
        odd = odd + 2;
    end
    
    for column = 1:length(arr_size)
        for row = 1:Number_of_stations
            if row == 1  
                Cl_normal_column_add(ver, column) = Cl_normal(row, column);
            else
                Cl_normal_column_add(ver, column) = Cl_normal_column_add(ver, column) + Cl_normal(row, column);
            end
        end
    end
        
    for column = 1:length(arr_size)
            Cl_normal_plot(ver,column) = ((2*(1+lambda)/(pi*X(1,1)))*(1/Chord_percentage(ver,column)))*(Cl_normal_column_add(ver,column));
    end
    
    for row = 1:Number_of_stations
        for column = 1:length(arr_size)
                Cl_normal_plot_comp(row,column) = ((2*(1+lambda)/(pi*X(1,1)))*(1/Chord_percentage(ver,column)))*(Cl_normal(row,column));
        end
    end 
end
y = 0;
for npts = 1:length(arr_size)
Cl_sweep_twist_plot(1,npts) = Cl_normal_plot(1,npts) - (1-(2*y/(0.5*b_span)))*(2*(1-cos(gamma))); % No Flap + Twist + Sweep
Cl_sweep_twist_plot(2,npts) = Cl_normal_plot(2,npts) - (1-(2*y/(0.5*b_span)))*(2*(1-cos(gamma))); % Flap + Twist + Sweep
Cl_sweep_twist_plot(3,npts) = Cl_normal_plot(3,npts) - (1-(2*y/(0.5*b_span)))*(2*(1-cos(gamma))) ; % No Flap + Sweep
Cl_sweep_twist_plot(4,npts) = Cl_normal_plot(4,npts) - (1-(2*y/(0.5*b_span)))*(2*(1-cos(gamma))); % Flap + Sweep
Cl_sweep_twist_plot(5,npts) = Cl_normal_plot(1,npts); % No Flap + Twist
Cl_sweep_twist_plot(6,npts) = Cl_normal_plot(2,npts); % Flap + Twist
Cl_sweep_twist_plot(7,npts) = Cl_normal_plot(3,npts); % No Flap
Cl_sweep_twist_plot(8,npts) = Cl_normal_plot(4,npts); % Flap

Cl_3D(1,npts) = Cl_sweep_twist_plot(1,npts)*CL(1,1); % No Flap + Twist + Sweep
Cl_3D(2,npts) = Cl_sweep_twist_plot(2,npts)*CL(1,2); % Flap + Twist + Sweep
Cl_3D(3,npts) = Cl_sweep_twist_plot(3,npts)*CL(1,1); % No Flap + Sweep
Cl_3D(4,npts) = Cl_sweep_twist_plot(4,npts)*CL(1,2); % Flap + Sweep
Cl_3D(5,npts) = Cl_sweep_twist_plot(5,npts)*CL(1,1); % No Flap + Twist
Cl_3D(6,npts) = Cl_sweep_twist_plot(6,npts)*CL(1,2); % Flap + Twist
Cl_3D(7,npts) = Cl_sweep_twist_plot(7,npts)*CL(1,1); % No Flap
Cl_3D(8,npts) = Cl_sweep_twist_plot(8,npts)*CL(1,2); % Flap
y = y + res;
end
spanwise = [-1:res:0];
subplot(2,3,3)
plot(spanwise, Cl_sweep_twist_plot(1,:), 'k')
hold on
plot(spanwise, Cl_sweep_twist_plot(3,:), '-.g')
plot(spanwise, Cl_sweep_twist_plot(5,:), '-.b')
plot(spanwise, Cl_sweep_twist_plot(7,:), '-.r')

plot(-spanwise, Cl_sweep_twist_plot(1,:), 'k')
plot(-spanwise, Cl_sweep_twist_plot(3,:), '-.g')
plot(-spanwise, Cl_sweep_twist_plot(5,:), '-.b')
plot(-spanwise, Cl_sweep_twist_plot(7,:), '-.r')

% plot(spanwise, Cl_normal_plot_comp)
% plot(-spanwise, Cl_normal_plot_comp)
legend('No Flap + Twist + Sweep','No Flap + Sweep','No Flap + Twist', 'No Flap','Interpreter','latex')
title('Monoplane Wing Loading')
xlabel('$\frac{y}{(\frac{b}{2})}$ [\%]','Interpreter','latex')
ylabel('$\frac{C_{l}}{C_{L}}$','Interpreter','latex')
axis([-1.2 1.2 0 2])

% CL vs Span
max_Cl_no_flap = max(Cl_3D(1,:));
max_Cl_flap = max(Cl_3D(2,:));
spanwise = [-1:res:0];
subplot(2,3,4)
plot(spanwise, Cl_3D(1,:), 'k') % No Flap + Twist + Sweep
hold on
plot(spanwise, Cl_3D(2,:), '-.k') % Flap + Twist + Sweep
plot(spanwise, Cl_3D(3,:), 'g') % No Flap + Sweep
plot(spanwise, Cl_3D(4,:), '-.g') % Flap + Sweep
plot(spanwise, Cl_3D(5,:), 'b') % No Flap + Twist
plot(spanwise, Cl_3D(6,:), '-.b') % Flap + Twist
plot(spanwise, Cl_3D(7,:), 'r') % No Flap
plot(spanwise, Cl_3D(8,:), '-.r') % Flap

xline(-1+find(Cl_3D(1,:)==max_Cl_no_flap)*res,'m')
xline(-1+find(Cl_3D(2,:)==max_Cl_flap)*res,'-.m')
xline(1-find(Cl_3D(1,:)==max_Cl_no_flap)*res,'m')
xline(1-find(Cl_3D(2,:)==max_Cl_flap)*res,'-.m')

plot(-spanwise, Cl_3D(1,:), 'k') % No Flap + Twist + Sweep
plot(-spanwise, Cl_3D(2,:), '-.k') % Flap + Twist + Sweep
plot(-spanwise, Cl_3D(3,:), 'g') % No Flap + Sweep
plot(-spanwise, Cl_3D(4,:), '-.g') % Flap + Sweep
plot(-spanwise, Cl_3D(5,:), 'b') % No Flap + Twist
plot(-spanwise, Cl_3D(6,:), '-.b') % Flap + Twist
plot(-spanwise, Cl_3D(7,:), 'r') % No Flap
plot(-spanwise, Cl_3D(8,:), '-.r') % Flap

legend('No Flap + Twist + Sweep','Flap + Twist + Sweep', 'No Flap + Sweep','Flap + Sweep', 'No Flap + Twist', 'Flap + Twist', 'No Flap', 'Flap','Stall Point - No Flap', 'Stall Point - Flap', 'Interpreter','latex')
title('C_{L} vs Span')
xlabel('$\frac{y}{(\frac{b}{2})}$ [\%]','Interpreter','latex')
ylabel('C_{l} [#]')
axis([-1.2 1.2 0 3])

% Span Efficiency
delta = 0;
for npts = length(X)
delta = delta + npts*(X(npts,1)/X(1,1))^2;
end
e = (1+delta)^-1;

% Induced Drag
CD_i = (CL_lifting_line^2)/(pi*e*AR);

% Total Drag
CD = CD_0 + CD_i;