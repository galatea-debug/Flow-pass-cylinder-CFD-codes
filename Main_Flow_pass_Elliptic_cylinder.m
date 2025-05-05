clc
clear
close all

set(groot,'defaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times New Roman');
addpath("..\Functions\")

%% Initialization
% Domain size and basic informations
Scale_factor = 5;
xmin = -5*Scale_factor;
xmax = 10*Scale_factor;
ymin = -5*Scale_factor;
ymax = 5*Scale_factor;
x_c = 0;
y_c = 0;
a = 0.5;
b = 1/6;
AOA = -30;

dt = 0.0025;
t_end = 40;
max_iter_burgers = 100;
max_iter_PPE = 100;
tol_burgers = 1e-10;
tol_PPE = 1e-10;
tol_steady = 1e-8;

% Viscosity and Reynolds number
u_in = 1;
nu = 1/300;
L = 2*a;
Re = u_in*L/nu;           % Reynolds number

% colormap
cmap = [0 0 1;  % Blue (Min)
        1 1 1;  % White (Middle)
        1 0 0]; % Red (Max)
num_color = 256; % Number of colors
customColormap = interp1([1, num_color/2, num_color], cmap, linspace(1, num_color, num_color));

%% Mesh grid
% Region of uniform mesh (Cylinder region)
x_uni_min = -1;
x_uni_max = 1;
y_uni_min = -1;
y_uni_max = 1;
Ncore = 200;
r = 1.03;

[x,y] = Get_1D_grid_coordinate(xmin,xmax,ymin,ymax,x_uni_min,x_uni_max,y_uni_min,y_uni_max,Ncore,r);

Nx_uni_min = find(x == x_uni_min);
Nx_uni_max = find(x == x_uni_max);
Ny_uni_min = find(y == y_uni_min);
Ny_uni_max = find(y == y_uni_max);

% dx1 = (xmax-xmin)/(Ncore);
% dy1 = (ymax-ymin)/(Ncore);
% x = linspace(xmin-dx1,xmax+dx1,Ncore+3);
% y = linspace(ymin-dy1,ymax+dy1,Ncore+3);

[X_edge,Y_edge] = meshgrid(x,y);
Nx = length(x)-3;
Ny = length(y)-3;

[X_cc,Y_cc,X_u_fc,Y_u_fc,X_v_fc,Y_v_fc,Delta_x,Delta_y,dx,dy] = Calculate_coordinate(X_edge,Y_edge,Nx,Ny);
dx_min=min(min(dx,[],"all"),min(Delta_x,[],"all"));
dy_min=min(min(dy,[],"all"),min(Delta_y,[],"all"));
dx_min=min(dx_min,dy_min);
% Delta_x: Distance between edge point in x dir
% Delta_y: Distance between edge point in y dir
% dx: Distance between face center point (v_face) in x dir
% dy: Distance between face center point (v_face) in y dir

%% Create cylinder color matrix
% Shift coordinates
X_shift = X_cc - x_c;
Y_shift = Y_cc - y_c;

% Rotated coordinates
theta_rad = deg2rad(AOA);
X_rot =  X_shift * cos(theta_rad) + Y_shift * sin(theta_rad);
Y_rot = -X_shift * sin(theta_rad) + Y_shift * cos(theta_rad);

Cylinder_mask = (X_rot./a).^2 + (Y_rot./b).^2 <= 1;

% Plot black-filled cylinder using the mask
Cylinder_display = nan(size(Cylinder_mask));  % Initialize with NaNs
Cylinder_display(Cylinder_mask) = 1;          % Mark cylinder interior

% Cylinder_mask = zeros(Ny+2,Nx+2); 
North = zeros(Ny+2,Nx+2); 
South = zeros(Ny+2,Nx+2); 
West = zeros(Ny+2,Nx+2); 
East = zeros(Ny+2,Nx+2);

for i1 = 2:Nx+1
    for j1 = 2:Ny+1
        if(Cylinder_mask(j1,i1) == 0)
            if (Cylinder_mask(j1+1,i1))
                North(j1,i1) = 1;
            end

            if (Cylinder_mask(j1-1,i1))
                South(j1,i1) = 1;
            end

            if (Cylinder_mask(j1,i1-1))
                West(j1,i1) = 1;
            end

            if (Cylinder_mask(j1,i1+1))
                East(j1,i1) = 1;
            end
        end
    end
end

color = struct();
color.cylinder = Cylinder_mask;
color.East = East;
color.West = West;
color.North = North;
color.South = South;

% Plot grid
figure(1)
plot(X_edge,Y_edge,'k'); hold on
plot(X_edge.',Y_edge.','k');hold on
scatter(X_cc,Y_cc,'x','filled','MarkerEdgeColor','k','MarkerFaceColor','k');hold on
surface(X_cc, Y_cc, zeros(size(Cylinder_display)), Cylinder_display, ...
        'FaceColor', 'flat', 'EdgeColor', 'none', ...
        'FaceAlpha', 1, 'CDataMapping', 'direct'); hold on
hold off
xlabel('x')
ylabel('y')
title('Grid plots')

%% Initial conditions and BCs
% Initial conditions: cell center
randomamp = 0.01;
u_cc = ones(size(X_cc));
v_cc = zeros(size(Y_cc));
u_cc(2:end-1,2:end-1) = 1 + randomamp.*rand(Ny,Nx);
v_cc(2:end-1,2:end-1) = randomamp.*rand(Ny,Nx);
[u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_in);

u_cc(Cylinder_mask) = 0;
v_cc(Cylinder_mask) = 0;

p_cc = zeros(size(X_cc));

[U_fc,~] = Get_face_variable(u_cc,Delta_x,Delta_y,Nx,Ny,color,"D"); 
[~,V_fc] = Get_face_variable(v_cc,Delta_x,Delta_y,Nx,Ny,color,"D"); 

u_fc = U_fc;                 % Average u on face center
v_fc = V_fc;                 % Average v on face center
U_fc_old = U_fc;
V_fc_old = V_fc;
u_cc_old=u_cc;
v_cc_old=v_cc;

% Calculate coefficient for diffusion term
coeff = calculate_diffusion_coefficient(dx,dy,Delta_x,Delta_y,Nx,Ny);

% Get Sparse matrix for PPE
Corelation = Corelation_Matrix_array(Nx,Ny);
A_Matrix = Get_sparse_matrix(coeff, color, Corelation, Nx, Ny);

%% Time marching
t = 0;
Nt = t_end/dt;
Res_steady = zeros(Nt,1);

% Force
Drag = zeros(Nt,1);
Lift = zeros(Nt,1);

% Save velocity for POD
% t_POD = [10,20,40,60,80,100];
% N_POD = t_POD./dt;
% k1 = 1;
% u_POD = cell(1,length(t_POD));
% v_POD = cell(1,length(t_POD));

% Plot Velocity
t_plot = [0.1,1,10,20,30,40];
N_plot = t_plot./dt;
k2 = 1;

for n = 1:Nt
    tic
    t = t + dt; % Update time

    % Burgers solver
    [u_cc_star,v_cc_star,uv_iter] = Burgers_solver(u_cc,v_cc,U_fc,V_fc,u_cc_old, v_cc_old, U_fc_old,V_fc_old, ...
                                    Delta_x,Delta_y,coeff,Nx,Ny,dt,Re,max_iter_burgers,tol_burgers,u_in,color);

    % Calculate divergence
    [U_fc_star,~] = Get_face_variable(u_cc_star,Delta_x,Delta_y,Nx,Ny,color,"D"); 
    [~,V_fc_star] = Get_face_variable(v_cc_star,Delta_x,Delta_y,Nx,Ny,color,"D"); 
    
    % Adjust net flow
    Q_in = U_fc_star(:,1).*Delta_y(1:end-1,1);
    Q_out = U_fc_star(:,end).*Delta_y(1:end-1,end);
    dQ = Q_in - Q_out;
    Q_out_modify = Q_out + dQ;
    U_fc_star(:,end) = Q_out_modify./Delta_y(1:end-1,end);

    Div_U_star = calculate_divergence(U_fc_star,V_fc_star,Delta_x,Delta_y,Nx,Ny,color);
    %sum(abs(Div_U_star).*Delta_x.*Delta_y,"all")

    % PPE solver
    RHS_p = Div_U_star./dt;
    % [p_cc,Iter_count,RMS_residual] = PPE_solver_GS(p_cc,RHS_p,coeff,color,Nx,Ny,max_iter_PPE,tol_PPE);      % GS solver
    p_cc = PPE_solver_Matlab_solver(p_cc,RHS_p,A_Matrix,Corelation,color,Nx,Ny,max_iter_PPE,tol_PPE);       % Matlab solver
    [Grad_p_cc_x,Grad_p_cc_y] = calculate_gradient_p_cc(p_cc,Delta_x,Delta_y,Nx,Ny,color);
    [Grad_p_fc_x,Grad_p_fc_y] = calculate_gradient_p_fc(p_cc,dx,dy,Nx,Ny,color);

    % Update velocity
    u_cc_new = u_cc_star - dt.*Grad_p_cc_x;
    v_cc_new = v_cc_star - dt.*Grad_p_cc_y;
    U_fc_new = U_fc_star - dt.*Grad_p_fc_x;
    V_fc_new = V_fc_star - dt.*Grad_p_fc_y;
    
    U_fc_old = U_fc;
    V_fc_old = V_fc;
    u_cc_old = u_cc;
    v_cc_old = v_cc;

    u_cc = u_cc_new;
    v_cc = v_cc_new;
    U_fc = U_fc_new;
    V_fc = V_fc_new;

    % Calculate force
    Drag_p = color.East(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max).*p_cc(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max) - ...
             color.West(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max).*p_cc(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max);
    Lift_p = color.North(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max).*p_cc(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max) - ...
             color.South(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max).*p_cc(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max);
    Drag(n) = sum(Drag_p.*Delta_y(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max),"all");
    Lift(n) = sum(Lift_p.*Delta_x(Ny_uni_min:Ny_uni_max,Nx_uni_min:Nx_uni_max),"all");

    % Calculate steady state
    diff_u=(u_cc-u_cc_old);
    diff_v=(v_cc-v_cc_old);
    conv_u=rms(diff_u,"all");
    conv_v=rms(diff_v,"all");
    conv=max(conv_u,conv_v);
    Res_steady(n) = conv;

    Div_U_new = calculate_divergence(U_fc_new,V_fc_new,Delta_x,Delta_y,Nx,Ny,color);
    err_in_div=max(abs(Div_U_new.*Delta_x.*Delta_y),[],"all");
    err_in_div_star = sum(abs((1-color.cylinder).*Div_U_star.*Delta_x.*Delta_y),"all");
    alltime=toc;
    fprintf("TS number=%d, (u,v)i=%d, All_time=%d, conv=%d, div=%d \n",...
            n,uv_iter,alltime,conv, err_in_div );

    % if(n == N_POD(k1))
    %     u_POD{k1} = u_cc;
    %     v_POD{k1} = v_cc;
    %     k1 = k1 + 1;
    % end

    if (n == N_plot(k2))
        Plot_variable(u_cc,v_cc,p_cc,X_cc,Y_cc,Nx,Ny,Delta_x,Delta_y,color,k2,t,customColormap,Re,xmin,xmax,ymin,ymax,Cylinder_display);
        k2 = k2 + 1;
    end

end

%% Force coefficient
Lref = 2*sqrt(a^2*cos(AOA)^2+b^2+sin(AOA)^2);
C_D = 2.*Drag./(u_in^2*Lref);
C_L = 2.*Lift./(u_in^2*Lref);

% filename_POD = "Velocity_POD_Re_" + num2str(Re) + ".mat";
filename_Force = "Elliptic_cylinder_with_AOA_Force_coefficient_Re_" + num2str(Re) + ".mat";
% save(filename_POD,"u_POD","v_POD");
save(filename_Force,"C_D","C_L");

%% Plots
u_cc_plot = u_cc;
v_cc_plot = v_cc;
p_cc_plot = p_cc;
u_cc_plot(color.cylinder) = nan;
v_cc_plot(color.cylinder) = nan;
p_cc_plot(color.cylinder) = nan;

w = calculate_vorticity(u_cc,v_cc,Nx,Ny,Delta_x,Delta_y,color);
w_plot = w;
w_plot(color.cylinder) = nan;

xmin_plot = xmin + 10;
xmax_plot = xmax - 20;
ymin_plot = ymin + 10;
ymax_plot = ymax - 10;

figure
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc_plot(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['u: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
xlim([xmin_plot,xmax_plot])
ylim([ymin_plot,ymax_plot])

figure
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),v_cc_plot(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['v: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
xlim([xmin_plot,xmax_plot])
ylim([ymin_plot,ymax_plot])

figure
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),w_plot(2:end-1,2:end-1),100,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['Vorticity: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
clim([-1,1])
xlim([xmin_plot,xmax_plot])
ylim([ymin_plot,ymax_plot])

figure
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),p_cc_plot(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['Pressure: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
xlim([xmin_plot,xmax_plot])
ylim([ymin_plot,ymax_plot])

figure
l = streamslice(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),10); hold on
surface(X_cc, Y_cc, zeros(size(Cylinder_display)), Cylinder_display, ...
    'FaceColor', 'flat', 'EdgeColor', 'none', ...
    'FaceAlpha', 1, 'CDataMapping', 'direct');
% colormap([0 0 0])  % Set colormap to black
% caxis([0 1])       % Set color range so 1 maps to black, NaN to transparent
hold off
axis equal
set(l,'LineWidth',1)
set(l, 'Color', "k")
axis equal
title(['Streamline: Re = ',num2str(Re), ' (t = ', num2str(t),' s)'],'FontWeight','bold')
xlim([xmin_plot,xmax_plot])
ylim([ymin_plot,ymax_plot])


%% Force plot
t = dt:dt:t_end;

% Force coefficient
C_D = 2.*Drag./(u_in^2*L);
C_L = 2.*Lift./(u_in^2*L);

figure(7)
plot(t,C_D,'-'); hold on
hold off
xlabel('t (s)')
ylabel('C_{D}')
grid on
xlim([1,t_end])

figure(8)
plot(t,C_L,'--'); hold on
hold off
xlabel('t (s)')
ylabel('C_{L}')
grid on
xlim([1,t_end])
