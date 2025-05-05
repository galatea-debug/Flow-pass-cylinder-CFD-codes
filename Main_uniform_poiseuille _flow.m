clc
clear
close all

set(groot,'defaultLineLineWidth',1.5)
set(0,'DefaultAxesFontSize',16);
set(0,'DefaultAxesFontName','Times New Roman');
addpath("../Functions/")

%% Initialization
% Domain size and basic informations
xmin = 0;
xmax = 20;
ymin = 0;
ymax = 1;
x_c = 0;
y_c = 0;
r_c = 0.5;

dt = 0.01;
t_end = 100;
max_iter_burgers = 1000;
max_iter_PPE = 100;
tol_burgers = 1e-10;
tol_PPE = 1e-10;
tol_steady = 1e-8;

% Viscosity and Reynolds number
u_in = 1;
D_c = ymax - ymin;
nu = 1/150;
Re = u_in*D_c/nu;           % Reynolds number

%% Mesh grid
% Region of uniform mesh (Cylinder region)
Nx = 100;
Ny = 100;
r = 1.01;

dx1 = (xmax-xmin)/(Nx);
dy1 = (ymax-ymin)/(Ny);
x = linspace(xmin-dx1,xmax+dx1,Nx+3);
y = linspace(ymin-dy1,ymax+dy1,Ny+3);

% Generate 1D interior cell faces for first half (dense near x=0)
M = Ny/2; % Number of cells in first half
delta_x0 = 0.5 * (1 - r) / (1 - r^M); % Initial cell width
y_half = zeros(1, M+1); % M cells -> M+1 faces
for i = 1:M
    y_half(i+1) = y_half(i) + delta_x0 * r^(i-1);
end

% Create interior 1D faces (symmetric, dense near x=0 and x=1)
% y_interior = zeros(1, Ny+1);
% y_interior(1:M+1) = y_half;
% y_interior(M+2:Ny+1) = 1 - fliplr(y_half(1:M)); % Mirror for second half
% 
% y = [-y_interior(2),y_interior,ymax+(ymax-y_interior(end-1))];
% x = linspace(xmin-dx1,xmax+dx1,Nx + 3);

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
% Cylinder_mask = (X_cc-x_c).^2 + (Y_cc-y_c).^2 <= r_c^2;

Cylinder_mask = zeros(Ny+2,Nx+2); 
North = zeros(Ny+2,Nx+2); 
South = zeros(Ny+2,Nx+2); 
West = zeros(Ny+2,Nx+2); 
East = zeros(Ny+2,Nx+2);

% for i1 = 2:Nx+1
%     for j1 = 2:Ny+1
%         if(Cylinder_mask(j1,i1) == 0)
%             if (Cylinder_mask(j1+1,i1))
%                 North(j1,i1) = 1;
%             end
% 
%             if (Cylinder_mask(j1-1,i1))
%                 South(j1,i1) = 1;
%             end
% 
%             if (Cylinder_mask(j1,i1-1))
%                 West(j1,i1) = 1;
%             end
% 
%             if (Cylinder_mask(j1,i1+1))
%                 East(j1,i1) = 1;
%             end
%         end
%     end
% end

color = struct();
color.cylinder = Cylinder_mask;
color.East = East;
color.West = West;
color.North = North;
color.South = South;

% Plot grid
figure(1)
set(gcf,'Position',[100 100 800 600])
plot(X_edge,Y_edge,'k'); hold on
plot(X_edge.',Y_edge.','k');hold on
scatter(X_cc,Y_cc,'.','MarkerEdgeColor','k','MarkerFaceColor','k');hold on
% surface(X_cc, Y_cc, zeros(size(Cylinder_display)), Cylinder_display, ...
%         'FaceColor', 'flat', 'EdgeColor', 'none', ...
%         'FaceAlpha', 1, 'CDataMapping', 'direct');
axis equal
hold off
xlabel('x')
ylabel('y')
title(['Grid plots: Uniform flow'])

%% Initial conditions and BCs
% Initial conditions: cell center
randomamp = 0.01;
u_cc = zeros(size(X_cc));
v_cc = zeros(size(Y_cc));
% u_cc(2:end-1,2:end-1) = 1 + randomamp.*rand(Ny,Nx);
% v_cc(2:end-1,2:end-1) = randomamp.*rand(Ny,Nx);
[u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_in);

% u_cc(Cylinder_mask) = 0;
% v_cc(Cylinder_mask) = 0;

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

%% Time marching
t = 0;
Nt = t_end/dt;
Res_steady = zeros(Nt,1);

for n = 1:Nt
    t = t + dt; % Update time
    tic
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
    tic
    [p_cc,Iter_count,RMS_residual] = PPE_solver_GS(p_cc,RHS_p,coeff,color,Nx,Ny,max_iter_PPE,tol_PPE,dx);
    
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
    ptime=toc;
    fprintf("TS number=%d, (u,v)i=%d, Alltime=%d, conv=%d, div=%d \n",...
            n,uv_iter,ptime,conv, err_in_div );

    if (conv < tol_steady)
        break;
    end
end


%% Contour plot
cmap = [0 0 1;  % Blue (Min)
        1 1 1;  % White (Middle)
        1 0 0]; % Red (Max)
num_color = 256; % Number of colors
customColormap = interp1([1, num_color/2, num_color], cmap, linspace(1, num_color, num_color));

figure(2)
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['u: Re = ',num2str(Re)],'FontWeight','bold')
xlim([9,11])
ylim([ymin,ymax])
% clim([0.5,1.5])

figure(3)
contourf(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),50,'LineStyle','none')
colorbar
colormap(customColormap)
xlabel('X (m)')
ylabel('Y (m)')
% zlabel('u')
axis equal
title(['v: Re = ',num2str(Re)],'FontWeight','bold')
xlim([9,11])
ylim([ymin,ymax])
% clim([-0.5,0.5])

figure(4);
l = streamslice(X_cc(2:end-1,2:end-1),Y_cc(2:end-1,2:end-1),u_cc(2:end-1,2:end-1),v_cc(2:end-1,2:end-1),10); hold on
hold off
axis equal
set(l,'LineWidth',1)
set(l, 'Color', "k")
% set(l, 'Color', "#0072BD")
title(['Streamline: Re = ',num2str(Re)],'FontWeight','bold')
xlim([xmin,xmax])
ylim([ymin,ymax])

%% u,v Profile at half
Nx_half = floor(Nx/2) + 0;
Ny_half = floor(Ny/2) + 0;

u_half = u_cc(:,Nx_half);
v_half = v_cc(Ny_half,:);
x_half = X_cc(Ny_half,:);
y_half = Y_cc(:,Nx_half);

u_analytical = 6*u_in.*y_half.*(D_c - y_half)./D_c^2;
% u_analytical = ones(length(y_half),1);

 % u
 figure(5)
 plot(u_half,y_half,'linestyle','none','Marker','o','color','k','DisplayName','Numerical solution'); hold on
 plot(u_analytical,y_half,'linestyle','-','color','k','DisplayName','Analytical solution'); hold on
 hold off
 xlabel('$ u \ (m/s) $','Interpreter','latex')
 ylabel('$ y \ (m) $','Interpreter','latex')
 grid on
 legend('Location','southeast')
 xlim([0,1.5])
 ylim([0,1])
 title(['Poiseuille flow (Re = ',num2str(Re),'): \it{u} \rm{profile}'])
 drawnow
