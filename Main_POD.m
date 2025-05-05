clc; 
clear;
close all; 
set(groot,'defaultLineLineWidth',1.5)
set(groot,'DefaultAxesFontSize',16,'DefaultAxesFontName','Times New Roman')

% diverging red–white–blue colormap
cmap = [0 0 1;  % Blue (Min)
        1 1 1;  % White (Middle)
        1 0 0]; % Red (Max)
num_color = 256; % Number of colors
customColormap = interp1([1, num_color/2, num_color], cmap, linspace(1, num_color, num_color));

%% ----------- 2. Load velocity snapshots & grid information -------------
load Velocity_POD_Re_150.mat            % provides: u_POD v_POD X_cc Y_cc ...

% strip ghost cells only once here if needed
x_c = 0;  y_c = 0;  r_c = 0.5;          % cylinder
X_cc_plot   = X_cc(Ny_start:Ny_end, Nx_start:Nx_end);
Y_cc_plot   = Y_cc(Ny_start:Ny_end, Nx_start:Nx_end);
cylMask     = (X_cc_plot-x_c).^2 + (Y_cc_plot-y_c).^2 <= r_c^2;
cylDisp     = nan(size(cylMask));  cylDisp(cylMask) = 1;

Ny = size(u_POD{1},1);  
Nx = size(u_POD{1},2);
Ns = Ny*Nx;                                     % spatial DOF per field

%% ----------- 3. Assemble POD data matrix (time × space*2) --------------
% convert cell array -> 3‑D array (Ny×Nx×Nt) without a for‑loop
u_stack = cat(3, u_POD{7:end});
v_stack = cat(3, v_POD{7:end});
Nt = size(u_stack,3);

% vectorise each field & place side‑by‑side
Xsnap   = [ reshape(u_stack, Ns, Nt).' , reshape(v_stack, Ns, Nt).' ];  % Nt × 2Ns
Xfluc = Xsnap - mean(Xsnap);

%% ----------- 4. Run POD -------------------------------------------------
[U,S,V]  = svd(Xfluc,'econ');                 % U: modes, S: singular values, V: coeffs
sing   = diag(S);        % vector of singular values
Phi    = V.';            % each ROW is now one spatial POD mode
energy = sing.^2 / sum(sing.^2);
energy   = sing.^2;
eFrac    = energy / sum(energy);
cumE     = cumsum(eFrac);

%% ----------- 5. Visualise first few spatial modes ----------------------
nShow = 4;                       % how many modes you want to display

for k = 1:nShow
    plotMode(Phi, Ns, Ny, Nx, X_cc_plot, Y_cc_plot, cylDisp, customColormap,...
             k, 1, 1);          % u‑component
    plotMode(Phi, Ns, Ny, Nx, X_cc_plot, Y_cc_plot, cylDisp, customColormap,...
             k, 2, 2);          % v‑component
end


set(figure(1),'Position',[100 100 800 600])
set(figure(2),'Position',[100 100 800 600])

% -----------------------------------------------------------------------
function plotMode(U, Ns, Ny, Nx, Xcc, Ycc, cylDisp, cmap, ...
                  modeIdx, fieldIdx, figNo)
% fieldIdx = 1 → u‑part,  2 → v‑part
    figure(figNo); subplot(2,2,modeIdx); hold on

    if fieldIdx == 1
        fld = U(modeIdx, 1:Ns);
        compChar = 'u';
    else
        fld = U(modeIdx, Ns+1 : 2*Ns);
        compChar = 'v';
    end
    
    % draw the black cylinder
    theta  = linspace(0,2*pi,120);
    xc     = 0.5 * cos(theta);          % centre (0,0) and r = 0.5 -> change if needed
    yc     = 0.5 * sin(theta);

    fldM = reshape(fld, Ny, Nx);
    % imagesc(Xcc(1,:), Ycc(:,1), fldM);
    contourf(Xcc, Ycc, fldM, 100,'LineStyle','none');
    patch('XData',xc,'YData',yc,...
          'FaceColor','k','EdgeColor','none');     % solid black disk

    % set(gca,'CLim',[-0.01 0.01]);   % <- re‑apply your limits
    hold off
    axis equal;  
    set(gca,'ydir','normal');
    clim([min(fldM,[],'all')  max(fldM,[],'all')])
    colormap("jet")
    % colormap(cmap); 
    colorbar
    xlabel('x (m)'); 
    ylabel('y (m)')
    xlim([-5,10])
    ylim([-5,5])
    title(sprintf('%s Mode %d', compChar, modeIdx))
end
