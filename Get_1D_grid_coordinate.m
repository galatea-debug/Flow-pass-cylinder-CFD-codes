function [x,y] = Get_1D_grid_coordinate(xmin,xmax,ymin,ymax,x_uni_min,x_uni_max,y_uni_min,y_uni_max,Ncore,r)

% Generate uniform part of mesh
x_uni = linspace(x_uni_min,x_uni_max,Ncore+1);
y_uni = linspace(y_uni_min,y_uni_max,Ncore+1);
dx_uni = x_uni(2) - x_uni(1);
dy_uni = y_uni(2) - y_uni(1);

% Generate 1D interior cell faces for Left half: x
Mx_left = floor(log(1-(x_uni_min-xmin)*(1-r)/dx_uni)/log(r)); % Number of cells in first half
delta_x_left = (x_uni_min - xmin) * (1 - r) / (1 - r^Mx_left); % Initial cell width
% delta_x_left = dx_uni;
x_left = zeros(1, Mx_left); % M cells -> M+1 faces
x_left(1) = xmin;

for i = 1:Mx_left
    x_left(i+1) = x_left(i) + delta_x_left * r^(Mx_left-i);
end

% Generate 1D interior cell faces for Left half: y
My_down = floor(log(1-(y_uni_min-ymin)*(1-r)/dy_uni)/log(r)); % Number of cells in first half
delta_y_down = (y_uni_min - ymin) * (1 - r) / (1 - r^My_down); % Initial cell width
% delta_y_down = dy_uni;
y_down = zeros(1, My_down); % M cells -> M+1 faces
y_down(1) = ymin;

for i = 1:My_down
    y_down(i+1) = y_down(i) + delta_y_down * r^(My_down-i);
end

% Generate 1D interior cell faces for Right half: x
Mx_right = floor(log(1-(xmax-x_uni_max)*(1-r)/dx_uni)/log(r)); % Number of cells in first half
delta_x_right = (xmax - x_uni_max) * (1 - r) / (1 - r^Mx_right); % Initial cell width
% delta_x_right = dx_uni;
x_right = zeros(1, Mx_right); % M cells -> M+1 faces
x_right(1) = x_uni_max;

for i = 1:Mx_right
    x_right(i+1) = x_right(i) + delta_x_right * r^(i-1);
end

% Generate 1D interior cell faces for Right half: y
My_up = floor(log(1-(ymax-y_uni_max)*(1-r)/dy_uni)/log(r)); % Number of cells in first half
delta_y_up = (ymax - y_uni_max) * (1 - r) / (1 - r^My_up); % Initial cell width
% delta_y_up = dy_uni;
y_up = zeros(1, My_up); % M cells -> M+1 faces
y_up(1) = y_uni_max;

for i = 1:My_up
    y_up(i+1) = y_up(i) + delta_y_up * r^(i-1);
end

dx_left = x_left(2) - x_left(1);
dy_left = y_down(2) - y_down(1);
dx_right = x_right(end) - x_right(end-1);
dy_right = y_up(end) - y_up(end-1);

x_left_node = xmin - dx_left;
y_left_node = ymin - dy_left;
x_right_node = xmax + dx_right;
y_right_node = ymax + dy_right;

x = [x_left_node,x_left(1:end-1),x_uni,x_right(2:end),x_right_node];
y = [y_left_node,y_down(1:end-1),y_uni,y_up(2:end),y_right_node];

end

