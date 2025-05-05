function [u_cc,v_cc] = Add_BCs_to_cc_velocity(u_cc,v_cc,u_in)

    % Poiseuille Flow
    u_cc(:,1) = 2*u_in-u_cc(:,2);              % BCs at u(0,y)
    u_cc(:,end) = u_cc(:,end-1);               % BCs at u(1,y)
    u_cc(1,:) = -u_cc(2,:);                     % BCs at u(x,0)
    u_cc(end,:) =  -u_cc(end-1,:);              % BCs at u(x,1)

    v_cc(:,1) = -v_cc(:,2);                     % BCs at v(0,y)
    v_cc(:,end) = v_cc(:,end-1);                 % BCs at v(1,y)
    v_cc(1,:) = -v_cc(2,:);                     % BCs at v(x,0)
    v_cc(end,:) = -v_cc(end-1,:);                 % BCs at v(x,1)

% Uniform flow
    % u_cc(:,1) = 2*u_in-u_cc(:,2);                     % BCs at u(0,y)
    % u_cc(:,end) = u_cc(:,end-1);               % BCs at u(1,y)
    % u_cc(1,:) = u_cc(2,:);                     % BCs at u(x,0)
    % u_cc(end,:) =  u_cc(end-1,:);      % BCs at u(x,1)
    % 
    % v_cc(:,1) = -v_cc(:,2);                     % BCs at v(0,y)
    % v_cc(:,end) = v_cc(:,end-1);                 % BCs at v(1,y)
    % v_cc(1,:) = -v_cc(2,:);                     % BCs at v(x,0)
    % v_cc(end,:) = -v_cc(end-1,:);                 % BCs at v(x,1)

% Lid-driven cavity
    % u_cc(:,1) = -u_cc(:,2);                     % BCs at u(0,y)
    % u_cc(:,end) = -u_cc(:,end-1);               % BCs at u(1,y)
    % u_cc(1,:) = -u_cc(2,:);                     % BCs at u(x,0)
    % u_cc(end,:) = 2*u_in - u_cc(end-1,:);      % BCs at u(x,1)
    % 
    % v_cc(:,1) = -v_cc(:,2);                     % BCs at v(0,y)
    % v_cc(:,end) = -v_cc(:,end-1);                 % BCs at v(1,y)
    % v_cc(1,:) = -v_cc(2,:);                     % BCs at v(x,0)
    % v_cc(end,:) = -v_cc(end-1,:);                 % BCs at v(x,1)

end