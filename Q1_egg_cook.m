% This code serves to calculate the time needed to cook a selection of eggs
% by means of boiling and hot water. The time is caculated by the numerical
% solution of the 1D heat equation in spherical coordinates.

arg_set = ... % Code Name, Text Name, Radius (mm), initial temp, Outer temp, 3D graph div(*@\label{code:q1_args_start}@*)
 {"reg_egg",     "20 deg Chicken Egg",               0.0255, 20, 100,  1000; % 45mm x 57mm
  "fridge_egg",  "5 deg Chicken Egg",                0.0255, 5,  100,  1000;
  "quail_egg",   "20 deg Quail Egg",                 0.0145, 20, 100,   500; % 32mm x 25mm
  "ostrich_egg", "20 deg Ostrich Egg",               0.0875, 20, 100, 10000; % 20cm x 15cm
  "slow_egg",    "20 deg Chicken Egg in warm water", 0.0255, 20, 81,   2500; }; % (*@\label{code:q1_args_end}@*)

size_ = size(arg_set);
len = size_(1);
for j = 1:len
    args = arg_set(j,:);
        
    % Parameter definitions ---------------------------------------------------
    safe_name = args{1};
    log = fopen(sprintf('out/q1/%s/output.log',safe_name),'w');
    human_name = args{2};
    R = args{3}; % Radius of egg [M]
    T_start = args{4}; % Room temperature in celsius
    T_water = args{5}; % Temperature of the cooking water
    graph_div = args{6};
    T_min = 80; % Minimum cooking temp throughout
    t_hold = 10; % Hold for 10 seconds
    k = .500; % conductivity of egg [W K^-1 m^-1]
    rho = 1035; % Density of egg [kg m^-3]
    c_p = 3200; % Specific heat of Egg [J kg^-1 K^-1]
    N = 100; % number of grid points
    dt=0.01; % Size of time step
    dr=R/N;  % grid spacing

    alpha = k / (rho*c_p); % Increment coefficient
    coeff = 1-2*alpha*dt/(dr^2);
    fprintf(log, "alpha = %2.10f, coeff = %2.10f \n", alpha, coeff);

    % Initialization ----------------------------------------------------------
    % solution grid
    x = linspace(0,R,N+1);
    T = ones(N+1,1);
    % Initial Condition (*@\label{code:q1_initial_cond_start}@*)
    T(:,1) = T_start;
    % Boundary conditions
    T(end,1) = T_water; %(*@\label{code:q1_initial_cond_end}@*)

    % PDE Solution ------------------------------------------------------------ (*@\label{code:q1_loop_start}@*)
    k = 1;
    at_temp_time = 0;
    while (at_temp_time < t_hold); % While holding for temp
        T(1,k+1) = max(T(2,k) - (T(3,k) - T(2,k)), T_start); % Set center to 1 'slope' lower
        T(end,k+1) = T(end,k); % Retain end value (*@\label{code:q1_boundary_cond_start}@*) (*@\label{code:q1_boundary_cond_end}@*)
        for i=2:N % Increment over all but the ends (*@\label{code:q1_increment_start}@*)
            r = (i-1) * dr; % Get the radius
            d2T_dr2 = (T(i+1,k)-2*T(i,k)+T(i-1,k))/(dr^2); % Get the instantaneous accel
            dT_dr = (T(i+1,k)-T(i-1,k))/(2*dr); % Get the instaneous slope
            T(i,k+1)=T(i,k) + alpha*dt*(d2T_dr2 + (2/r)*dT_dr); % Increment Temp 
        end % (*@\label{code:q1_increment_end}@*)
        % At every N seconds, Output some temp values.
        time = k * dt;
        if not(mod(time, 10))
            fprintf(log, "Time: %5.6f - Temp: %3.14f, %3.14f, %3.14f, %3.14f, %3.14f, %3.14f, %3.14f \n", ...
            time, T(1, k), T(N/100, k), T(N/10, k), T(N*5/10,k), T(N*9/10,k), T(N*99/100,k), T(N-1,k));
        end
        k = k + 1; % Increment the counter
        % Increment the time if we are at temp
        if all(T(:,k) > T_min)
            at_temp_time = at_temp_time + dt;
        end
    end % (*@\label{code:q1_loop_end}@*)
    
    
    % Final Tally --------------------------------------------------------------
    
    size_t = length(T(1,:));
    size_r = length(T(:,1));
    t_set = linspace(0, time, size_t);
    r_set = linspace(0, R, size_r);
    
    % Cook Time
    time = size_t * dt;
    min = (time - mod(time, 60)) / 60;
    sec = round(time - min*60);
    fprintf(log, "Time to cook %s: %f sec - %im%is\n", human_name, time, min, sec);

    % Calculate heat energy absorbed
    vol = 4 * pi .* r_set.^3 .* dr;
    energy = (T(:,end) - T_start) .* vol' .* c_p .* rho;
    fprintf(log, "Total energy absorbed for %s: %f Joules \n", human_name, sum(energy));

    % Visualization ------------------------------------------------------------
    f = figure('visible','off');
    % line plot
    indexes = [1, round(size_t/4), round(size_t/2), round(size_t*3/4), round(size_t)];
    times = indexes * dt;
    hold on;
    plot(T(:,indexes(1)));
    plot(T(:,indexes(2)));
    plot(T(:,indexes(3)));
    plot(T(:,indexes(4)));
    plot(T(:,indexes(5)));
    xlabel('r (mm)'); ylabel('T(x,t) (C)');
    legend(split(sprintf('t = %i sec,', round(times(:))), ","), 'Location', "east");
    exportgraphics(f, sprintf("out/q1/%s/2D_Plot.png",safe_name), 'Resolution', 300);
    clf(f);

    % surface plot
    [X,Y] = meshgrid(r_set,t_set(1:graph_div:end));
    mesh(X,Y,T(:,1:graph_div:end)');
    colormap('parula');
    xlabel('r (mm)'); ylabel('t (sec)'); zlabel('T(x,t) (C)');
    colorbar;
    caxis([5 100]);
    exportgraphics(f, sprintf("out/q1/%s/3D_Plot.png",safe_name), 'Resolution', 300);
    clf(f);
    fclose(log);
end
exit
