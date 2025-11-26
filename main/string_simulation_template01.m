function string_simulation_template01(string_params, V0, tspan, record_status, file_name)
    
    %load string_params into rate function
    my_rate_func = @(t_in, V_in) string_rate_func01(t_in, V_in, string_params);

    %run the integration
    [tlist, Vlist] = ode45(my_rate_func, tspan, V0);
    
    % visualization
    xlist = linspace(0, string_params.L, string_params.n + 2);
    figure();
    hold on;
    axis equal; axis square;
    axis([xlist(1), xlist(end), -0.5, 0.5]);
    xlabel('x (m)'); ylabel('y (m)');

    % adjustments in order to create a real time animation
    fps = 60; % desired frames-per-second
    % uniformly spaced animation times
    t_anim = (tspan(1) : 1/fps : tspan(2)); 
    V_anim = interp1(tlist, Vlist, t_anim, 'linear'); % linear interp
    dt_real = 1/fps; % seconds per frame

    % initialize video
    if record_status == true
        myVideo = VideoWriter(file_name); %open video file
        myVideo.FrameRate = 60;
        open(myVideo)
    end

    for i = 1:length(t_anim)
        title(sprintf('Mass-String System (t = %.2f s)', t_anim(i)));

        % extract interpolated state
        U = V_anim(i, 1:string_params.n);

        U_edited = [0, U, string_params.Uf_func(t_anim(i))];

        plot(xlist, U_edited, '.-', 'LineWidth', 2, 'MarkerSize', 20);

        drawnow;

        if record_status == true
            frame = getframe(gcf); %get frame
            writeVideo(myVideo, frame);
        end

        pause(dt_real); % keeps playback consistent with physical time

        cla;

    end

    if record_status == true
        close(myVideo);
    end
    
    clc;

end