%% Function to Execute the Analytical Solution.


function T_bound = analytical_solver(f_anal,Nx,Ny)
    
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    N = Nx*Ny;
    x = linspace(hx,1-hx,Nx);     % x coordinates of interior points
    y = linspace(hy,1-hy,Ny);     % y coordinates of interior points
    T = zeros(Ny,Nx);
    for j = 1:Ny
        for i = 1:Nx
            T(j,i) = f_anal(x(i),y(j));
        end
    end
    
    % Add zero-boundaries and plotting
    T_bound = [zeros(1,Nx); T; zeros(1,Nx)];
    T_bound = [zeros(Ny+2,1), T_bound, zeros(Ny+2,1)];
    title = strcat('Analytical solution for Nx = Ny = ', num2str(Nx));
    plotter(T_bound,Nx,Ny,title);
    
end
