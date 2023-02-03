% INPUTS
% T: Temperature values, including boundaries
% what: Title for plots


function plotter(T,Nx,Ny,what)
    x = linspace(0,1,Nx+2);
    y = linspace(0,1,Ny+2);

    % Create Countour Plot
    figure();
    [C,h] = contourf(x,y,T);
    colormap(jet);
    clabel(C,h);
    xlabel("X axis"); ylabel("Y axis");
    title("Contour plot: ",what);

    % Create surface plot
    figure();
    surf(x,y,T);
    xlabel("X axis"); ylabel("Y axis");
    title("Surface Plot: ",what);
end
    
