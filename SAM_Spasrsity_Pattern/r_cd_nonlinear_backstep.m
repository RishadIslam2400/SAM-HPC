% script to run cd2d_nonlinear

% interior grid points - actual grid (nx+2) x (ny+2)
nx = 10; % increase nx and ny to create larger systems
ny = 10;
n = nx*ny;
u0 = ones(n,1);
c = 0; 
u0 = c*u0; % Initial guess (currenty just zeros, but play with previous 
            % two lines to experiment with different initial guesses)

% This next line will call gmres with no preconditioner - it's commented 
% because convergence is so poor, but you can uncomment and compare against 
% the preconditioned solves

%[iterations, soltime, prectime,fdtime,gmresinfo,backs,sol] = noprec_gmres(nx,ny,n,u0);

[iterations1, soltime1, prectime1,fdtime1,gmresinfo1,back1,sol1] = ILUTPall_gmres(nx,ny,n,u0); % recompute

[iterations2, soltime2, prectime2,fdtime2,gmresinfo2,back2,sol2] = ILUTPonce_gmres(nx,ny,n,u0); % reuse

[iterations3, soltime3, prectime3,fdtime3,gmresinfo3,back3,sol3] = ILUTPSAMall_gmres(nx,ny,n,u0); % approximate
