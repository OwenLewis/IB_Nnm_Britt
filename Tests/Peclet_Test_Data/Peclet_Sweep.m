Pe = 1.2e3;
Refine = [1,2,4,8];

%some path management. 
addpath("../")
addpath("../../src/")

% These parameters set the diffusive timescale
% & total simulation time
D = 2.5e-3;
L = 3;
Tmax = 50;

%Now calculate advection velocity based on Peclet
v = Pe*D/L;

%This is the cfl I was running my old tests at. 
cfl = 0.1*0.1*64/3;

%Lowest resolution grid
baseN = 64;



for i = 1:length(Refine)


    N = baseN*Refine(i)
    dx = L/N

    dt = cfl*dx/v

    [u,t] = Diffusion_translate_single(N,v,D,Tmax,dt,0,0,0);


    filename = sprintf("diff_translate_N%i_Pe%0.1e.mat",N,Pe)
    save(filename,'u','t')
end

rmpath("../")
rmpath("../../src/")
