%This is the script to generate data for the test of advection-diffusion in
%a circle. We are advecting along the coordinate axis (x in this case), at
%speed v = 0.1. The time-step is chosen to exactly enforce CFL = 1. In this
%case, our simple upwinding advection scheme 'should' be perfect. This
%should allow us to compare with the stationary circle, and exactly analyze
%the errors introduced by the movement of the boundary (i.e. uncovering new
%cells with the extended solution. 

%One thing to think about here, the diffusion coefficient is very small, so
%this is a high peclet number test. Does that matter? Unclear. Should
%definitely return to this. 


refine = 1
[u,t] = Diffusion_translate_single(64*refine,0.1,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/064_adv_axis_align.mat','u','t','-v7.3')
refine = 2
[u,t] = Diffusion_translate_single(64*refine,0.1,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/128_adv_axis_align.mat','u','t','-v7.3')
refine = 4
[u,t] = Diffusion_translate_single(64*refine,0.1,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/256_adv_axis_align.mat','u','t','-v7.3')
refine = 8
[u,t] = Diffusion_translate_single(64*refine,0.1,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/512_adv_axis_align.mat','u','t','-v7.3')
refine = 16
[u,t] = Diffusion_translate_single(64*refine,0.1,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/1024_adv_axis_align.mat','u','t','-v7.3')



%Now, lets run some stationary domain diffusion problems with the exact
%same parameters for comparison. 

refine = 1
[u,t] = Diffusion_circle_single(64*refine,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/064_diff_shorttime.mat','u','t','-v7.3')
refine = 2
[u,t] = Diffusion_circle_single(64*refine,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/128_diff_shorttime.mat','u','t','-v7.3')
refine = 4
[u,t] = Diffusion_circle_single(64*refine,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/256_diff_shorttime.mat','u','t','-v7.3')
refine = 8
[u,t] = Diffusion_circle_single(64*refine,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/512_diff_shorttime.mat','u','t','-v7.3')
refine = 16
[u,t] = Diffusion_circle_single(64*refine,2.5e-3,5,((4/64)/0.1)/refine,1,0,0);
save('../TimeSplitting/1024_diff_shorttime.mat','u','t','-v7.3')