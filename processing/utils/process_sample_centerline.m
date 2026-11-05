function [x_raw, x] = process_sample_centerline(centerline_loc, n_sample, w, N_s)
% makeFramedCurve
% Makes the anatomical frame from the raw gut centerline and mycoyte curves.
% requires smooth_curve3D 
% inputs: 
%  - centerline_loc : full filename of combined centerline file
%  - n_samples: int number of equi-arclength sampled points
%  - w: float for smoothing routine
%  - N_s : number of Gaussian smoothing iterations
% outputs:
%  x_raw : (n_time)x(n_samples)x3 unsmoothed centerline
%  x : (n_time)x(n_samples)x3 smoothed centerline

%%%% Load data and format them

% load centerlines.
listing = dir(fullfile(centerline_loc, 'centerline_from_DVhoops_*.mat'));
num_t = length(listing);

x_raw = zeros(num_t, n_sample, 3);
for i = 1:num_t
    % listing is ascii dictionary order - which should be good
    loader = load(fullfile(centerline_loc,listing(i).name));
    avgpts_ss = loader.avgpts_ss;
    avgpts = loader.avgpts;
    % arclength resample
    sample_ss = linspace(0, avgpts_ss(end), n_sample);
    x_raw(i,:,:) = interp1(avgpts_ss, avgpts,sample_ss);
end

% smooth the centerline in space and time

ell = w;
tau = w;
x = smooth_stdata(x_raw, ell, tau, N_s);

end