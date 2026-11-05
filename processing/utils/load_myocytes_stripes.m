function [myocytes_raw, myocytes] = load_myocytes_stripes(folder, num_sample, w, N_s)
%%% load_myocyte_stripes
% load the myocytes from muscleNucleiStripes.mat files
%
% Input:
% folder : path to folder containing muscleNucleiStripes.mat
%  - num_sample : int - number of desired points on the myocyte curves
%  - w: float for smoothing routine
%  - N_s : number of Gaussian smoothing iterations
% Output:
% myocytes : struct containing all the myocytes at all available times

% color to position map: Red = VL; Green =  DL ; Yellow = DR ; Blue = VR

loader = load(fullfile(folder, 'muscleNucleiStripes.mat'));
dat = loader.dat;
n_times = length(dat);
times  = zeros(1, n_times);

myocytes_raw = struct;
myocytes_raw.red = zeros(n_times,num_sample, 3);
myocytes_raw.green = zeros(n_times,num_sample, 3);
myocytes_raw.yellow = zeros(n_times,num_sample, 3);
myocytes_raw.blue = zeros(n_times,num_sample, 3);

s_interp = linspace(0,1,num_sample);

for t = 1:n_times
    cdat = dat{t};
    % resample each curve to num_sample points

    % Red
    
    myocytes_raw.red(t,:,:) = interp1(linspace(0,1,size(cdat.curv1,1)), cdat.curv1, s_interp);
    % Green
    myocytes_raw.green(t,:,:) = interp1(linspace(0,1,size(cdat.curv2,1)), cdat.curv2, s_interp);
    % Yellow
    myocytes_raw.yellow(t,:,:) = interp1(linspace(0,1,size(cdat.curv3,1)), cdat.curv3, s_interp);
    % Blue
    myocytes_raw.blue(t,:,:) = interp1(linspace(0,1,size(cdat.curv4,1)), cdat.curv4, s_interp);
end

% Space-time smoothing
myocytes = struct;
myocytes.red = smooth_stdata(myocytes_raw.red, w,w,N_s);
myocytes.green = smooth_stdata(myocytes_raw.green, w,w,N_s);
myocytes.yellow = smooth_stdata(myocytes_raw.yellow, w,w,N_s);
myocytes.blue = smooth_stdata(myocytes_raw.blue, w,w,N_s);

end
