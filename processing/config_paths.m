function P = config_paths()
% config_paths  All paths and per-sample settings of the processing pipeline.
%
%   P = config_paths();
%
% Edit this file (only) to point the pipeline at your copies of the raw
% data, the TubULAR code and the intermediate data. Every script in
% processing/ calls it and adds processing/utils to the path.
%
% Sample numbering: samples 1-10 in the order below. Odd samples are
% byn>Myo1C, even samples are WT. This numbering is used everywhere,
% including the time alignment.

here = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(here, 'utils')));

%%% code
P.tubular = '/mnt/data/code/nromeo/tubular/tubular';   % TubULAR (github.com/npmitchell/tubular)

%%% raw data: one TubULAR project per embryo (the folder holding xp.mat and tubular_output/)
raw = '/mnt/data/Chris';
myo = 'HandGFPbynGAL4klar_UASMyo1CRFPHiRFP';
wt  = 'HandGFPbynGAL4klar_UASmChCAAXH2AviRFP';
P.samples = struct( ...
  'id',       num2cell(1:10), ...
  'is_WT',    num2cell(logical(mod(1:10, 2) == 0)), ...
  'raw_dir',  { ...
    fullfile(raw, myo, '20240528_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/Deconvolved_2024-05-28_150855/All_timepoint_channels/deconvolved_16bit_0810'), ...
    fullfile(raw, wt,  '20240527_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/Deconvolved/All_timepoint_channels/deconvolved_16bit_0810'), ...
    fullfile(raw, myo, '2024-06-26_175116_crisp_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/deconvolved32bit/deconvolved_16bit_0810'), ...
    fullfile(raw, wt,  '20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/deconvolved/deconvolved_16bit_0810'), ...
    fullfile(raw, myo, '2024-06-28_180134_good_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpacked_flipped_rotated_TIFFs_new/deconvolved/deconvolved_16bit_0810'), ...
    fullfile(raw, wt,  '2024-06-07_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/deconvolved32bit/deconvolved_16bit_0810'), ...
    fullfile(raw, myo, '2024-06-29_232230_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/deconvolved32bit/deconvolved_16bit_0810'), ...
    fullfile(raw, wt,  '2024-07-01_3part_HandGFPbynGAL4klar_UASmCHCAAXHiRFP/deconvolved32bit/deconvolved_16bit_0810'), ...
    fullfile(raw, myo, '2024-05-29_182916_22C_excellent_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/deconvolved32bit/deconvolved_16bit_0810'), ...
    fullfile(raw, wt,  '2024-07-03_111720_HandGFPbynGAL4klar_UASmCHCAAXHiRFP/deconvolved32bit/deconvolved_16bit_0810')}, ...
  ... % number of timepoints used when exporting meshes (tubi_to_meshrs) and in the
  ... % time alignment (timelag_fromRadii); [] = all. These are the crops the
  ... % published data were made with; the two steps differ for samples 4 and 10.
  'n_t_meshes',    {[], [], [], 74, [], [], [], 63, 71, 91}, ...
  'n_t_timealign', {[], [], [], [], [], [], [], 63, 71, 90}, ...
  ... % sample 8 was re-processed into a separate grid folder
  'spcutmesh_override', {'', '', '', '', '', '', '', ...
    fullfile(raw, wt, '2024-07-01_3part_HandGFPbynGAL4klar_UASmCHCAAXHiRFP/deconvolved32bit/deconvolved_16bit_0810/tubular_output/gridCoords_nU0100_nV0100_forNico/sphi_cutMesh_000step/spcutMesh_%06d.mat'), ...
    '', ''});
for k = 1:numel(P.samples)
    P.samples(k).tubular_output = fullfile(P.samples(k).raw_dir, 'tubular_output');
    if P.samples(k).is_WT, P.samples(k).genotype = 'WT'; else, P.samples(k).genotype = 'byn>Myo1C'; end
end

%%% intermediate data (not in the repository)
P.data       = '/mnt/data/code/nromeo/data';                     % framed-curve dataset
P.framesV6   = fullfile(P.data, 'all_framesInterpV6_10samples.mat');
P.meshes     = '/mnt/data/code/nromeo/meshes_spRS';              % meshes_sample<N>.mat (tubi_to_meshrs)
P.timeshifts = fullfile(here, 'timeAlign', 'timeshifts_fromRadiusKymos_10samples.mat');

%%% inputs of individual figures
% stripe data of the embryo shown in the q-measure illustration (Fig. 1)
P.qmeasure_example = '/mnt/data/code/nromeo/curve_functions/stare_tests/muscleNucleiStripes.mat';
% selected Cosserat runs (output of cosserat/prod/make_selected.sh), for fig2_wr_twist_sim
P.sims = '/mnt/data/code/nromeo/Fig2/selected_20260922';

%%% outputs
P.out_fig1       = fullfile(here, 'Fig1', 'output');
P.out_fig2       = fullfile(here, 'Fig2', 'output');
P.figshare_fig1  = fullfile(here, 'Fig1', 'data_figshare_fig1');
P.figshare_fig2  = fullfile(here, 'Fig2', 'data_figshare_fig2');
P.bounding_contact = fullfile(here, '..', 'cosserat', 'data');    % refIV6_sample<N>_walls3_cball.mat
P.canalization   = fullfile(here, '..', 'cosserat', 'canalization');  % SRVbox3f_s<WT>s<Myo>/ (prep_interp_data)

%%% constants
P.dt = 2;                 % frame interval [min]
P.time_convention = 'published';   % see aligned_time.m

%%% optional machine-specific overrides: processing/config_paths_local.m
%%% (not tracked), a function P = config_paths_local(P)
if exist(fullfile(here, 'config_paths_local.m'), 'file')
    addpath(here);
    P = config_paths_local(P);
end
end
