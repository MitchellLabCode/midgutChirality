Data for Fig. 1: midgut centerline, material frame and chirality of 10 Drosophila embryos
=======================================================================================

Files in this folder
--------------------
  chiral_kymographs.h5         coiling (chirality) kymographs q(s, t), one group per embryo
  wr_tw_10samples_V6.csv       writhe and twist per embryo and timepoint
  readme.txt                   this file

The source dataset all_framesInterpV6_10samples.mat (~240 MB) is described in
section 3; the two files above are derived from it.


Samples, genotypes and time alignment
-------------------------------------
10 embryos, imaged every 2 min. Odd sample numbers are byn>Myo1C, even are WT:
  WT        : HandGFP; byn-GAL4, klar > UAS-mCh-CAAX, H2Av-iRFP   (samples 2, 4, 6, 8, 10)
  byn>Myo1C : HandGFP; byn-GAL4, klar > UAS-Myo1C-RFP, H2Av-iRFP     (samples 1, 3, 5, 7, 9)

All Fig. 1 / Fig. 2 plots use a common, aligned time axis. Each embryo has a
time lag (an integer number of frames) found by cross-correlating radius
kymographs: within each genotype against a reference embryo, plus one
cross-genotype offset, and finally shifted so that t = 0 is the onset of the
2nd fold of the reference WT embryo (sample 4). The aligned time of frame f
(1-based) of sample i is

    t_aligned [min] = 2 * (f - 1) + 2 * time_lag_frames(i)      (kymographs, h5)

The Wr/Tw plots and the csv use t = 2*f + 2*time_lag_frames, i.e. the same
lag with the frame placed one interval (2 min) later.

  sample  genotype    time_lag_frames  time_lag_min  n_frames
     1    byn>Myo1C        -16            -32           76
     2    WT                -5            -10           76
     3    byn>Myo1C        -21            -42           78
     4    WT                -9            -18           74 (75 centerlines)
     5    byn>Myo1C         +1             +2           51
     6    WT               -26            -52           78
     7    byn>Myo1C         +9            +18           53
     8    WT                -6            -12           63
     9    byn>Myo1C        -18            -36           71
    10    WT               -32            -64           91


1. chiral_kymographs.h5
-----------------------
Quasi-local coiling measure q(s, t): at each position s along the centerline,
the signed area (solid angle, steradians) on the unit sphere swept by the unit
tangent over a window of 2w+1 = 11 consecutive tangents, computed with the
spherical-triangle sum of qmeasure.m. Positive and negative q are the two
handednesses of local coiling. The centerline used is the smoothed centerline
`x` of the .mat file (section 3), resampled to equal arclength; the first 9
and last 1 points are skipped.

Root attributes:
  description, data_file, window_half_width_w (5), n_skip_start (9),
  n_skip_end (1), frame_interval_min (2), time_convention, created

One group per embryo, /sample_01 ... /sample_10, with attributes
  sample_id        int     1..10
  is_WT            uint8   1 = WT, 0 = byn>Myo1C
  genotype         string
  time_lag_frames  double  shift used to align this embryo in all Fig1/Fig2 plots [frames]
  time_lag_min     double  same, in minutes (= 2 * time_lag_frames)
  n_frames         int     number of timepoints
  source_dir       string  original TubULAR output folder
and datasets
  q                (n_s x n_frames)  coiling measure [sr]; n_s = 84
  s_over_L         (n_s)             arclength fraction of each row of q (0 = anterior end)
  time_min         (n_frames)        aligned time of each column [min]
  fold_s_over_L    (n_frames x 3)    s/L of folds 1-3; NaN before each fold's onset
  max_x_s_over_L   (n_frames)        s/L of the centerline's max-x point, used to split
                                     anterior / posterior (dashed line in the kymographs)

Array orientation: MATLAB writes HDF5 in column-major order, so in Python
(h5py, numpy) every 2D dataset appears transposed, e.g. q has shape
(n_frames, n_s): q[t, s].

  import h5py
  with h5py.File('chiral_kymographs.h5') as f:
      g = f['sample_04']
      q, s, t = g['q'][()], g['s_over_L'][()], g['time_min'][()]   # q: (n_frames, n_s)
      is_wt, lag = g.attrs['is_WT'][0], g.attrs['time_lag_min'][0]

In MATLAB, h5read returns the (n_s x n_frames) orientation described above.


2. wr_tw_10samples_V6.csv
-------------------------
One row per embryo and timepoint (711 rows), columns
  aligned_time_min  aligned time, t = 2*f + 2*time_lag_frames [min]
  sample_id         1..10
  is_WT             1 = WT, 0 = byn>Myo1C
  Wr                writhe of the centerline `cent` (Klenin & Langowski 2000, method 1a)
  Tw                twist of the material director d1 about `cent` (same paper, method 2a)
Note: sample 6, last timepoint (104 min): Tw drops from 0.94 to 0.08 in one
frame. This is an artifact of the discrete twist (one node wrapping by a full
turn at a very tight anterior bend), not a real change.


3. all_framesInterpV6_10samples.mat  (MATLAB v7.3 / HDF5)
---------------------------------------------------------
Contains one variable, all_data: a 1 x 10 cell array, one struct per embryo
(sample i = all_data{i}). Lengths are in micrometers, in the TubULAR lab
frame; the anterior-posterior axis runs roughly along +x (anterior at low x).
n_t is the number of timepoints of that embryo (see the table above). Index
order is always (time, position along the gut, ...). The per-sample fields:

  dir           string            TubULAR output folder of the raw data.

  x_raw         (n_t, 105, 3)     Gut centerline: centroids of the TubULAR mesh hoops,
                                  resampled to 105 points equally spaced in arclength.
  x             (n_t, 105, 3)     x_raw smoothed in space (Taubin: lambda = 0.3,
                                  mu = -0.32, 20 iterations) and time (diffusive,
                                  tau = 0.3, 20 iterations).
  cent          (n_t, 100, 3)     x with its first 5 (anterior-most) points removed;
                                  the centerline on which the frame is built.

  myocytes_raw  struct            The four longitudinal myocyte stripes, traced on the
                                  surface, each (n_t, 100, 3), resampled to 100 points:
                                    red = ventral-left (VL), green = dorsal-left (DL),
                                    yellow = dorsal-right (DR), blue = ventral-right (VR).
  myocytes      struct            Same stripes, smoothed as x.

  Q             (n_t, 99, 3, 3)   Material frame on each edge k of cent (between points
                                  k and k+1): Q(t, k, i, :) is director d_i,
                                    d1 = Q(:,:,1,:)  left -> right,
                                    d2 = Q(:,:,2,:)  ventral -> dorsal,
                                    d3 = Q(:,:,3,:)  unit tangent, anterior -> posterior.
                                  d1 is built from the four stripes (where each stripe
                                  crosses the cross-section at that point, matched by
                                  arclength), made orthogonal to d3; d2 = d3 x d1.
                                  Orthonormal to 1e-16.
  kappa         (n_t, 98, 3)      Darboux vector at the interior points of cent, in the
                                  lab frame [rad/um]: rotation from frame k to frame k+1,
                                  log(Q_k' Q_{k+1}), divided by the local (Voronoi)
                                  segment length. Components in the material frame,
                                  kappa_i = dot(kappa, d_i): kappa_1, kappa_2 bending
                                  about d1 (dorso-ventral bending) and d2 (left-right
                                  bending), kappa_3 twisting.
  omega         (n_t-1, 99, 3)    Frame spin: rotation of the frame on each edge from
                                  timepoint t to t+1, lab frame [rad per frame;
                                  divide by 2 for rad/min]. Components in the material
                                  frame as for kappa.

  meshes        (n_t, 1) cell     Gut surface at each timepoint, struct with
                                    v (10000 x 3) vertices [um], f (19602 x 3) faces
                                    (1-based vertex indices).
                                  The vertices form a 100 x 100 grid:
                                  reshape(v, [100, 100, 3]) is (hoop, position around
                                  the hoop, xyz); hoop 1 is anterior, hoop 100 posterior.
  radii         (n_t, 100)        Mean radius of each mesh hoop [um].

  folds         (n_t, 3) int16    Mesh-hoop index (1..100) of folds 1-3 at each
                                  timepoint. Meaningful only after each fold's onset;
                                  before it, the value is held constant at the hoop
                                  where the fold will form.
  ssfold        (n_t, 3)          Arclength position of folds 1-3 along the TubULAR
                                  centerline [um]; ssfold ./ ssmax is the fold
                                  position as a fraction of gut length.
  ssmax         (n_t, 1)          Total length of the TubULAR centerline [um].
  fold_onset    (1, 3)            Frame index at which folds 1, 2 and 3 form.
  timeInterval  scalar            Frame interval [min] (2).

Known small inconsistencies:
  - sample 4: x and x_raw have 75 timepoints, all other fields 74.
  - sample 10: folds, ssfold and ssmax have 90 rows, all other fields 91.

Python access (h5py): each element of all_data is an HDF5 object reference;
arrays come out transposed (reversed axis order) relative to the MATLAB
shapes above. In MATLAB: all_data = load('all_framesInterpV6_10samples.mat').all_data.
