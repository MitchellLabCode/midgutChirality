Data for Fig. 2: shape kymographs and deformation measures of the framed midgut
==============================================================================

Files in this folder
--------------------
  fig2_kymographs.h5            per embryo: radius R(s,t), contact radius R_contact(s,t),
                                edge growth rate, fold and anterior/posterior split positions
  fig2_curvature_measures.csv   curvature measures E_b1, E_b2, E_b, E_tw     (Sec. IV, Eqs. S16-S18)
  fig2_rotation_rates.csv       average rotation speeds omega_b, omega_tw   (Sec. IV, Eq. S15)
  fig2_curvature_change.csv     average curvature change mu_b, mu_tw         (Sec. IV, Eq. S19)
  readme.txt                    this file

Everything is derived from all_framesInterpV6_10samples.mat; its structure
(centerline x, frame Q, Darboux vector kappa, spin omega, meshes, folds, ...)
is described in the Fig. 1 data readme (data_figshare_fig1/readme.txt).
Symbols follow Sec. IV ("Framed curves") of the Supplementary Information.


Samples, genotypes and time alignment
-------------------------------------
10 embryos imaged every 2 min; odd sample numbers are byn>Myo1C, even are WT.
Each embryo is placed on a common time axis by its time lag (integer frames;
from cross-correlation of radius kymographs, with t = 0 at the onset of the
2nd fold of the reference WT embryo, sample 4). The lag of each sample is
stored in the h5 attributes (time_lag_frames, time_lag_min) and is the same
as in the Fig. 1 data:

  sample        1    2    3    4    5    6    7    8    9   10
  genotype     Myo  WT  Myo  WT  Myo  WT  Myo  WT  Myo  WT
  lag [frames] -16   -5  -21   -9   +1  -26   +9   -6  -18  -32

Each quantity lives either on frames or on intervals between frames; the
time column (time_min) gives its aligned time in minutes, using the same
convention as the corresponding Fig. 2 panel (f = 1-based frame index):
  curvature measures E        t = 2 (f - 1)   + 2 lag      (frame f)
  rotation rates omega        t = 2 (f - 1)   + 2 lag      (interval f -> f+1, at its start)
  curvature change mu         t = 2 (f - 1/2) + 2 lag      (interval f -> f+1, at its midpoint)
  radius R, contact radius    t = 2 f         + 2 lag      (frame f)
  edge growth rate            t = 2 (f + 1/2) + 2 lag      (interval f -> f+1, at its midpoint)


Regions
-------
Every deformation measure is reported for three regions, as in the
"Anterior only", "Posterior only" and "All midgut" rows of Fig. S8:
  region = anterior    centerline from the anterior end to the point of
                       maximal x (the apex of the first loop; this split moves
                       with time; max_x_s_over_L in the h5 file)
  region = posterior   from that point to the posterior end
  region = all_midgut  the whole centerline


Deformation measures (the three csv files)
------------------------------------------
Columns common to all three: time_min, sample_id (1-10), is_WT (1 = WT,
0 = byn>Myo1C), region. One row per sample, region and timepoint.

Notation: kappa^i = kappa . d_i and omega^i = omega . d_i are the components
of the Darboux vector and of the spin vector in the material frame (d1 = L-R,
d2 = D-V, d3 = tangent); 1, 2 are bending, 3 is twisting. All integrals run
over the region; "mean" = integral / region length.

  fig2_curvature_measures.csv  (on frames)                                   [units]
    E_b1        mean of (kappa^1)^2 / 2                                     [um^-2]
    E_b2        mean of (kappa^2)^2 / 2                                     [um^-2]
    E_b         E_b1 + E_b2  (bending)                                      [um^-2]
    E_tw        mean of (kappa^3)^2 / 2  (twisting)                         [um^-2]
    E_*_total   the same, integrated instead of averaged over the region,
                e.g. E_b_total = integral (1/2)((kappa^1)^2 + (kappa^2)^2) ds  [um^-1]
                These are the "Total energies" plotted in Fig. S8.

  fig2_rotation_rates.csv  (on intervals; rad/min)
    omega_b     mean of sqrt((omega^1)^2 + (omega^2)^2)  (bending rotation speed)
    omega_tw    mean of |omega^3|                        (twisting rotation speed)
    omega_1     mean of |omega^1|
    omega_2     mean of |omega^2|

  fig2_curvature_change.csv  (on intervals; rad / (um min))
    mu_b        mean of sqrt((d_t kappa^1)^2 + (d_t kappa^2)^2)
    mu_tw       mean of |d_t kappa^3|
    mu_1        mean of |d_t kappa^1|
    mu_2        mean of |d_t kappa^2|

Numerics, identical for all measures and matching the Fig. 2 panels:
  - omega and d_t kappa are frame-to-frame rotations / differences divided by
    the 2-min frame interval.
  - a 3x3 median filter is applied in the (time, arclength) plane; the first and
    last time row and arclength column, where the filter is undefined, are
    dropped (hence the end points of the centerline are not included).
  - integrals are weighted by the local length: Voronoi lengths (average of the
    two adjacent edge lengths) for kappa and d_t kappa, which live on vertices;
    edge lengths for omega, which lives on edges.

Relation to the Supplement:
  - the bending columns (omega_b, mu_b) use the Euclidean norm of the two bending
    components, sqrt(x1^2 + x2^2); Eq. S15 writes omega_b with |omega^1| + |omega^2|.
    omega_1, omega_2 (and mu_1, mu_2) are provided so that either can be formed.
  - Eqs. S16-S18 define E per unit length (the E_b1 ... E_tw columns); Fig. S8 plots
    the integrated values (E_*_total), whose unit is um^-1.


fig2_kymographs.h5
------------------
Root attributes: description, data_file, frame_interval_min (2), created.
One group per embryo, /sample_01 ... /sample_10, with attributes
  sample_id, is_WT (uint8), genotype, time_lag_frames, time_lag_min, n_frames, source_dir
and datasets (MATLAB shapes; see the note on orientation below)

  R                     (n_t x 100)  mean radius of each surface-mesh hoop [um]
  R_s_over_L            (n_t x 100)  arclength fraction of each hoop along the hoop
                                     centroids, per frame (hoops are not evenly spaced)
  R_time_min            (n_t)        aligned time of each row of R and R_contact [min]

  R_contact             (n_t x 100)  contact radius at each centerline point [um]:
                                     1.2 x the distance from the centerline point to
                                     the nearest surface vertex (the radius of the
                                     sphere, centred on the centerline, that touches
                                     the surface, inflated by 20%; used as an excluded
                                     volume / confinement in the rod simulations)
  cent_s_over_L         (100)        arclength fraction of each centerline point

  edge_growth_rate      (n_t-1 x 99) rate of change of each centerline edge length,
                                     (l_k(t+1) - l_k(t)) / 2 min  [um/min]
  edge_s_over_L         (99)         arclength fraction of each edge midpoint
  edge_growth_time_min  (n_t-1)      aligned time of each row [min]

  fold_hoop_index       (n_t x 3)    mesh hoop (1-100) of folds 1-3; NaN before onset
  fold_s_over_L         (n_t x 3)    position of folds 1-3 as arclength fraction along
                                     the centerline (TubULAR fold arclength); NaN before onset
  fold_onset_frame      (3)          frame at which folds 1, 2, 3 form
  max_x_s_over_L        (n_t)        position of the anterior/posterior split (max-x point)

Orientation: MATLAB writes HDF5 column-major, so in Python (h5py / numpy) all
2D datasets appear transposed, e.g. R has shape (100, n_t) -> use R.T for
(time, position). In MATLAB, h5read returns the shapes listed above.

  import h5py
  with h5py.File('fig2_kymographs.h5') as f:
      g = f['sample_04']
      R  = g['R'][()].T                 # (n_t, 100)
      sR = g['R_s_over_L'][()].T        # (n_t, 100)
      t  = g['R_time_min'][()]
      print(g.attrs['is_WT'][0], g.attrs['time_lag_min'][0])

Code: https://github.com/MitchellLabCode/chirality_processing
(processing/Fig2/export_figshare_fig2.m writes this folder.)
