
datdir="/Users/npmitchell/Dropbox/Soft_Matter_Biophysics/PAPER/gut_chirality/supplementary_videos/20240531_142410_22C_HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/";
infn = datdir+"Time_000005_WT_masked_embryo_rgb.tif"; 
outfn = datdir + "Time_000005_WT_masked_embryo_rgb-anatomicalFiji.tif";
open(infn);
run("Rotate... ", "angle=13 grid=1 interpolation=None fill enlarge");
run("Reslice [/]...", "output=1.000 start=Top avoid");
run("Rotate... ", "angle=-8 grid=1 interpolation=None fill enlarge");
run("Reslice [/]...", "output=1.000 start=Left avoid");
run("Rotate 90 Degrees Right");
run("Rotate... ", "angle=49 grid=1 interpolation=None fill enlarge");
save(outfn);
