"""3D renderings of the tube_coil.jl runs: rows = initial shapes, then the final shapes for each
confinement and η τE; columns = 4 values of α (−1, −1/3, 1/3, 1). The centerline is drawn as a
shaded tube, the hard wall as a translucent cylinder.

  --view persp   oblique perspective view (rod axis horizontal)
  --view above   from above in perspective, az −15°, el 82°
  --colour alpha     each rod in the colour of its α (vik)
  --colour k2tau     coloured along the rod by κ²τ          } managua_br, ±the 95th percentile of
  --colour qmeasure  coloured along the rod by the Q measure } |value| over the final shapes shown
  --overlay      two columns: the α < 0 shapes drawn on top of each other (translucent), and the α > 0
                 shapes; Wr of each shape in its colour
  --at L/L0      take the snapshot at this L/L0 as the final shape
  --pct          values as % of max |value| over the final shapes shown

usage: python plot_3d.py <hard_folder> <free_folder> [options]
  perspective view:  python plot_3d.py H F --view persp
  overlay:           python plot_3d.py H F --view above --overlay --at 3
  chirality renders: python plot_3d.py H F --view above --colour k2tau --at 3
                     python plot_3d.py H F --view above --colour qmeasure --pct --at 3 --etas 10
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import tubeplots as tp  # noqa: E402

p = argparse.ArgumentParser()
p.add_argument("hard"); p.add_argument("free")
p.add_argument("--view", default="persp", choices=list(tp.VIEWS))
p.add_argument("--colour", default="alpha", choices=["alpha", "k2tau", "qmeasure"])
p.add_argument("--overlay", action="store_true")
p.add_argument("--opacity", type=float, default=0.55, help="opacity of the overlaid shapes")
p.add_argument("--at", type=float, default=None, help="final shape at L/L0 = AT")
p.add_argument("--etas", default=None, help="comma-separated η τE to show (default all)")
p.add_argument("--alphas", default="-1,-0.3333,0.3333,1")
p.add_argument("--confs", default="hard,free")
p.add_argument("--pct", action="store_true")
p.add_argument("--out", default=None)
a = p.parse_args()

data = {"hard": tp.load_runs(a.hard, a.at), "free": tp.load_runs(a.free, a.at)}
tp.orient(*data.values())
confs = a.confs.split(",")
etas = [float(e) for e in a.etas.split(",")] if a.etas else [e for e in tp.ETAS if any(r["eta"] == e for r in data["hard"])]
shown = tp.nearest_alphas(data["hard"], [float(v) for v in a.alphas.split(",")])
R_wall = data["hard"][0]["tube_radius"]
cmap = tp.vik()
run = lambda c, eta, al: next((r for r in data[c] if r["eta"] == eta and r["a"] == al), None)
finals = [(c, eta, al, run(c, eta, al)) for c in confs for eta in etas for al in shown]
finals = [f for f in finals if f[3] is not None]
ref = {al: run("hard", etas[0], al) or run("free", etas[0], al) for al in shown}   # initial shapes

# geometry of the panels
ext = {c: 1.05 * max([np.abs(r["xsaves"][[0, -1], :, :2]).max() for cc, _, _, r in finals if cc == c] + [0.2]) for c in confs}
if "hard" in ext:
    ext["hard"] = max(ext["hard"], 1.15 * R_wall)
zall = np.concatenate([r["xsaves"][[0, -1], :, 2].ravel() for *_, r in finals])
zlim = (min(0.0, zall.min()) - 0.05, max(1.0, zall.max()) + 0.05)

# colour scale along the rod
if a.colour != "alpha":
    measure, clabel = tp.MEASURES[a.colour]
    vals = {(c, eta, al): measure(r["xsaves"][-1]) for c, eta, al, r in finals}
    vals.update({("init", al): measure(r["xsaves"][0]) for al, r in ref.items()})
    scale = 1.0
    if a.pct:
        scale = 100.0 / max(np.abs(v).max() for k, v in vals.items() if k[0] != "init")
        print(f"max |{a.colour}| = {1 / scale * 100:.3g} (= 100%)")
        clabel = clabel.rsplit("(", 1)[0].rstrip() + "\n(% of max)"
    vals = {k: v * scale for k, v in vals.items()}
    vmax = float(np.percentile(np.abs(np.concatenate([v for k, v in vals.items() if k[0] != "init"])), 95))
    vnorm, vcmap = plt.Normalize(-vmax, vmax), tp.managua_br()
    print(f"{a.colour} colour range ±{vmax:.3g} (95th percentile of the final shapes shown)")

groups = [g for g in ([al for al in shown if al < 0], [al for al in shown if al > 0]) if g] if a.overlay else [[al] for al in shown]
opacity = a.opacity if a.overlay else 1.0
rows = [("init", None)] + [(c, eta) for c in confs for eta in etas]
Lf = ref[shown[0]]["Lrest"]
F = plt.figure(figsize=(1.9 * len(groups) + 1.6, 1.9 * len(rows) + 0.6))
gs = F.add_gridspec(len(rows), len(groups) + 1, width_ratios=[0.6] + [1] * len(groups), hspace=0.08, wspace=0.0)
for i, (c, eta) in enumerate(rows):
    lab = F.add_subplot(gs[i, 0]); lab.axis("off")
    lab.text(1.0, 0.5, f"initial, $L/L_0$ = 1" if c == "init" else
             ("hard-walled tube" if c == "hard" else "no confinement") + f", {tp.eta_label(eta)}\n"
             + rf"$L/L_0$ = {Lf[-1] / Lf[0]:.2f}", ha="right", va="center", fontsize=7, transform=lab.transAxes)
    for k, group in enumerate(groups):
        ax = F.add_subplot(gs[i, k + 1], projection="3d")
        wr_txt = []
        for jg, al in enumerate(group):
            r = ref[al] if c == "init" else run(c, eta, al)
            if r is None:
                continue
            xs = r["xsaves"][0 if c == "init" else -1]
            wall = R_wall if (c in ("init", "hard") and "hard" in confs and jg == 0) else None   # wall once per panel
            e = ext[c] if c != "init" else ext.get("hard", ext[confs[0]])
            col = tp.alpha_colour(cmap, al)
            if a.colour == "alpha":
                tp.draw_rod(ax, xs, a.view, zlim, e, color=col, R_wall=wall, opacity=opacity)
            else:
                tp.draw_rod(ax, xs, a.view, zlim, e, values=vals[("init", al) if c == "init" else (c, eta, al)],
                            vnorm=vnorm, vcmap=vcmap, R_wall=wall, opacity=opacity)
            if c != "init":
                wr_txt.append((f"{r['Wr'][-1]:+.2f}", col))
        if c == "init":
            ax.set_title(r"$\alpha$ = " + ", ".join(f"{al:+.2f}" for al in group), fontsize=8, pad=0)
        for jj, (txt, col) in enumerate(wr_txt):    # one Wr per shape, in its colour when overlaid
            n = len(wr_txt)
            ax.text2D(0.5 + (jj - (n - 1) / 2) * 0.22, 0.95, ("Wr = " if n == 1 else "") + txt, ha="center",
                      fontsize=6, color="k" if n == 1 else col, transform=ax.transAxes)
if a.colour != "alpha":
    cax = F.add_axes([0.93, 0.88, 0.012, 0.07])
    F.colorbar(plt.cm.ScalarMappable(cmap=vcmap, norm=vnorm), cax=cax, extend="both")
    cax.set_title(clabel, fontsize=7); cax.tick_params(labelsize=6)
out = a.out or (os.path.basename(os.path.normpath(a.hard)) + "_free" + (f"_at{a.at:g}" if a.at else "")
                + f"_3d_{a.view}" + ("" if a.colour == "alpha" else f"_{a.colour}" + ("_pct" if a.pct else ""))
                + ("_overlay" if a.overlay else "") + (f"_eta{a.etas.replace(',', '-')}" if a.etas else "") + ".pdf")
tp.savefig(F, out)
