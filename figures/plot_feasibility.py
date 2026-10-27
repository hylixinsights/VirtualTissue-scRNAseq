"""Preliminary lineage-count screening; no state-outcome visualization."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
d = json.loads((ROOT / "reports/evidence/cohort-audit.json").read_text())
fig, (flow, ax) = plt.subplots(1, 2, figsize=(12, 5.2), gridspec_kw={"width_ratios": [1, 1.2]})
fig.patch.set_facecolor("white")
flow.set_xlim(0, 1); flow.set_ylim(0, 1); flow.axis("off")
boxes = [
    (0.75, f"Clinical cohort\n{d['clinical_patients']} patients · {d['clinical_biopsies']} biopsies"),
    (0.44, f"Deposited immune object\n{d['object_patients']} patients · {d['object_biopsies']} biopsies\n{d['object_cells']:,} cells"),
    (0.13, "Anti-PD-1 monotherapy\nNivolumab or pembrolizumab\nCD8 totals only; no state changes"),
]
for y, label in boxes:
    flow.add_patch(FancyBboxPatch((0.05, y), .9, .19, boxstyle="round,pad=0.015", facecolor="#EEF3F7", edgecolor="#748494"))
    flow.text(.5, y+.095, label, ha="center", va="center", fontsize=11, color="#1C344B")
for y in [.75, .44]:
    flow.annotate("", xy=(.5, y-.10), xytext=(.5, y-.02), arrowprops={"arrowstyle":"->", "color":"#748494"})
keys = ["30", "15", "50"]
values = [d["count_eligible_upper_bounds"][k] for k in keys]
ax.barh([2, 1, 0], values, color=["#AC4436", "#487D95", "#8A949B"], height=.48)
ax.set_yticks([2,1,0], ["30 CD8 / biopsy\nPrimary", "15 CD8 / biopsy\nSensitivity", "50 CD8 / biopsy\nSensitivity"])
for y, value in zip([2,1,0], values): ax.text(value+.25,y,str(value),va="center",fontsize=12)
ax.axvline(14, color="#AC4436", linewidth=1.5, linestyle="--")
ax.text(13.8, 2.48, "Needed for G1:\n10 development + 4 lockbox", ha="right", va="bottom", fontsize=10, color="#AC4436")
ax.set_xlim(0,16); ax.set_ylim(-.5,3.1); ax.set_xticks([0,3,6,10,14])
ax.set_xlabel("Maximum patients passing CD8 counts")
ax.spines[["top","right","left"]].set_visible(False)
ax.tick_params(axis="y",length=0)
fig.suptitle("SPECIAL cannot meet the current cohort-size requirement", fontsize=15, x=.05, ha="left")
fig.text(.05,.025,"Preliminary upper bounds: intervals remain unavailable and pairing is unresolved. No lockbox has been assigned.",fontsize=10,color="#555555")
fig.tight_layout(rect=[.02,.08,1,.93])
fig.savefig(ROOT / "figures/S1_flow_preliminary.svg", bbox_inches="tight")
fig.savefig(ROOT / "figures/S1_flow_preliminary.png", dpi=160, bbox_inches="tight")
