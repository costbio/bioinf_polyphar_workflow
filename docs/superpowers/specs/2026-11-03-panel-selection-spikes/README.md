# Panel-selection investigation spikes (2026-11-03)

Throwaway scripts behind the evidence in
`../2026-11-03-panel-selection-redesign-design.md`. They are **not pipeline code**: no
tests, hardcoded assumptions, run by hand against a v2 run folder on POSEIDON.
Kept so the numbers in the design can be reproduced.

Run order (each needs the previous one's cached JSON in `~/scratch_clust`):

| script | what it established |
|---|---|
| `cluster_spike.py` | Louvain on the unfiltered Pmax>=0.7 network: purity 0.19, no families |
| `cluster_spike2.py` | druggability/quality filters do not rescue clustering (internal Pmax flat) |
| `cluster_spike3.py` | size-conditional threshold on size *disparity* — WRONG, admits ~1% overlap pairs |
| `cluster_spike4.py` | corrected size-conditional threshold on absolute size; per-bin table; still purity 0.19 |
| `plddt.py` | pocket pLDDT does not detect AF mega-cavity artifacts (84.8 vs 83.3) |
| `span.py` | artifacts are oversized pockets: 90 residues / 7 segments vs 16 / 3 |
| `panel_grid.py` | TM annotation fetch; drugScore/volume are inert as panel filters |
| `panel_final.py` | volume-ceiling x segment-cap sweep; the 241-protein panel (pre-subpocket fix) |
| `tm_variants.py` | TM policies A/B/C; policy C removes only 2 proteins |
| `filtered_clust.py` | **clustering works on the filter-eligible set**; subpocket hole found |
| `families.py` | the 7 clusters at Pmax>=0.7 and panel coverage; corrected 209-protein panel |
| `mw_check.py` | eDrug3D heavy-atom distribution; HCC drugs at the 64th-88th percentile |
| `build_calib.py` | builds the 48-pocket / 250-ligand calibration run folder |
| `ligsize.py` | eDrug3D size/flexibility profile; the `rot>15 OR MW>900` dockability filter |

`build_calib.py` is the only one that writes anything outside `~/scratch_clust`.
