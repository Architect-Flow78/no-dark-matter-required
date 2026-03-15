# No Dark Matter Required

**Effective Gravitational Dimension D(R) Transitions to 2 in 116 of 171 SPARC Disk Galaxies**

Nicolae Pascal, 2026

## Key Result

We analyze 171 galaxy rotation curves from the SPARC database and find that **68% show effective gravitational dimension D → 2 at outer radii** (median = 1.87). This explains flat rotation curves as a geometric effect of disk structure, without dark matter.

The baryonic Tully-Fisher exponent of **4** follows analytically from the D: 3→2 transition.

![All 171 galaxies](results/fig1_all_D_curves.png)

## How It Works

One formula: **D(R) = 2(1 - α)**, where α = d(log V)/d(log R)

- D = 3 → Keplerian decline (spherical geometry)
- D = 2 → Flat rotation curve (disk geometry)
- D = 0 → Solid body rotation

No fitting. No free parameters. Just observed velocities and radii.

## Quick Start

```bash
# 1. Clone this repo
git clone https://github.com/YOUR_USERNAME/no-dark-matter-required.git
cd no-dark-matter-required

# 2. Download SPARC data
wget https://zenodo.org/records/16284118/files/Rotmod_LTG.zip
unzip Rotmod_LTG.zip -d sparc_data

# 3. Run analysis
pip install numpy matplotlib
python sparc_dimension_analysis.py
```

## Paper

Full paper with derivations: [Zenodo DOI](https://doi.org/10.5281/zenodo.XXXXXXX)

## Results Summary

| Metric | Value |
|--------|-------|
| Galaxies analyzed | 171 |
| D_outer in [1.5, 2.5] | 116 (68%) |
| Median D_outer | 1.87 |
| BTFR exponent (predicted) | 4.00 |
| BTFR exponent (measured) | 4.09 |
| a_0 predicted (cH_0/2pi) | 1.1 × 10⁻¹⁰ m/s² |
| a_0 observed | 1.2 × 10⁻¹⁰ m/s² |

## License

MIT License. Data from SPARC (Lelli, McGaugh, Schombert 2016).
