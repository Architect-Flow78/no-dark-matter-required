"""
===========================================================================
No Dark Matter Required:
Effective Gravitational Dimension D(R) in SPARC Disk Galaxies

Author: Nicolae Pascal
Date: March 2026
Paper: https://doi.org/10.5281/zenodo.XXXXXXX  (update after upload)
===========================================================================

This script analyzes rotation curves of 171 disk galaxies from the SPARC
database and extracts the effective gravitational dimension D(R).

Key result: 68% of galaxies show D -> 2 at outer radii, explaining flat
rotation curves without dark matter.

USAGE:
  1. Download SPARC data: https://zenodo.org/records/16284118/files/Rotmod_LTG.zip
  2. Unzip into a folder called 'sparc_data/'
  3. Run: python sparc_dimension_analysis.py

REQUIREMENTS:
  pip install numpy matplotlib
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import sys

# ============================================================
# CONFIGURATION
# ============================================================
SPARC_DIR = "sparc_data"  # folder with *_rotmod.dat files
OUTPUT_DIR = "results"
MIN_POINTS = 5            # minimum data points per galaxy

# ============================================================
# CORE FUNCTION: Extract effective dimension from rotation curve
# ============================================================
def compute_D(R, V):
    """
    Compute effective gravitational dimension D(R).
    
    D(R) = 2(1 - alpha), where alpha = d(log V)/d(log R)
    
    D = 3: Keplerian (spherical geometry, V ~ R^{-1/2})
    D = 2: Flat curve (disk geometry, V = const)
    D = 0: Solid body rotation (V ~ R)
    
    Parameters
    ----------
    R : array - radii (any units)
    V : array - rotation velocities (any units)
    
    Returns
    -------
    D : array - effective dimension at each radius
    alpha : array - local power-law index
    """
    log_R = np.log(R)
    log_V = np.log(V)
    alpha = np.gradient(log_V, log_R)
    D = 2.0 * (1.0 - alpha)
    return D, alpha


# ============================================================
# LOAD AND PROCESS ALL SPARC GALAXIES
# ============================================================
def process_sparc(sparc_dir):
    """Load all SPARC rotation curves and compute D(R)."""
    
    files = sorted(glob.glob(os.path.join(sparc_dir, "*_rotmod.dat")))
    
    if len(files) == 0:
        print(f"ERROR: No *_rotmod.dat files found in '{sparc_dir}/'")
        print(f"Download from: https://zenodo.org/records/16284118/files/Rotmod_LTG.zip")
        sys.exit(1)
    
    print(f"Found {len(files)} galaxy files")
    
    results = []
    
    for f in files:
        name = os.path.basename(f).replace("_rotmod.dat", "")
        try:
            # SPARC format: Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul
            data = np.loadtxt(f, comments='#')
            R = data[:, 0]    # kpc
            V = data[:, 1]    # km/s
            errV = data[:, 2] # km/s
            
            # Filter valid points
            mask = (R > 0) & (V > 0)
            R, V, errV = R[mask], V[mask], errV[mask]
            
            if len(R) < MIN_POINTS:
                continue
            
            # Compute D(R)
            D, alpha = compute_D(R, V)
            
            # Summary statistics
            D_outer = np.mean(D[-3:])      # mean of outer 3 points
            V_flat = np.mean(V[-3:])       # flat velocity
            R_max = R[-1]                   # outermost radius
            alpha_outer = np.mean(alpha[-3:])
            
            results.append({
                'name': name,
                'D_outer': D_outer,
                'V_flat': V_flat,
                'R_max': R_max,
                'alpha_outer': alpha_outer,
                'n_points': len(R),
                'R': R, 'V': V, 'errV': errV,
                'D': D, 'alpha': alpha,
            })
            
        except Exception as e:
            print(f"  Warning: could not process {name}: {e}")
    
    print(f"Successfully processed: {len(results)} galaxies\n")
    return results


# ============================================================
# ANALYSIS AND STATISTICS
# ============================================================
def analyze(results):
    """Print summary statistics."""
    
    D_outers = np.array([r['D_outer'] for r in results])
    N = len(results)
    
    n_near_2 = np.sum((D_outers > 1.5) & (D_outers < 2.5))
    
    print("=" * 60)
    print("  RESULTS: EFFECTIVE GRAVITATIONAL DIMENSION")
    print("=" * 60)
    print(f"\n  Total galaxies: {N}")
    print(f"  D_outer in [1.5, 2.5]: {n_near_2} ({100*n_near_2/N:.1f}%)")
    print(f"  Mean D_outer:   {np.mean(D_outers):.3f}")
    print(f"  Median D_outer: {np.median(D_outers):.3f}")
    print(f"  Std D_outer:    {np.std(D_outers):.3f}")
    
    print(f"\n  Distribution:")
    for lo, hi, label in [(-99,1,"D < 1"), (1,1.5,"1-1.5"),
                           (1.5,2.0,"1.5-2.0"), (2.0,2.5,"2.0-2.5"),
                           (2.5,3.0,"2.5-3.0"), (3.0,99,"D > 3")]:
        n = np.sum((D_outers >= lo) & (D_outers < hi))
        bar = "#" * int(n * 40 / N)
        print(f"    {label:>8}: {n:3d} ({100*n/N:4.1f}%) {bar}")
    
    print(f"\n  PREDICTION: D -> 2 at outer radii")
    print(f"  RESULT: Median = {np.median(D_outers):.2f} -- ", end="")
    if 1.5 < np.median(D_outers) < 2.5:
        print("CONFIRMED")
    else:
        print("NEEDS INVESTIGATION")
    print()
    
    return D_outers


# ============================================================
# PLOTTING
# ============================================================
def make_plots(results, D_outers, output_dir):
    """Generate all figures."""
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Figure 1: All D(R) curves + median
    fig, ax = plt.subplots(figsize=(8, 6))
    for r in results:
        Rn = r['R'] / r['R_max']
        ax.plot(Rn, r['D'], 'b-', alpha=0.06, lw=0.8)
    
    # Compute median curve
    bins = np.linspace(0.05, 1.0, 25)
    D_med, D_25, D_75 = [], [], []
    for i in range(len(bins)-1):
        vals = []
        for r in results:
            Rn = r['R'] / r['R_max']
            m = (Rn >= bins[i]) & (Rn < bins[i+1])
            if np.any(m): vals.extend(r['D'][m])
        if vals:
            D_med.append(np.median(vals))
            D_25.append(np.percentile(vals, 25))
            D_75.append(np.percentile(vals, 75))
        else:
            D_med.append(np.nan); D_25.append(np.nan); D_75.append(np.nan)
    
    Rm = 0.5 * (bins[:-1] + bins[1:])
    ax.fill_between(Rm, D_25, D_75, color='red', alpha=0.2, label='25-75 percentile')
    ax.plot(Rm, D_med, 'r-', lw=3, label='Median D(R)')
    ax.axhline(2.0, color='green', ls='--', lw=2, alpha=0.7, label='D = 2 (flat curve)')
    ax.axhline(3.0, color='gray', ls='--', lw=1.5, alpha=0.4, label='D = 3 (Keplerian)')
    ax.set_xlabel('R / R_max', fontsize=14)
    ax.set_ylabel('D(R) = 2(1 - alpha)', fontsize=14)
    ax.set_title(f'Effective Dimension for {len(results)} SPARC Galaxies', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.set_xlim(0, 1); ax.set_ylim(-1, 5)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig1_all_D_curves.png', dpi=150)
    plt.close()
    
    # Figure 2: Histogram
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(D_outers, bins=30, range=(-1, 5), color='steelblue', edgecolor='black')
    ax.axvline(2.0, color='red', lw=3, ls='--', label='D = 2 (prediction)')
    ax.axvline(np.median(D_outers), color='orange', lw=2, label=f'Median = {np.median(D_outers):.2f}')
    ax.set_xlabel('D at outer radii', fontsize=14)
    ax.set_ylabel('Count', fontsize=14)
    ax.set_title('Distribution of D_outer', fontsize=14, fontweight='bold')
    ax.legend(fontsize=12)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{output_dir}/fig2_histogram.png', dpi=150)
    plt.close()
    
    print(f"  Figures saved to {output_dir}/")


# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    
    # Allow custom data directory
    data_dir = sys.argv[1] if len(sys.argv) > 1 else SPARC_DIR
    
    print("\n  No Dark Matter Required:")
    print("  Effective Gravitational Dimension in SPARC Galaxies")
    print("  Nicolae Pascal, 2026\n")
    
    results = process_sparc(data_dir)
    D_outers = analyze(results)
    make_plots(results, D_outers, OUTPUT_DIR)
    
    print("  Done. Check results/ folder for figures.\n")
