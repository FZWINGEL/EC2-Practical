#!/usr/bin/env python3
"""
Smoothing comparison script using actual dQdE data.

This script loads your actual coin cell data and demonstrates different
smoothing methods on the real dQdE curves.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plot import (
    SmoothingConfig, DEFAULT_SMOOTHING, smooth, 
    read_coin_csv, identify_cc_segments, integrate_capacity_per_segment, 
    estimate_cycles_from_segments, MASS_MG
)

def load_and_process_data(csv_path: Path, sample_id: str):
    """Load and process coin cell data."""
    print(f"Loading data from {csv_path.name}...")
    
    # Load and process data
    df = read_coin_csv(csv_path)
    df = identify_cc_segments(df)
    df = integrate_capacity_per_segment(df)
    df = estimate_cycles_from_segments(df)
    
    # Get mass
    mass_mg = MASS_MG.get(sample_id, 9.70)  # Default fallback
    mass_g = mass_mg / 1000.0
    
    # Normalize capacity
    df["cap_mAh_g_seg"] = df["cap_mAh_seg"] / mass_g
    
    return df, mass_g

def extract_dqdE_data(df: pd.DataFrame):
    """Extract dQdE data from processed DataFrame."""
    dqdE_data = []
    
    for cyc, gcyc in df.groupby("cycle_est", sort=True):
        if cyc == 0:
            continue
        for sid, gseg in gcyc.groupby("seg_id", sort=True):
            if not bool(gseg["is_discharge"].iloc[0]):
                continue
            
            E = gseg["voltage_v"].values
            Q = gseg["cap_mAh_g_seg"].abs().values
            
            if len(E) < 5:
                continue
                
            dqdE_data.append({
                'cycle': cyc,
                'segment': sid,
                'E': E,
                'Q': Q,
                'label': f'Cycle {cyc}'
            })
    
    return dqdE_data

def compute_dqdE_with_smoothing(E: np.ndarray, Q: np.ndarray, smoothing_config: SmoothingConfig):
    """Compute dQdE with specified smoothing."""
    # Smooth E and Q
    Es = smooth(E, smoothing_config)
    Qs = smooth(Q, smoothing_config)
    
    # Compute derivative dQ/dE
    with np.errstate(divide='ignore', invalid='ignore'):
        dQ_dE = np.gradient(Qs, Es)
    
    # Apply guards against tiny dE and spikes
    dE_local = np.gradient(Es)
    tiny = max(1e-6, np.nanmedian(np.abs(dE_local)) * 0.01)
    dQ_dE = np.where(np.abs(dE_local) < tiny, np.nan, dQ_dE)
    
    # Clip extreme outliers
    finite_mask = np.isfinite(dQ_dE)
    if finite_mask.sum() > 20:
        lo, hi = np.nanpercentile(dQ_dE[finite_mask], [1.0, 99.0])
        dQ_dE = np.clip(dQ_dE, lo, hi)
    
    # Lightly smooth the derivative
    derivative_config = SmoothingConfig(
        method=smoothing_config.method,
        window_length=min(101, max(7, (len(dQ_dE)//50)*2 + 1)),
        polyorder=smoothing_config.polyorder
    )
    dQ_dE = smooth(dQ_dE, derivative_config)
    
    return Es, dQ_dE

def create_smoothing_comparison(dqdE_data: list, output_path: Path):
    """Create comparison plot of different smoothing methods."""
    
    # Different smoothing configurations
    configs = {
        'No Smoothing': SmoothingConfig(method=SmoothingConfig.NONE),
        'Savitzky-Golay (light)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=15, polyorder=2),
        'Savitzky-Golay (default)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=21, polyorder=3),
        'Savitzky-Golay (heavy)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=31, polyorder=4),
        'Moving Average': SmoothingConfig(method=SmoothingConfig.MOVING_AVERAGE, window_length=15),
    }
    
    # Create subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    axes = axes.flatten()
    
    colors = plt.cm.tab10.colors
    
    # Plot each smoothing method
    for i, (name, config) in enumerate(configs.items()):
        ax = axes[i]
        
        # Plot first few cycles for clarity
        cycles_to_plot = min(3, len(dqdE_data))
        
        for j, data in enumerate(dqdE_data[:cycles_to_plot]):
            try:
                Es, dQ_dE = compute_dqdE_with_smoothing(data['E'], data['Q'], config)
                
                ax.plot(Es, dQ_dE, 
                       color=colors[j % len(colors)], 
                       label=data['label'],
                       linewidth=1.5, alpha=0.8)
            except Exception as e:
                print(f"Warning: Failed to process {data['label']} with {name}: {e}")
                continue
        
        ax.set_xlabel("E [V]", fontsize=12)
        ax.set_ylabel("dQ/dE [mAh/(g·V)]", fontsize=12)
        ax.set_title(f"{name}", fontsize=14, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        if i == 0:  # Only show legend on first plot
            ax.legend(fontsize=10)
    
    # Remove the last empty subplot if odd number
    if len(configs) < len(axes):
        fig.delaxes(axes[-1])
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Smoothing comparison saved to: {output_path}")

def create_overlay_comparison(dqdE_data: list, output_path: Path):
    """Create overlay comparison showing all methods on same plot."""
    
    configs = {
        'No Smoothing': SmoothingConfig(method=SmoothingConfig.NONE),
        'Savitzky-Golay (light)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=15, polyorder=2),
        'Savitzky-Golay (default)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=21, polyorder=3),
        'Savitzky-Golay (heavy)': SmoothingConfig(method=SmoothingConfig.SAVITZKY_GOLAY, window_length=31, polyorder=4),
        'Moving Average': SmoothingConfig(method=SmoothingConfig.MOVING_AVERAGE, window_length=15),
    }
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    colors = ['red', 'blue', 'green', 'orange', 'purple']
    linestyles = ['-', '--', '-.', ':', '-']
    
    # Use first cycle for overlay comparison
    if dqdE_data:
        data = dqdE_data[0]
        
        for i, (name, config) in enumerate(configs.items()):
            try:
                Es, dQ_dE = compute_dqdE_with_smoothing(data['E'], data['Q'], config)
                
                ax.plot(Es, dQ_dE, 
                       color=colors[i], 
                       linestyle=linestyles[i],
                       label=name,
                       linewidth=2, alpha=0.8)
            except Exception as e:
                print(f"Warning: Failed to process with {name}: {e}")
                continue
    
    ax.set_xlabel("E [V]", fontsize=14)
    ax.set_ylabel("dQ/dE [mAh/(g·V)]", fontsize=14)
    ax.set_title(f"dQdE Smoothing Comparison - {data['label']}", fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    print(f"Overlay comparison saved to: {output_path}")

def main():
    """Main function to run smoothing comparison."""
    
    # Find CSV files
    csv_files = list(Path('.').glob('*.csv'))
    if not csv_files:
        print("No CSV files found in current directory!")
        return 1
    
    print(f"Found {len(csv_files)} CSV files:")
    for csv_file in csv_files:
        print(f"  - {csv_file.name}")
    
    # Process each CSV file
    for csv_path in csv_files:
        # Extract sample ID from filename
        sample_id = csv_path.stem.split('_')[-1]  # Get last part after underscore
        
        try:
            print(f"\n{'='*60}")
            print(f"Processing: {csv_path.name}")
            print(f"Sample ID: {sample_id}")
            
            # Load and process data
            df, mass_g = load_and_process_data(csv_path, sample_id)
            
            # Extract dQdE data
            dqdE_data = extract_dqdE_data(df)
            
            if not dqdE_data:
                print("No discharge segments found!")
                continue
            
            print(f"Found {len(dqdE_data)} discharge segments")
            
            # Create comparison plots
            comparison_path = Path(f"../figures/T3.1_dQdE_smoothing_comparison_{sample_id}.png")
            overlay_path = Path(f"../figures/T3.1_dQdE_smoothing_overlay_{sample_id}.png")
            
            create_smoothing_comparison(dqdE_data, comparison_path)
            create_overlay_comparison(dqdE_data, overlay_path)
            
        except Exception as e:
            print(f"Error processing {csv_path.name}: {e}")
            continue
    
    print(f"\n{'='*60}")
    print("Smoothing comparison complete!")
    print("\nTo change smoothing method, edit DEFAULT_SMOOTHING in plot.py:")
    print("  - SmoothingConfig.NONE: No smoothing")
    print("  - SmoothingConfig.SAVITZKY_GOLAY: Savitzky-Golay filter")
    print("  - SmoothingConfig.MOVING_AVERAGE: Moving average")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
