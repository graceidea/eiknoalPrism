#!/usr/bin/env python3
"""
Visualization script for HDF5 results from the Eikonal solver
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

def visualize_travel_times(filename="complete_results.h5", analytical_filename="analytical_solution.h5", error_filename="error_analysis.h5"):
    """
    Visualize travel time results from HDF5 file with analytical comparison
    """
    try:
        with h5py.File(filename, 'r') as f:
            # Read data
            xx = f['xx'][:]
            yy = f['yy'][:]
            T = f['T'][:]
            
            # Optional: read other data if available
            if 'a' in f:
                a = f['a'][:]
            if 'b' in f:
                b = f['b'][:]
            if 'c' in f:
                c = f['c'][:]
            if 'slowness' in f:
                slowness = f['slowness'][:]
        
        # Try to read analytical solution
        T_analytical = None
        try:
            with h5py.File(analytical_filename, 'r') as f_analytical:
                T_analytical = f_analytical['T_analytical'][:]
        except FileNotFoundError:
            print(f"Analytical solution file {analytical_filename} not found. Skipping comparison.")
        
        # Try to read error analysis
        error_data = {}
        try:
            with h5py.File(error_filename, 'r') as f_error:
                error_data['error'] = f_error['error'][:]
                error_data['abs_error'] = f_error['abs_error'][:]
                error_data['rel_error'] = f_error['rel_error'][:]
        except FileNotFoundError:
            print(f"Error analysis file {error_filename} not found. Skipping error plots.")
        
        print(f"Loaded data from {filename}")
        print(f"Grid size: {len(xx)} x {len(yy)}")
        print(f"X range: [{xx.min():.2f}, {xx.max():.2f}]")
        print(f"Y range: [{yy.min():.2f}, {yy.max():.2f}]")
        print(f"Travel time range: [{T.min():.4f}, {T.max():.4f}]")
        
        if T_analytical is not None:
            print(f"Analytical solution range: [{T_analytical.min():.4f}, {T_analytical.max():.4f}]")
            print(f"Max absolute error: {np.abs(T - T_analytical).max():.6f}")
            print(f"Mean absolute error: {np.abs(T - T_analytical).mean():.6f}")
        
        # Create meshgrid for plotting
        X, Y = np.meshgrid(xx, yy)
        
        # Determine number of subplots
        n_plots = 6 if T_analytical is None else 9
        fig = plt.figure(figsize=(20, 15) if T_analytical is not None else (15, 10))
        
        # Plot 1: Numerical solution contours
        ax1 = plt.subplot(3, 3, 1) if T_analytical is not None else plt.subplot(2, 3, 1)
        contour = ax1.contour(X, Y, T.T, levels=20, colors='black', alpha=0.6, linewidths=0.5)
        contourf = ax1.contourf(X, Y, T.T, levels=50, cmap='viridis')
        plt.colorbar(contourf, ax=ax1, label='Travel Time')
        ax1.set_title('Numerical Solution')
        ax1.set_xlabel('X')
        ax1.set_ylabel('Y')
        ax1.set_aspect('equal')
        
        # Plot 2: 3D surface of numerical solution
        ax2 = fig.add_subplot(3, 3, 2, projection='3d') if T_analytical is not None else fig.add_subplot(2, 3, 2, projection='3d')
        surf = ax2.plot_surface(X, Y, T.T, cmap='viridis', alpha=0.8)
        ax2.set_title('Numerical Solution (3D)')
        ax2.set_xlabel('X')
        ax2.set_ylabel('Y')
        ax2.set_zlabel('Travel Time')
        
        # Plot 3: Travel time along center lines
        ax3 = plt.subplot(3, 3, 3) if T_analytical is not None else plt.subplot(2, 3, 3)
        center_idx = len(yy) // 2
        ax3.plot(xx, T[:, center_idx], 'b-', linewidth=2, label='Numerical (Y=0)')
        center_idx = len(xx) // 2
        ax3.plot(yy, T[center_idx, :], 'r-', linewidth=2, label='Numerical (X=0)')
        
        if T_analytical is not None:
            center_idx = len(yy) // 2
            ax3.plot(xx, T_analytical[:, center_idx], 'b--', linewidth=2, label='Analytical (Y=0)')
            center_idx = len(xx) // 2
            ax3.plot(yy, T_analytical[center_idx, :], 'r--', linewidth=2, label='Analytical (X=0)')
        
        ax3.set_title('Travel Time Profiles')
        ax3.set_xlabel('Distance')
        ax3.set_ylabel('Travel Time')
        ax3.legend()
        ax3.grid(True)
        
        # If analytical solution is available, add comparison plots
        if T_analytical is not None:
            # Plot 4: Analytical solution contours
            ax4 = plt.subplot(3, 3, 4)
            contour_analytical = ax4.contour(X, Y, T_analytical.T, levels=20, colors='black', alpha=0.6, linewidths=0.5)
            contourf_analytical = ax4.contourf(X, Y, T_analytical.T, levels=50, cmap='viridis')
            plt.colorbar(contourf_analytical, ax=ax4, label='Travel Time')
            ax4.set_title('Analytical Solution')
            ax4.set_xlabel('X')
            ax4.set_ylabel('Y')
            ax4.set_aspect('equal')
            
            # Plot 5: Error contours
            if 'abs_error' in error_data:
                ax5 = plt.subplot(3, 3, 5)
                error_contourf = ax5.contourf(X, Y, error_data['abs_error'].T, levels=50, cmap='Reds')
                plt.colorbar(error_contourf, ax=ax5, label='Absolute Error')
                ax5.set_title('Absolute Error')
                ax5.set_xlabel('X')
                ax5.set_ylabel('Y')
                ax5.set_aspect('equal')
            
            # Plot 6: Relative error
            if 'rel_error' in error_data:
                ax6 = plt.subplot(3, 3, 6)
                rel_error_percent = error_data['rel_error'] * 100
                rel_contourf = ax6.contourf(X, Y, rel_error_percent.T, levels=50, cmap='plasma')
                plt.colorbar(rel_contourf, ax=ax6, label='Relative Error (%)')
                ax6.set_title('Relative Error')
                ax6.set_xlabel('X')
                ax6.set_ylabel('Y')
                ax6.set_aspect('equal')
            
            # Plot 7: Error along center lines
            ax7 = plt.subplot(3, 3, 7)
            if 'abs_error' in error_data:
                center_idx = len(yy) // 2
                ax7.plot(xx, error_data['abs_error'][:, center_idx], 'b-', linewidth=2, label='Error along Y=0')
                center_idx = len(xx) // 2
                ax7.plot(yy, error_data['abs_error'][center_idx, :], 'r-', linewidth=2, label='Error along X=0')
                ax7.set_title('Error Profiles')
                ax7.set_xlabel('Distance')
                ax7.set_ylabel('Absolute Error')
                ax7.legend()
                ax7.grid(True)
            
            # Plot 8: Statistics comparison
            ax8 = plt.subplot(3, 3, 8)
            stats_data = [
                ['Min', T.min(), T_analytical.min(), np.abs(T.min() - T_analytical.min())],
                ['Max', T.max(), T_analytical.max(), np.abs(T.max() - T_analytical.max())],
                ['Mean', T.mean(), T_analytical.mean(), np.abs(T.mean() - T_analytical.mean())],
                ['Std', T.std(), T_analytical.std(), np.abs(T.std() - T_analytical.std())]
            ]
            
            table_text = []
            for row in stats_data:
                table_text.append([f"{row[0]}", f"{row[1]:.4f}", f"{row[2]:.4f}", f"{row[3]:.6f}"])
            
            table = ax8.table(cellText=table_text,
                             colLabels=['Statistic', 'Numerical', 'Analytical', 'Abs Diff'],
                             cellLoc='center',
                             loc='center')
            table.auto_set_font_size(False)
            table.set_fontsize(9)
            table.scale(1.2, 1.5)
            ax8.set_title('Statistics Comparison')
            ax8.axis('off')
            
            # Plot 9: Error statistics
            ax9 = plt.subplot(3, 3, 9)
            if 'abs_error' in error_data and 'rel_error' in error_data:
                error_stats_text = [
                    f"Max abs error: {error_data['abs_error'].max():.6f}",
                    f"Mean abs error: {error_data['abs_error'].mean():.6f}",
                    f"RMS error: {np.sqrt(np.mean(error_data['error']**2)):.6f}",
                    f"Max rel error: {error_data['rel_error'].max()*100:.4f}%",
                    f"Mean rel error: {error_data['rel_error'].mean()*100:.4f}%"
                ]
                
                for i, text in enumerate(error_stats_text):
                    ax9.text(0.1, 0.8 - i*0.15, text, transform=ax9.transAxes, fontsize=12)
                
                ax9.set_title('Error Statistics')
                ax9.set_xlim(0, 1)
                ax9.set_ylim(0, 1)
                ax9.axis('off')
        
        else:
            # Original plots for when no analytical solution is available
            ax4 = plt.subplot(2, 3, 4)
            dy, dx = np.gradient(T.T)
            grad_mag = np.sqrt(dx**2 + dy**2)
            im = ax4.imshow(grad_mag, extent=[xx.min(), xx.max(), yy.min(), yy.max()], 
                           origin='lower', cmap='plasma')
            plt.colorbar(im, ax=ax4, label='Gradient Magnitude')
            ax4.set_title('Travel Time Gradient')
            ax4.set_xlabel('X')
            ax4.set_ylabel('Y')
            
            # Plot 5: Slowness field (if available)
            if 'slowness' in locals():
                ax5 = plt.subplot(2, 3, 5)
                im = ax5.imshow(slowness.T, extent=[xx.min(), xx.max(), yy.min(), yy.max()], 
                               origin='lower', cmap='coolwarm')
                plt.colorbar(im, ax=ax5, label='Slowness')
                ax5.set_title('Slowness Field')
                ax5.set_xlabel('X')
                ax5.set_ylabel('Y')
            
            # Plot 6: Statistics
            ax6 = plt.subplot(2, 3, 6)
            ax6.text(0.1, 0.8, f"Grid: {len(xx)} × {len(yy)}", transform=ax6.transAxes, fontsize=12)
            ax6.text(0.1, 0.7, f"Min T: {T.min():.4f}", transform=ax6.transAxes, fontsize=12)
            ax6.text(0.1, 0.6, f"Max T: {T.max():.4f}", transform=ax6.transAxes, fontsize=12)
            ax6.text(0.1, 0.5, f"Mean T: {T.mean():.4f}", transform=ax6.transAxes, fontsize=12)
            ax6.text(0.1, 0.4, f"Std T: {T.std():.4f}", transform=ax6.transAxes, fontsize=12)
            ax6.set_title('Statistics')
            ax6.set_xlim(0, 1)
            ax6.set_ylim(0, 1)
            ax6.axis('off')
        
        plt.tight_layout()
        plt.savefig('travel_time_visualization.png', dpi=300, bbox_inches='tight')
        plt.show()
        
    except FileNotFoundError:
        print(f"Error: File {filename} not found. Please run the C++ solver first.")
    except Exception as e:
        print(f"Error reading HDF5 file: {e}")

def compare_h5_files(file1="travel_times.h5", file2="complete_results.h5"):
    """
    Compare travel time data from two HDF5 files
    """
    try:
        with h5py.File(file1, 'r') as f1, h5py.File(file2, 'r') as f2:
            T1 = f1['T'][:]
            T2 = f2['T'][:]
            
            diff = np.abs(T1 - T2)
            print(f"Maximum difference: {diff.max():.2e}")
            print(f"Mean difference: {diff.mean():.2e}")
            print(f"Files are {'identical' if diff.max() < 1e-12 else 'different'}")
            
    except Exception as e:
        print(f"Error comparing files: {e}")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "complete_results.h5"
    
    print("Visualizing Eikonal solver results...")
    visualize_travel_times(filename)
    
    # Optional: compare files if both exist
    try:
        compare_h5_files()
    except:
        pass
