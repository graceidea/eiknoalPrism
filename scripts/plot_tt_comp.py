# %%
import matplotlib.pyplot as plt
import h5py
import numpy as np

# %%
with h5py.File('build/travel_time_gpu_result.h5', 'r') as f:
    x = f['xx'][:]
    y = f['yy'][:]
    tgpu = f['T'][:]

with h5py.File('build/travel_time_cpu_result.h5', 'r') as f:
    tcpu = f['T'][:]

# %%
plt.figure(figsize=(12, 5))
plt.subplot(1, 3, 1)
plt.title('GPU Travel Time')
plt.imshow(tgpu.T, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower')
plt.colorbar(label='Time (s)')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.subplot(1, 3, 2)
plt.title('CPU Travel Time')
plt.imshow(tcpu.T, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower')
plt.colorbar(label='Time (s)')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.subplot(1, 3, 3)
plt.title('Difference (GPU - CPU)')
diff = tgpu - tcpu
vmax = np.max(np.abs(diff))
plt.imshow(diff.T, extent=(x.min(), x.max(), y.min(), y.max()), origin='lower', cmap='bwr', vmin=-vmax, vmax=vmax)
plt.colorbar(label='Time Difference (s)')
plt.xlabel('X (m)')
plt.ylabel('Y (m)')
plt.tight_layout()
plt.savefig('travel_time_comparison.png', dpi=300)


