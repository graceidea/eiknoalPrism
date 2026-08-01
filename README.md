### ⚙️ Core Pipeline (What It Does)

The `eiknoalPrism` engine executes a 5-step geometric and differential pipeline to generate prism meshes from a closed polygon boundary using the Eikonal equation.

| Step | Operation Phase | Input Asset | Output Generated File | Output Visualization |
| :---: | :--- | :--- | :--- | :---: |
| **1** | **Source Function (SF) Processing** | `input.obj` | `polygon.msh` | [View Geometry](#1-source-function) |
| **2** | **Background Cartesian Grid (BCG)** | `polygon.msh` | `grid.msh` | [View Grid](#2-background-cartesian-grid) |
| **3** | **Eikonal Distance Field Computation** | `grid.msh` | `distance_exact.vtk` | [View Field](#3-distance-field) |
| **4** | **Gradient of Distance Function (GDF)**| `distance_exact.vtk` | `gradient_magnitude.vtk`| [View Gradient](#4-gradient-field) |
| **5** | **Prism Mesh Generation** | `gradient_magnitude.vtk`| `prism.msh` | [View Mesh](#5-prism-mesh) |

---

#### 1. Source Function (SF)
Reads the initial closed boundary polygon setup.
<img src="https://github.com/user-attachments/assets/60c52924-c825-404e-9e0e-464f4000b77d" width="450" alt="Source Function Polygon" />

#### 2. Background Cartesian Grid (BCG)
Generates the background bounding grid tracking the boundary path.
<img src="https://github.com/user-attachments/assets/3de02070-47e9-49e8-b4d7-212434b02843" width="450" alt="Background Cartesian Grid" />

#### 3. Distance Field
Solves the Eikonal equation to yield a smooth, continuous scalar distance field map.
<img src="https://github.com/user-attachments/assets/ba4d7b17-d050-4af9-8cf4-bfb7cf9b7eb0" width="450" alt="Eikonal Distance Field" />

#### 4. Gradient Field (GDF)
Computes the vector gradients and magnitude spikes along the medial axis partitions.
<img src="https://github.com/user-attachments/assets/bbf6c3c2-ff62-4d83-ac45-c895d6ec39a6" width="450" alt="Gradient of Distance Function" />

#### 5. Prism Mesh
Constructs structured layers of anisotropic prism elements radiating away from the geometry boundary.
<img src="https://github.com/user-attachments/assets/2c1f973c-178f-4e40-9a9a-59cc7c8848ab" width="450" alt="Prism Mesh Generation" />






















# eiknoalPrism

A modern project implementation for processing Source Functions.

---

### 🔧 Installation

```bash
# Clone the repository
git clone git@github.com:graceidea/eiknoalPrism.git
cd eiknoalPrism/

# Build the binaries
mkdir build && cd build
cmake ..
make
```

### 🧪 Running Tests

```bash
cd eiknoalPrism/build/bin/
./test
```

### ⚙️ Core Functionality

The application processes **Source Functions (SF)** represented as closed polygons.

| Phase | Target File Path | Format Description |
| :--- | :--- | :--- |
| **Input** | `eiknoalPrism/tests/example/input.obj` | Source closed polygon file |
| **Output** | `eiknoalPrism/tests/example/polygon.msh` | Generated mesh function file |

> ℹ️ **Note:** The example dataset can be verified locally inside the `tests/example` directory immediately after running `./test`.



=============WHAT IT DOES
The code can do the following:

1, Read in the Source Function (SF), which is a closed polygon in this example. The file is in eiknoalPrism/tests/example folder, named as input.obj, and saved a new file named as polygon.msh in the eiknoalPrism/tests/example folder.

<img width="573" height="567" alt="image" src="https://github.com/user-attachments/assets/60c52924-c825-404e-9e0e-464f4000b77d" />

2, Generate a Background Cartesian Grid (BCG, the file will be saved in eiknoalPrism/tests/example folder, named as grid.msh after running ./test) based on input Source Function (SF). The following picture shows the BCG & SF.

<img width="547" height="550" alt="image" src="https://github.com/user-attachments/assets/3de02070-47e9-49e8-b4d7-212434b02843" />

3, Compute the distance function using Eikonal Equation (the file will be saved in eiknoalPrism/tests/example folder, named as distance_exact.vtk after running ./test), please see the attached distance field.

<img width="530" height="532" alt="image" src="https://github.com/user-attachments/assets/ba4d7b17-d050-4af9-8cf4-bfb7cf9b7eb0" />

4, Compute Gradient of Distance Function (GDF, the file will be saved in eiknoalPrism/tests/example folder, named as  gradient_magnitude.vtk after running ./test)), please see the attached example of GDF.

<img width="529" height="533" alt="image" src="https://github.com/user-attachments/assets/bbf6c3c2-ff62-4d83-ac45-c895d6ec39a6" />

5, Generate Prism mesh for SF (the file will be saved in eiknoalPrism/tests/example folder, named as prism.msh after running ./test),.

<img width="505" height="534" alt="image" src="https://github.com/user-attachments/assets/2c1f973c-178f-4e40-9a9a-59cc7c8848ab" />
