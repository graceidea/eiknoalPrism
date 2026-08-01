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
