# Nematic Order Parameter (NOP)

This repository provides a Python function to calculate the **nematic order parameter (NOP)** of points in a two-dimensional binary matrix (1 = point, 0 = no point).

The function was originally developed to analyze the spatial distribution of macroalgal patch centers represented in a 2D binary matrix (`center_2d`).

The NOP quantifies the degree of directional alignment of point pairs within a specified distance range. It ranges from 0 (no preferred orientation) to 1 (perfect alignment). The function also calculates the preferred orientation angle of the point pairs. Unlike the conventional NOP based on individual object orientations, this implementation uses the orientations between pairs of points within a specified distance range. Therefore, even a nearly perfectly banded point distribution may yield NOP values well below 1.

## Requirements

- Python 3
- NumPy
- SciPy

## Function

`compute_nematic_order_from_binary(center_2d, dx, dy, r_min, r_max)`

**Inputs:**

- `center_2d`: Two-dimensional binary matrix containing point centers (1 = point, 0 = no point).
- `dx`: Horizontal distances between adjacent grid positions, in meters.
- `dy`: Vertical distances between adjacent grid positions, in meters.
- `r_min`: Minimum pair separation, in meters (inclusive).
- `r_max`: Maximum pair separation, in meters (inclusive).

**Outputs:**

- **`S`**: Nematic order parameter, ranging from 0 to 1.
  \[
  S = \sqrt{C^2 + Q^2}
  \]

- **`psi_deg`**: Preferred axial orientation in degrees, ranging from 0° to less than 180°.
  \[
  \psi = \frac{1}{2}\operatorname{atan2}(Q,C)\pmod{180^\circ}
  \]

- **`C`**: Mean cosine of twice the pair orientation angle.
  \[
  C = \frac{1}{N}\sum_{i=1}^{N}\cos(2\phi_i)
  \]

- **`Q`**: Mean sine of twice the pair orientation angle.
  \[
  Q = \frac{1}{N}\sum_{i=1}^{N}\sin(2\phi_i)
  \]

- **`n_pairs`**: Number of point pairs contributing to the calculation (\(N\)).
## Example usage

```python
import numpy as np
from nematic_order import compute_nematic_order_from_binary

# Example binary matrix
center_2d = np.array([
    [1, 0, 0, 1],
    [0, 0, 1, 0],
    [1, 0, 0, 1],
    [0, 1, 0, 0]
])

# Grid spacing in meters
dx = np.full((4, 3), 100.0)
dy = np.full((3, 4), 100.0)

# Define the distance range in meters
r_min = 0
r_max = 1000

# Calculate nematic order for point pairs separated by 0–1000 m
S, psi_deg, C, Q, n_pairs = compute_nematic_order_from_binary(
    center_2d, dx, dy, r_min, r_max
)

print("Nematic order parameter:", S)
print("Preferred orientation (degrees):", psi_deg)
print("Number of point pairs:", n_pairs)
```

## Notes

The function uses the actual horizontal and vertical grid spacings to calculate point coordinates in meters. Horizontal positions are calculated separately for each row using `dx`, while vertical positions are calculated using the first column of `dy`. This assumes that vertical grid spacing is approximately constant across columns.

The preferred orientation is measured counterclockwise from the positive x-axis. Orientations separated by 180° are equivalent.

No periodic boundary conditions or edge corrections are applied. Boundary effects should therefore be considered when selecting the maximum pair separation relative to the observation domain.

If no valid point pairs exist within the specified distance range, the function returns `NaN` for `S`, `psi_deg`, `C`, and `Q`, and zero for `n_pairs`. When `S` is close to zero, the preferred orientation is poorly defined.

## License

This project is licensed under the MIT License.