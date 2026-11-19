"""Scale-dependent 2D nematic order of points in a binary grid.


The input is a binary 2D grid of 0 and 1 in which each marked cell (1) represents an algae patch center. 
Grid spacings (dx and dy) are supplied in meters. """


import numpy as np
from scipy.spatial import cKDTree


def compute_nematic_order_from_binary(center_2d, dx, dy, r_min, r_max):
    """Compute nematic order from all distinct point pairs in a distance range.

    Parameters
    ----------
    center_2d : 2D NumPy array, shape (ny, nx)
        Binary 2D matrix (1= Algae center; 0= no algae center).
    dx : 2D NumPy array, shape (ny, nx - 1)
        Horizontal distances between adjacent grid positions, in meters.
    dy : 2D NumPy array, shape (ny - 1, nx)
        Vertical distances between adjacent grid positions, in meters.
    r_min, r_max : float 
        Minimum (inclusive) and maximum (inclusive) pair separation, in meters.
        

    Returns
    -------
    S : float
        Nematic order parameter between 0 and 1.
    psi_deg : float
        Preferred axial orientation in degrees in [0, 180).
    C, Q : float
        Mean cosine and sine of twice the pair angle.
    n_pairs : int
        Number of contributing unordered pairs.

    Notes
    -----

    The grid coordinates horizontal positions are accumulated separately along each row, while vertical positions
    use the first column of ``dy``. This assumes that vertical spacings are
    effectively constant across columns, as in the reference dataset.

    The orientation is measured counterclockwise from the positive x-axis.
    When no valid pairs exist, S, psi_deg, C and Q are NaN, and n_pairs is 0.
    For nearly isotropic configurations (S near zero), psi_deg is not reliable.
    """

    # Convert input matrices to NumPy arrays if they are not already
    center_2d = np.asarray(center_2d)
    dx = np.asarray(dx, dtype=float)
    dy = np.asarray(dy, dtype=float)
    #---------------------------------------------------------------


    # Validate input dimensions and values
    if center_2d.ndim != 2:
        raise ValueError("center_2d must be a 2D array.")
    ny, nx = center_2d.shape
    if dx.shape != (ny, nx - 1) or dy.shape != (ny - 1, nx):
        raise ValueError("Expected dx.shape == (ny, nx-1) and dy.shape == (ny-1, nx).")
    if not np.all((center_2d == 0) | (center_2d == 1)):
        raise ValueError("center_2d must contain only 0 and 1.")
    if not (np.all(np.isfinite(dx)) and np.all(np.isfinite(dy)) and np.all(dx > 0) and np.all(dy > 0)):
         raise ValueError('dx and dy must be positive and finite.')
    
    if not np.isfinite(r_min) or r_min < 0:
        raise ValueError("r_min must be finite and nonnegative.")
    if not np.isfinite(r_max) or r_max <= r_min:
        raise ValueError("r_max must be finite and greater than r_min.")

    iy, ix = np.where(center_2d == 1)
    n = len(ix)
    if n < 2:
        raise ValueError('At least two patch centers are required.')

    #---------------------------------------------------------------


    # Calculate spatial coordinates using cumulative distances in meters
    x_grid = np.concatenate((np.zeros((ny, 1)), np.cumsum(dx, axis=1)), axis=1)
    y_axis = np.concatenate(([0.0], np.cumsum(dy[:, 0])))
    x = x_grid[iy, ix]
    y = y_axis[iy]
    coords = np.column_stack((x, y))
    #------------------------------------------------------------------


    # Find all point pairs within r_max
    tree = cKDTree(coords)
    pairs = tree.query_pairs(r=r_max, output_type="ndarray")
    #------------------------------------------------------------------

    # Calculate pairwise displacements and keep only pairs with distances above the minimum
    delta = coords[pairs[:, 1]] - coords[pairs[:, 0]]
    distance = np.hypot(delta[:, 0], delta[:, 1])
    valid = (distance >= r_min) & (distance > 0)
    delta = delta[valid]
    #------------------------------------------------------------------


    if len(delta) == 0:
        return np.nan, np.nan, np.nan, np.nan, 0

    # Calculate the nematic order parameter and preferred orientation
    phi = np.arctan2(delta[:, 1], delta[:, 0])
    C = float(np.mean(np.cos(2 * phi)))
    Q = float(np.mean(np.sin(2 * phi)))
    S = float(np.hypot(C, Q))
    psi_deg = float(np.degrees(0.5 * np.arctan2(Q, C)) % 180.0)
    #------------------------------------------------------------------

    return S, psi_deg, C, Q, len(delta)
