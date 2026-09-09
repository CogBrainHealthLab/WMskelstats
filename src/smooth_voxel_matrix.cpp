#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;

// Helper function to build 1D Gaussian kernel
std::vector<double> make_gaussian_kernel(double sigma) {
  if (sigma <= 1e-6) {
    return {1.0};
  }
  int radius = static_cast<int>(std::ceil(3.0 * sigma));
  int size = 2 * radius + 1;
  std::vector<double> kernel(size);
  double sum = 0.0;

  for (int i = -radius; i <= radius; ++i) {
    double val = std::exp(-0.5 * (i * i) / (sigma * sigma));
    kernel[i + radius] = val;
    sum += val;
  }
  for (int i = 0; i < size; ++i) {
    kernel[i] /= sum;
  }
  return kernel;
}

// Fast 3D separable Gaussian convolution on flat 1D array
void smooth_3d_separable(
    const std::vector<double>& input,
    std::vector<double>& output,
    int dimX, int dimY, int dimZ,
    const std::vector<double>& kX,
    const std::vector<double>& kY,
    const std::vector<double>& kZ) 
{
  int total_size = dimX * dimY * dimZ;
  std::vector<double> tmp1(total_size, 0.0);
  std::vector<double> tmp2(total_size, 0.0);

  int radX = (kX.size() - 1) / 2;
  int radY = (kY.size() - 1) / 2;
  int radZ = (kZ.size() - 1) / 2;

  // Pass 1: X direction
  for (int z = 0; z < dimZ; ++z) {
    for (int y = 0; y < dimY; ++y) {
      int base_yz = dimX * (y + dimY * z);
      for (int x = 0; x < dimX; ++x) {
        double val = 0.0;
        for (int k = -radX; k <= radX; ++k) {
          int nx = x + k;
          if (nx >= 0 && nx < dimX) {
            val += input[base_yz + nx] * kX[k + radX];
          }
        }
        tmp1[base_yz + x] = val;
      }
    }
  }

  // Pass 2: Y direction
  for (int z = 0; z < dimZ; ++z) {
    for (int x = 0; x < dimX; ++x) {
      for (int y = 0; y < dimY; ++y) {
        double val = 0.0;
        for (int k = -radY; k <= radY; ++k) {
          int ny = y + k;
          if (ny >= 0 && ny < dimY) {
            val += tmp1[x + dimX * (ny + dimY * z)] * kY[k + radY];
          }
        }
        tmp2[x + dimX * (y + dimY * z)] = val;
      }
    }
  }

  // Pass 3: Z direction
  for (int y = 0; y < dimY; ++y) {
    for (int x = 0; x < dimX; ++x) {
      for (int z = 0; z < dimZ; ++z) {
        double val = 0.0;
        for (int k = -radZ; k <= radZ; ++k) {
          int nz = z + k;
          if (nz >= 0 && nz < dimZ) {
            val += tmp2[x + dimX * (y + dimY * nz)] * kZ[k + radZ];
          }
        }
        output[x + dimX * (y + dimY * z)] = val;
      }
    }
  }
}

// [[Rcpp::export]]
NumericMatrix cpp_smooth_voxel_matrix(
    NumericMatrix data_mat,
    IntegerMatrix coords,
    NumericVector sigma_voxels) 
{
  int n_subjects = data_mat.nrow();
  int n_voxels = data_mat.ncol();

  // Detect index offset (1-based from R vs 0-based)
  int min_val = coords(0, 0);
  for (int i = 0; i < n_voxels; ++i) {
    for (int j = 0; j < 3; ++j) {
      if (coords(i, j) < min_val) min_val = coords(i, j);
    }
  }
  int offset = (min_val == 1) ? 1 : 0;

  // Determine maximum dimensions for 3D bounding box
  int dimX = 0, dimY = 0, dimZ = 0;
  for (int i = 0; i < n_voxels; ++i) {
    int x = coords(i, 0) - offset;
    int y = coords(i, 1) - offset;
    int z = coords(i, 2) - offset;

    if (x + 1 > dimX) dimX = x + 1;
    if (y + 1 > dimY) dimY = y + 1;
    if (z + 1 > dimZ) dimZ = z + 1;
  }

  // Pre-calculate 1D flat indices for linear memory access
  std::vector<int> flat_indices(n_voxels);
  for (int i = 0; i < n_voxels; ++i) {
    int x = coords(i, 0) - offset;
    int y = coords(i, 1) - offset;
    int z = coords(i, 2) - offset;
    flat_indices[i] = x + dimX * (y + dimY * z);
  }

  int total_grid_size = dimX * dimY * dimZ;

  // Build 1D Gaussian kernels per axis
  std::vector<double> kX = make_gaussian_kernel(sigma_voxels[0]);
  std::vector<double> kY = make_gaussian_kernel(sigma_voxels[1]);
  std::vector<double> kZ = make_gaussian_kernel(sigma_voxels[2]);

  // Precompute normalisation weight grid once
  std::vector<double> mask(total_grid_size, 0.0);
  for (int i = 0; i < n_voxels; ++i) {
    mask[flat_indices[i]] = 1.0;
  }

  std::vector<double> weight_grid(total_grid_size, 0.0);
  smooth_3d_separable(mask, weight_grid, dimX, dimY, dimZ, kX, kY, kZ);

  NumericMatrix out_mat(n_subjects, n_voxels);

  // Sequential iteration over subjects
  for (int s = 0; s < n_subjects; ++s) {
    std::vector<double> vol_in(total_grid_size, 0.0);
    std::vector<double> vol_out(total_grid_size, 0.0);

    for (int v = 0; v < n_voxels; ++v) {
      vol_in[flat_indices[v]] = data_mat(s, v);
    }

    smooth_3d_separable(vol_in, vol_out, dimX, dimY, dimZ, kX, kY, kZ);

    for (int v = 0; v < n_voxels; ++v) {
      int idx = flat_indices[v];
      double w = weight_grid[idx];
      if (w > 1e-12) {
        out_mat(s, v) = vol_out[idx] / w;
      } else {
        out_mat(s, v) = NA_REAL;
      }
    }
  }

  return out_mat;
}