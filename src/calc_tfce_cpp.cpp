#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <queue>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calc_tfce_cpp(NumericVector t_stat,
                            IntegerVector dims,
                            double E = 0.5,
                            double H = 2.0,
                            double dh = 0.1,
                            int connectivity = 26) {
    int nx = dims[0];
    int ny = dims[1];
    int nz = dims[2];
    int n_voxels = nx * ny * nz;

    NumericVector tfce_map(n_voxels);
    double max_val = 0.0;
    std::vector<bool> is_na(n_voxels, false);

    for (int i = 0; i < n_voxels; ++i) {
        if (NumericVector::is_na(t_stat[i]) || std::isnan(t_stat[i])) {
            is_na[i] = true;
            tfce_map[i] = NA_REAL;
        } else {
            tfce_map[i] = 0.0;
            if (t_stat[i] > max_val) {
                max_val = t_stat[i];
            }
        }
    }

    if (max_val <= 0.0) {
        tfce_map.attr("dim") = dims;
        return tfce_map;
    }

    std::vector<int> dx, dy, dz;
    for (int x = -1; x <= 1; ++x) {
        for (int y = -1; y <= 1; ++y) {
            for (int z = -1; z <= 1; ++z) {
                if (x == 0 && y == 0 && z == 0) continue;
                int dist_sq = x * x + y * y + z * z;
                if ((connectivity == 6 && dist_sq == 1) ||
                    (connectivity == 18 && dist_sq <= 2) ||
                    (connectivity == 26)) {
                    dx.push_back(x);
                    dy.push_back(y);
                    dz.push_back(z);
                }
            }
        }
    }
    int n_neighbors = dx.size();

    std::vector<int> visited(n_voxels, 0);
    int current_visited_flag = 0;

    std::vector<int> cluster_indices;
    cluster_indices.reserve(10000);
    std::queue<int> q;

    for (double h = dh; h <= max_val + 1e-7; h += dh) {
        current_visited_flag++;
        double h_factor = std::pow(h, H) * dh;

        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                for (int x = 0; x < nx; ++x) {
                    int idx = x + nx * (y + ny * z);

                    if (is_na[idx] || t_stat[idx] < h || visited[idx] == current_visited_flag) {
                        continue;
                    }

                    cluster_indices.clear();
                    q.push(idx);
                    visited[idx] = current_visited_flag;

                    while (!q.empty()) {
                        int curr = q.front();
                        q.pop();
                        cluster_indices.push_back(curr);

                        int cz = curr / (nx * ny);
                        int rem = curr % (nx * ny);
                        int cy = rem / nx;
                        int cx = rem % nx;

                        for (int k = 0; k < n_neighbors; ++k) {
                            int nx_pos = cx + dx[k];
                            int ny_pos = cy + dy[k];
                            int nz_pos = cz + dz[k];

                            if (nx_pos >= 0 && nx_pos < nx &&
                                ny_pos >= 0 && ny_pos < ny &&
                                nz_pos >= 0 && nz_pos < nz) {

                                int nbr_idx = nx_pos + nx * (ny_pos + ny * nz_pos);

                                if (!is_na[nbr_idx] && t_stat[nbr_idx] >= h && visited[nbr_idx] != current_visited_flag) {
                                    visited[nbr_idx] = current_visited_flag;
                                    q.push(nbr_idx);
                                }
                            }
                        }
                    }

                    double extent = static_cast<double>(cluster_indices.size());
                    double tfce_inc = std::pow(extent, E) * h_factor;

                    for (int member_idx : cluster_indices) {
                        tfce_map[member_idx] += tfce_inc;
                    }
                }
            }
        }
    }

    tfce_map.attr("dim") = dims;
    return tfce_map;
}