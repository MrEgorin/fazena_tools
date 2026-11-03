/*
 * FAZENA Tools — ComTarget 3D C3DD Molecular Shape Descriptor Engine
 * High-performance C implementation of 3D Hu Moments, Surface, Volume, and Shape Invariants.
 * Reference: Li et al. ComTarget: Small-Molecule Target Prediction with Combinatorial Modeling.
 *            Pharmaceuticals 2026, 19, 715.
 * Author: Yehor Franchuk (@MrEgorin), FAZENA (fazena.org)
 * Affiliation: Institute of Cell Biology NASU
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * S. Alvarez vdW Radii (Dalton Trans. 2013, 42, 8617-8636)
 * Indexed directly by atomic number (Z = 0..92)
 */
static const float ALVAREZ_VDW[100] = {
    0.0f,  1.20f, 1.43f, 2.12f, 1.98f, 1.91f, 1.77f, 1.66f, 1.50f, 1.46f, 1.58f, /* 0-10 */
    2.50f, 2.51f, 2.25f, 2.19f, 1.90f, 1.89f, 1.82f, 1.83f, 2.73f, 2.62f, /* 11-20 */
    2.58f, 2.46f, 2.42f, 2.45f, 2.45f, 2.44f, 2.40f, 2.40f, 2.38f, 2.39f, /* 21-30 */
    2.32f, 2.29f, 1.88f, 1.82f, 1.86f, 2.25f, 2.78f, 2.52f, 2.49f, 2.40f, /* 31-40 */
    2.33f, 2.30f, 2.35f, 2.34f, 2.37f, 2.38f, 2.44f, 2.46f, 2.42f, 2.42f, /* 41-50 */
    2.47f, 1.99f, 2.04f, 2.06f, 2.82f, 2.68f, 2.58f, 2.53f, 2.49f, 2.48f, /* 51-60 */
    2.47f, 2.45f, 2.42f, 2.41f, 2.41f, 2.42f, 2.40f, 2.39f, 2.38f, 2.37f, /* 61-70 */
    2.37f, 2.34f, 2.31f, 2.33f, 2.34f, 2.37f, 2.38f, 2.39f, 2.37f, 2.38f, /* 71-80 */
    2.41f, 2.41f, 2.42f, 2.40f, 2.47f, 2.60f, 2.54f, 2.52f, 2.50f, 2.47f  /* 81-90 */
};

struct AtomCoord {
    float x, y, z;
    int ele;
    float rvdw;
};

struct GridPoint {
    float x, y, z;
    int state;
};

static inline int cal_point_state(const struct AtomCoord *atoms, int n_atoms, float gx, float gy, float gz) {
    int fit = 0;
    int index_min = 0;
    float min_dis = 1e9f;

    for (int i = 0; i < n_atoms; i++) {
        if (atoms[i].ele > 1) {
            float dx = atoms[i].x - gx;
            float dy = atoms[i].y - gy;
            float dz = atoms[i].z - gz;
            min_dis = sqrtf(dx*dx + dy*dy + dz*dz);
            index_min = i;
            break;
        }
    }

    for (int i = 0; i < n_atoms; i++) {
        if (atoms[i].ele > 1) {
            float dx = atoms[i].x - gx;
            float dy = atoms[i].y - gy;
            float dz = atoms[i].z - gz;
            float d = sqrtf(dx*dx + dy*dy + dz*dz);
            if (d < min_dis) {
                min_dis = d;
                index_min = i;
            }
            if (d <= atoms[i].rvdw) {
                fit++;
            }
        }
    }

    if (fit == 0) {
        if (min_dis <= atoms[index_min].rvdw + 2.5f) return 0;
        else return 100;
    }
    if (fit == 1) return 1;
    if (fit > 1) return 2;
    return -1;
}

/*
 * Computes the 10-dimensional ComTarget C3DD descriptors for an input 3D molecule:
 * v1: Surface area
 * v2: Volume
 * v3: Surface / Volume ratio
 * v4: Hu moment invariant tt1
 * v5: Hu moment invariant tt2
 * v6: Hu moment invariant tt3
 * v7: Hu moment invariant tt4
 * v8: Shape ratio ks
 * v9: Sphericity kst
 * v10: Sphere fitting ratio kfit
 */
int compute_c3dd_descriptors(
    int n_atoms,
    const float *x,
    const float *y,
    const float *z,
    const int *atomic_num,
    float *out_descriptors_10
) {
    if (n_atoms <= 0 || !x || !y || !z || !atomic_num || !out_descriptors_10) {
        return -1;
    }

    struct AtomCoord *atoms = (struct AtomCoord*)malloc(n_atoms * sizeof(struct AtomCoord));
    if (!atoms) return -2;

    float xmin = 1e9f, xmax = -1e9f;
    float ymin = 1e9f, ymax = -1e9f;
    float zmin = 1e9f, zmax = -1e9f;
    int heavy_count = 0;

    for (int i = 0; i < n_atoms; i++) {
        atoms[i].x = x[i];
        atoms[i].y = y[i];
        atoms[i].z = z[i];
        atoms[i].ele = atomic_num[i];
        atoms[i].rvdw = (atomic_num[i] < 100 && atomic_num[i] > 0) ? ALVAREZ_VDW[atomic_num[i]] : 2.0f;

        if (atomic_num[i] > 1) {
            heavy_count++;
            if (x[i] < xmin) xmin = x[i];
            if (x[i] > xmax) xmax = x[i];
            if (y[i] < ymin) ymin = y[i];
            if (y[i] > ymax) ymax = y[i];
            if (z[i] < zmin) zmin = z[i];
            if (z[i] > zmax) zmax = z[i];
        }
    }

    if (heavy_count == 0) {
        free(atoms);
        return -3;
    }

    const float grid_gap = 0.3800f;
    int box_x = (int)((xmax - xmin + 16.0f) / grid_gap);
    int box_y = (int)((ymax - ymin + 16.0f) / grid_gap);
    int box_z = (int)((zmax - zmin + 16.0f) / grid_gap);
    long long total_grid = (long long)box_x * box_y * box_z;

    struct GridPoint *grid = (struct GridPoint*)malloc(total_grid * sizeof(struct GridPoint));
    int *wsp = (int*)calloc(total_grid, sizeof(int));
    if (!grid || !wsp) {
        if (grid) free(grid);
        if (wsp) free(wsp);
        free(atoms);
        return -4;
    }

    long long count = 0;
    float cur_x = xmin - 8.0f;
    for (int i = 0; i < box_x; i++) {
        float cur_y = ymin - 8.0f;
        for (int j = 0; j < box_y; j++) {
            float cur_z = zmin - 8.0f;
            for (int k = 0; k < box_z; k++) {
                grid[count].x = cur_x;
                grid[count].y = cur_y;
                grid[count].z = cur_z;
                count++;
                cur_z += grid_gap;
            }
            cur_y += grid_gap;
        }
        cur_x += grid_gap;
    }

    int number_of_no_blank = 0;
    for (long long i = 0; i < total_grid; i++) {
        grid[i].state = cal_point_state(atoms, n_atoms, grid[i].x, grid[i].y, grid[i].z);
        if (grid[i].state <= 2 && grid[i].state >= 0) {
            wsp[i] = 1;
            number_of_no_blank++;
        }
    }

    int *noblank_index = (int*)malloc(number_of_no_blank * sizeof(int));
    if (!noblank_index) {
        free(wsp);
        free(grid);
        free(atoms);
        return -5;
    }

    int nb_j = 0;
    for (long long i = 0; i < total_grid; i++) {
        if (wsp[i] == 1) {
            noblank_index[nb_j++] = (int)i;
        }
        if (grid[i].state == 100) {
            grid[i].state = 0;
        }
    }

    /* Convert candidate boundary points to interior */
    for (int i = 0; i < number_of_no_blank; i++) {
        if (grid[noblank_index[i]].state == 1) {
            grid[noblank_index[i]].state = 2;
        }
    }

    /* Boundary surface detection along X axis */
    for (int i = 0; i < box_z; i++) {
        for (int j = 0; j < box_y; j++) {
            for (int k = 0; k < box_x; k++) {
                long long cal_n = ((long long)k * box_y * box_z) + ((long long)j * box_z) + i;
                long long cal_nadd1 = ((long long)(k + 1) * box_y * box_z) + ((long long)j * box_z) + i;
                long long cal_nsub1 = ((long long)(k - 1) * box_y * box_z) + ((long long)j * box_z) + i;
                if (cal_nsub1 >= 0 && cal_nsub1 < total_grid && cal_nadd1 >= 0 && cal_nadd1 < total_grid) {
                    if (grid[cal_nadd1].state != 1 && grid[cal_nsub1].state != 1) {
                        if (grid[cal_nadd1].state != grid[cal_nsub1].state) {
                            grid[cal_n].state = 1;
                        }
                    }
                }
            }
        }
    }

    /* Boundary surface detection along Y axis */
    for (int i = 0; i < box_x; i++) {
        for (int j = 0; j < box_z; j++) {
            for (int k = 0; k < box_y; k++) {
                long long cal_n = ((long long)i * box_y * box_z) + ((long long)k * box_z) + j;
                long long cal_nadd1 = ((long long)i * box_y * box_z) + ((long long)(k + 1) * box_z) + j;
                long long cal_nsub1 = ((long long)i * box_y * box_z) + ((long long)(k - 1) * box_z) + j;
                if (cal_nsub1 >= 0 && cal_nsub1 < total_grid && cal_nadd1 >= 0 && cal_nadd1 < total_grid) {
                    if (grid[cal_nadd1].state != 1 && grid[cal_nsub1].state != 1) {
                        if (grid[cal_nadd1].state != grid[cal_nsub1].state) {
                            grid[cal_n].state = 1;
                        }
                    }
                }
            }
        }
    }

    /* Boundary surface detection along Z axis */
    for (int i = 0; i < box_x; i++) {
        for (int j = 0; j < box_y; j++) {
            for (int k = 0; k < box_z; k++) {
                long long cal_n = ((long long)i * box_y * box_z) + ((long long)j * box_z) + k;
                long long cal_nadd1 = ((long long)i * box_y * box_z) + ((long long)j * box_z) + k + 1;
                long long cal_nsub1 = ((long long)i * box_y * box_z) + ((long long)j * box_z) + k - 1;
                if (cal_nsub1 >= 0 && cal_nsub1 < total_grid && cal_nadd1 >= 0 && cal_nadd1 < total_grid) {
                    if (grid[cal_nadd1].state != 1 && grid[cal_nsub1].state != 1) {
                        if (grid[cal_nadd1].state != grid[cal_nsub1].state) {
                            grid[cal_n].state = 1;
                        }
                    }
                }
            }
        }
    }

    int surf_count = 0;
    for (int i = 0; i < number_of_no_blank; i++) {
        if (grid[noblank_index[i]].state == 1) {
            surf_count++;
        }
    }

    float surface = 1.26100f * (surf_count * grid_gap * grid_gap) + 13.50600f;
    float volume = number_of_no_blank * grid_gap * grid_gap * grid_gap;
    out_descriptors_10[0] = surface;
    out_descriptors_10[1] = volume;
    out_descriptors_10[2] = (volume > 0.0f) ? (surface / volume) : 0.0f;

    /* 3D Hu Moments on surface points (state == 1) */
    double m000 = 0.0, m100 = 0.0, m010 = 0.0, m001 = 0.0;
    for (long long i = 0; i < total_grid; i++) {
        if (grid[i].state == 1) {
            m000 += 1.0;
            m100 += grid[i].x;
            m010 += grid[i].y;
            m001 += grid[i].z;
        }
    }

    if (m000 < 1.0) m000 = 1.0;
    double xc = m100 / m000;
    double yc = m010 / m000;
    double zc = m001 / m000;

    double u200 = 0, u020 = 0, u002 = 0, u110 = 0, u101 = 0, u011 = 0;
    double u300 = 0, u030 = 0, u003 = 0, u120 = 0, u102 = 0, u012 = 0, u210 = 0, u021 = 0, u201 = 0, u111 = 0;

    for (long long i = 0; i < total_grid; i++) {
        if (grid[i].state == 1) {
            double dx = grid[i].x - xc;
            double dy = grid[i].y - yc;
            double dz = grid[i].z - zc;
            u200 += dx * dx;
            u020 += dy * dy;
            u002 += dz * dz;
            u110 += dx * dy;
            u101 += dx * dz;
            u011 += dy * dz;
            u300 += dx * dx * dx;
            u030 += dy * dy * dy;
            u003 += dz * dz * dz;
            u120 += dx * dy * dy;
            u102 += dx * dz * dz;
            u012 += dy * dz * dz;
            u210 += dx * dx * dy;
            u021 += dy * dy * dz;
            u201 += dx * dx * dz;
            u111 += dx * dy * dz;
        }
    }

    double u000 = m000;
    double norm2 = pow(u000, 5.0 / 3.0);
    double norm3 = pow(u000, 2.0);

    u200 /= norm2; u020 /= norm2; u002 /= norm2;
    u110 /= norm2; u101 /= norm2; u011 /= norm2;
    u300 /= norm3; u030 /= norm3; u003 /= norm3;
    u120 /= norm3; u102 /= norm3; u012 /= norm3;
    u210 /= norm3; u021 /= norm3; u201 /= norm3;
    u111 /= norm3;

    double tt1 = u200 + u020 + u002;
    double tt2 = (u200*u020) + (u200*u002) + (u020*u002) - (u110*u110) - (u101*u101) - (u011*u011);
    double tt3 = (u200*u020*u002) + (2.0*u110*u101*u011) - (u200*u011*u011) - (u020*u101*u101) - (u002*u110*u110);
    double tt4 = (u300*u300) + (u030*u030) + (u003*u003) + 3.0*(u120*u120 + u102*u102 + u012*u012 + u210*u210 + u021*u021 + u201*u201) + 6.0*(u111*u111);

    out_descriptors_10[3] = (float)tt1;
    out_descriptors_10[4] = (float)tt2;
    out_descriptors_10[5] = (float)tt3;
    out_descriptors_10[6] = (float)tt4;

    double ks = pow(3.0, 5.0/3.0) / (5.0 * pow(4.0 * M_PI, 2.0/3.0)) * (1.0 / (tt1 != 0.0 ? tt1 : 1.0));
    out_descriptors_10[7] = (float)ks;

    double kst = (36.0 * M_PI * volume * volume) / (surface * surface * surface);
    out_descriptors_10[8] = (float)kst;

    /* Centroid of entire non-blank molecular space for equivalent sphere fitting */
    long long dot_in_mol = 0;
    float mcx = 0, mcy = 0, mcz = 0;
    for (long long i = 0; i < total_grid; i++) {
        if (grid[i].state == 2 || grid[i].state == 1) {
            dot_in_mol++;
            mcx += grid[i].x;
            mcy += grid[i].y;
            mcz += grid[i].z;
        }
    }

    if (dot_in_mol > 0) {
        mcx /= dot_in_mol;
        mcy /= dot_in_mol;
        mcz /= dot_in_mol;
    }

    double vr = pow((volume * 3.0 / (4.0 * M_PI)), (1.0 / 3.0));
    long long dot_in_sphere = 0;
    for (long long i = 0; i < total_grid; i++) {
        if (grid[i].state == 2 || grid[i].state == 1) {
            float dx = grid[i].x - mcx;
            float dy = grid[i].y - mcy;
            float dz = grid[i].z - mcz;
            if (sqrtf(dx*dx + dy*dy + dz*dz) <= vr) {
                dot_in_sphere++;
            }
        }
    }

    double vinsph = dot_in_sphere * grid_gap * grid_gap * grid_gap;
    double denom = (2.0 * volume) - vinsph;
    double kfit = (denom != 0.0) ? (vinsph / denom) : 0.0;
    out_descriptors_10[9] = (float)kfit;

    free(noblank_index);
    free(wsp);
    free(grid);
    free(atoms);
    return 0;
}
