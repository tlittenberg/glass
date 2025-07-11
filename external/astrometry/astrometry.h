/*
 # The below functions are from the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license
 */

#include <stdbool.h>
#include <stdint.h>

#ifndef ASTROMETRY_H
#define ASTROMETRY_H

// Internal type
struct hp_s {
    int bighp;
    int x;
    int y;
};
typedef struct hp_s hp_t;

int64_t healpixl_compose_xy(int bighp, int x, int y, int Nside);
void healpixl_decompose_ring(int64_t hp, int Nside, int* p_ring, int* p_longind);
const int64_t healpixl_ring_to_xy(int64_t ring, int Nside);
void healpixl_decompose_xy(int64_t finehp,
                           int* pbighp, int* px, int* py,
                           int Nside);
void healpixl_to_radec(int64_t ihp, int Nside, double dx, double dy,
                       double* ra, double* dec);
void xyzarr2radec(const double* xyz, double *ra, double *dec);
const double xy2ra(double x, double y);
const double z2dec(double z);
#endif // ASTROMETRY_H
