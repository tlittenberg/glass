/*
 # The below functions are from the Astrometry.net suite.
 # Licensed under a 3-clause BSD style license - see LICENSE
 */

#include <stdbool.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include "astrometry.h"

static inline void swap_double(double* i1, double* i2) {
    double tmp;
    tmp = *i1;
    *i1 = *i2;
    *i2 = tmp;
}

/*
static inline bool isequatorial(int healpix)
{
    // the north polar healpixes are 0,1,2,3
    // the south polar healpixes are 8,9,10,11
    return (healpix >= 4) && (healpix <= 7);
}
*/

static inline bool isnorthpolar(int healpix)
{
    return (healpix <= 3);
}

static inline bool issouthpolar(int healpix)
{
    return (healpix >= 8);
}

static void hp_decompose(hp_t* hp, int* php, int* px, int* py) {
    if (php)
        *php = hp->bighp;
    if (px)
        *px = hp->x;
    if (py)
        *py = hp->y;
}
static void hp_to_xyz(hp_t* hp, int Nside,
                      double dx, double dy,
                      double* rx, double *ry, double *rz) {
    int chp;
    bool equatorial = true;
    double zfactor = 1.0;
    int xp, yp;
    double x, y, z;
    double pi = M_PI, phi;
    double rad;

    hp_decompose(hp, &chp, &xp, &yp);

    // this is x,y position in the healpix reference frame
    x = xp + dx;
    y = yp + dy;

    if (isnorthpolar(chp)) {
        if ((x + y) > Nside) {
            equatorial = false;
            zfactor = 1.0;
        }
    }
    if (issouthpolar(chp)) {
        if ((x + y) < Nside) {
            equatorial = false;
            zfactor = -1.0;
        }
    }

    if (equatorial) {
        double zoff=0;
        double phioff=0;
        x /= (double)Nside;
        y /= (double)Nside;

        if (chp <= 3) {
            // north
            phioff = 1.0;
        } else if (chp <= 7) {
            // equator
            zoff = -1.0;
            chp -= 4;
        } else if (chp <= 11) {
            // south
            phioff = 1.0;
            zoff = -2.0;
            chp -= 8;
        } else {
            // should never get here
            assert(0);
        }

        z = 2.0/3.0*(x + y + zoff);
        phi = pi/4*(x - y + phioff + 2*chp);
        rad = sqrt(1.0 - z*z);

    } else {
        /*
         Rearrange eqns (19) and (20) to find phi_t in terms of x,y.

         y = Ns - k in eq(19)
         x - Ns - k in eq(20)

         (Ns - y)^2 / (Ns - x)^2 = (2 phi_t)^2 / (2 phi_t - pi)^2

         Recall than y<=Ns, x<=Ns and 0<=phi_t<pi/2, so we can choose the
         root we want by taking square roots:

         (Ns - y) (pi - 2 phi_t) = 2 phi_t (Ns - x)
         (Ns - y) pi = 2 phi_t (Ns - x + Ns - y)
         phi_t = pi (Ns-y) / (2 (Ns - x) + (Ns - y))
         */
        double phi_t;
        double vv;

        if (zfactor == -1.0) {
            swap_double(&x, &y);
            x = (Nside - x);
            y = (Nside - y);
        }

        if (y == Nside && x == Nside)
            phi_t = 0.0;
        else
            phi_t = pi * (Nside-y) / (2.0 * ((Nside-x) + (Nside-y)));

        if (phi_t < pi/4.) {
            // z = 1.0 - mysquare(pi * (Nside - x) / ((2.0 * phi_t - pi) * Nside)) / 3.0;
            vv = fabs(pi * (Nside - x) / ((2.0 * phi_t - pi) * Nside) / sqrt(3));
        } else {
            // z = 1.0 - mysquare(pi * (Nside - y) / (2.0 * phi_t * Nside)) / 3.0;
            vv = fabs(pi * (Nside-y) / (2. * phi_t * Nside) / sqrt(3));
        }
        z = (1 - vv) * (1 + vv);
        rad = sqrt(1.0 + z) * vv;

        assert(0.0 <= fabs(z) && fabs(z) <= 1.0);
        z *= zfactor;
        assert(0.0 <= fabs(z) && fabs(z) <= 1.0);

        // The big healpix determines the phi offset
        if (issouthpolar(chp))
            phi = pi/2.0* (chp-8) + phi_t;
        else
            phi = pi/2.0 * chp + phi_t;
    }

    if (phi < 0.0)
        phi += 2*pi;

    *rx = rad * cos(phi);
    *ry = rad * sin(phi);
    *rz = z;
}

static void intltohp(int64_t pix, hp_t* hp, int Nside) {
    healpixl_decompose_xy(pix, &hp->bighp, &hp->x, &hp->y, Nside);
}

void healpixl_decompose_xy(int64_t finehp,
                           int* pbighp, int* px, int* py,
                           int Nside) {
    int64_t hp;
    int64_t ns2 = (int64_t)Nside * (int64_t)Nside;
    assert(Nside > 0);
    assert(finehp < ((int64_t)12 * ns2));
    assert(finehp >= 0);
    if (pbighp) {
        int bighp = (int)(finehp / ns2);
        assert(bighp >= 0);
        assert(bighp < 12);
        *pbighp = bighp;
    }
    hp = finehp % ns2;
    if (px) {
        *px = (int)(hp / Nside);
        assert(*px >= 0);
        assert(*px < Nside);
    }
    if (py) {
        *py = hp % Nside;
        assert(*py >= 0);
        assert(*py < Nside);
    }
}

const double xy2ra(double x, double y) {
	double a = atan2(y, x);
	if (a < 0)
		a += 2.0 * M_PI;
	return a;
}

const double z2dec(double z) {
	return asin(z);
}

void xyz2radec(double x, double y, double z, double *ra, double *dec) {
    if (ra)
    	*ra = xy2ra(x, y);
    if (dec) {
        if (fabs(z) > 0.9)
            *dec = M_PI / 2.0 - atan2(hypot(x, y), z);
        else
            *dec = z2dec(z);
    }
}

void xyzarr2radec(const double* xyz, double *ra, double *dec) {
	xyz2radec(xyz[0], xyz[1], xyz[2], ra, dec);
}

void healpixl_to_radec(int64_t ihp, int Nside, double dx, double dy,
                       double* ra, double* dec) {
    hp_t hp;
    double xyz[3];
    intltohp(ihp, &hp, Nside);
    hp_to_xyz(&hp, Nside, dx, dy, xyz, xyz+1, xyz+2);
    xyzarr2radec(xyz, ra, dec);
}

int64_t healpixl_compose_xy(int bighp, int x, int y, int Nside) {
    int64_t ns = Nside;
    assert(Nside > 0);
    assert(bighp >= 0);
    assert(bighp < 12);
    assert(x >= 0);
    assert(x < Nside);
    assert(y >= 0);
    assert(y < Nside);
    return ((((int64_t)bighp * ns) + x) * ns) + y;
}

void healpixl_decompose_ring(int64_t hp, int Nside, int* p_ring, int* p_longind) {
    int64_t longind;
    int64_t offset = 0;
    int64_t Nside64;
    int64_t ns2;
    int ring;
    double x;
    Nside64 = (int64_t)Nside;
    ns2 = Nside64 * Nside64;
    if (hp < 2 * ns2) {
        ring = (int)(0.5 + sqrt(0.25 + 0.5 * hp));
        offset = 2 * (int64_t)ring * ((int64_t)ring - 1);
        // The sqrt above can introduce precision issues that can cause ring to
        // be off by 1, so we check whether the offset is now larger than the HEALPix
        // value, and if so we need to adjust ring and offset accordingly
        if (offset > hp) {
            ring -= 1;
            offset = 2 * (int64_t)ring * ((int64_t)ring - 1);
        }
        longind = hp - offset;
    } else {
        offset = 2 * Nside64 * (Nside64 - 1);
        if (hp < 10 * ns2) {
            ring = (int)((hp - offset) / ((int64_t)Nside * 4) + (int64_t)Nside);
            offset += 4 * (ring - Nside64) * Nside64;
            longind = hp - offset;
        } else {
            offset += 8 * ns2;
            x = (2 * Nside64 + 1 - sqrt((2 * Nside64 + 1) * (2 * Nside64 + 1) - 2 * (hp - offset)))*0.5;
            ring = (int)x;
            offset += 2 * (int64_t)ring * (2 * Nside64 + 1 - (int64_t)ring);
            // The sqrt above can introduce precision issues that can cause ring to
            // be off by 1, so we check whether the offset is now larger than the HEALPix
            // value, and if so we need to adjust ring and offset accordingly
            if (offset > hp) {
                ring -= 1;
                offset -= 4 * Nside64 - 4 * (int64_t)ring;
            }
            longind = (int)(hp - offset);
            ring += 3 * Nside;
        }
    }
    if (p_ring)
        *p_ring = ring;
    if (p_longind)
        *p_longind = (int)longind;
}

const int64_t healpixl_ring_to_xy(int64_t ring, int Nside) {
    int bighp, x, y;
    int ringind, longind;
    healpixl_decompose_ring(ring, Nside, &ringind, &longind);
    if (ring < 0 || Nside < 0) {
        return -1;
    } else if (ringind <= Nside) {
        int64_t ind;
        int v;
        int F1;
        int frow;
        bighp = longind / ringind;
        ind = (int64_t)longind - (int64_t)bighp * (int64_t)ringind;
        y = (Nside - 1) - (int)ind;
        frow = bighp / 4;
        F1 = frow + 2;
        v = F1*Nside - ringind - 1;
        x = v - y;
        return healpixl_compose_xy(bighp, x, y, Nside);
    } else if (ringind < (int64_t)3*Nside) {
        int panel;
        int ind;
        int bottomleft;
        int topleft;
        int frow, F1, F2, s, v, h;
        int bighp = -1;
        int x, y;
        int R = 0;

        panel = longind / Nside;
        ind = longind % Nside;
        bottomleft = ind < (ringind - Nside + 1) / 2;
        topleft = ind < ((int64_t)3*Nside - ringind + 1)/2;

        if (!bottomleft && topleft) {
            // top row.
            bighp = panel;
        } else if (bottomleft && !topleft) {
            // bottom row.
            bighp = 8 + panel;
        } else if (bottomleft && topleft) {
            // left side.
            bighp = 4 + panel;
        } else if (!bottomleft && !topleft) {
            // right side.
            bighp = 4 + (panel + 1) % 4;
            if (bighp == 4) {
                longind -= ((int64_t)4*Nside - 1);
                // Gah!  Wacky hack - it seems that since
                // "longind" is negative in this case, the
                // rounding behaves differently, so we end up
                // computing the wrong "h" and have to correct
                // for it.
                R = 1;
            }
        }

        frow = bighp / 4;
        F1 = frow + 2;
        F2 = 2*(bighp % 4) - (frow % 2) + 1;
        s = (ringind - Nside) % 2;
        v = F1*Nside - ringind - 1;
        h = 2*longind - s - F2*Nside;
        if (R)
            h--;
        x = (v + h) / 2;
        y = (v - h) / 2;
        //fprintf(stderr, "bighp=%i, frow=%i, F1=%i, F2=%i, s=%i, v=%i, h=%i, x=%i, y=%i.\n", bighp, frow, F1, F2, s, v, h, x, y);

        if ((v != (x+y)) || (h != (x-y))) {
            h++;
            x = (v + h) / 2;
            y = (v - h) / 2;
            //fprintf(stderr, "tweak h=%i, x=%i, y=%i\n", h, x, y);

            if ((v != (x+y)) || (h != (x-y))) {
                //fprintf(stderr, "still not right.\n");
            }
        }
        return healpixl_compose_xy(bighp, x, y, Nside);
    } else {
        int ind;
        int v;
        int F1;
        int frow;
        int ri;
        ri = 4*Nside - ringind;
        bighp = 8 + longind / ri;
        ind = longind - (bighp%4) * ri;
        y = (ri-1) - ind;
        frow = bighp / 4;
        F1 = frow + 2;
        v = F1*Nside - ringind - 1;
        x = v - y;
        return healpixl_compose_xy(bighp, x, y, Nside);
    }
}
