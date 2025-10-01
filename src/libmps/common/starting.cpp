/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Dario Andrea Bini <bini@dm.unipi.it>
*   Giuseppe Fiorentino <fiorent@dm.unipi.it>
*   Leonardo Robol <leonardo.robol@unipi.it>
*/


#include <math.h>
#include <float.h>
#include <limits.h>
#include <assert.h>
#include <mps/mps.h>

#define MPS_STARTING_SIGMA (0.66 * (PI / ctx->n))
#define pi2 6.283184

/**
* @brief Compute the greatest common divisor of <code>a</code>
* and <code>b</code>.
*
* @param a first integer
* @param b second integer
* @return <code>GCD(a,b)</code>
*/
static int
    mps_gcd(int a, int b)
{
    int temp;

    do
    {
        temp = b;
        b = a % b;
        a = temp;
    } while (b != 0);
    return a;
}

/**
* @brief Find the sigma that maximize distance between
* starting approximation in the last annulus and the one
* in the new annulus.
*
* This function also set <code>ctx->last_sigma</code> to the value
* that is computed.
*
* @param ctx the mps_context struct pointer.
* @param last_sigma the last value of sigma.
* @param cluster_item The element of the <code>mps_clusterization</code> of which
* we are computing the starting points, or NULL if we are computing the starting
* points for all the approximations.
* @param n the number of roots in the cluster.
* @return the shift advised for the starting approximation in the annulus.
*/
static double
    mps_maximize_distance(mps_context* ctx, double last_sigma,
        mps_cluster_item* cluster_item, int n)
{
    double delta_sigma;

    /* Find number of roots in the last cluster */
    if (!cluster_item || cluster_item->prev_cluster_item == NULL)
        return ctx->last_sigma;
    int old_clust_n = (cluster_item->prev_cluster_item->cluster->n);

    /* Compute right shifting angle for the new approximations, i.e.
    * pi / [m,n] where [m,n] is the least common multiply of m and n.
    * This is done by computing the gcd and then dividing m*n by it. */
    delta_sigma = PI * (old_clust_n * mps_gcd(old_clust_n, n)) / (4 * n);

    /* Return shifted value, archiving it for the next pass */
    ctx->last_sigma = last_sigma + delta_sigma;
    return ctx->last_sigma;
}

/**
* @brief Compute radii of the circles where the initial approximation
* will be disposed by mps_fstart()
*
* @param ctx mps_context* stuct pointer.
* @param n number of roots in the cluster.
* @param cluster_item The element of the <code>mps_clusterization</code> of which
* we are computing the starting points, or NULL if we are computing the starting
* points for all the approximations.
* @param clust_rad radius of the cluster.
* @param g new gravity center where the polynomial has been shifted.
* @param eps out epsilon.
* @param fap[] Array with the moduli of the coefficients.
*
* @return A mps_starting_configuration containing the necessary information
* to dispose the initial approximations on the circles.
*
* @see mps_fstart()
*/
MPS_PRIVATE mps_starting_configuration
    mps_fcompute_starting_radii(mps_context* ctx, int n, mps_cluster_item* cluster_item,
        double clust_rad, double g, rdpe_t eps,
        double fap[])
{
    MPS_DEBUG_THIS_CALL(ctx);

    const double big_double = DBL_MAX, small_double = DBL_MIN;
    const double xbig = ::log(big_double), xsmall = ::log(small_double);

    int i, j, k, nzeros, iold, ni, offset;
    double temp, r;

    mps_starting_configuration cfg;
    cfg.n_radii = 0;
    cfg.partitioning = NULL;
    cfg.fradii = NULL;
    cfg.dradii = NULL;

    ni = 0;
    nzeros = 0;
    r = 0.0;

    /**********************************************
        check for possible null entries in the trailing
        coefficients only in the case where the polynomial
        has been shifted in g, replace null coefficients
        with small numbers according to the working precision
        and to the number of null coefficients
    **********************************************/
    if (g != 0.0)
    {
        for (i = 0; i <= n; i++)
            if (fap[i] != 0.0)
            {
                ni = i;
                break;
            }
        if (ni == 0)
            temp = 2 * xsmall;
        else
            temp = ::log(fap[ni]) + ni * (::log(DBL_EPSILON) + ::log(g * ni * 10.0));
    }
    else
        temp = 2 * xsmall;

    for (i = 0; i <= n; i++)
    {
        if (fap[i] != 0.0)
            ctx->fap2[i] = ::log(fap[i]);
        else
            ctx->fap2[i] = temp;
    }

    /* Compute convex hull */
    int* h = mps_fconvex(ctx, n, ctx->fap2);

    j = 0;
    for (i = 1; i <= n; i++)
        if (h[i])
            j++;

    cfg.fradii = mps_newv(double, j + 2);
    cfg.partitioning = mps_newv(int, j + 2);

    /* compute the radii of the circles containing starting approximations  */
    cfg.n_radii = 0;
    cfg.partitioning[0] = 0;

    for (i = 1; i <= n; i++)
        if (h[i])
        {
            iold = cfg.partitioning[cfg.n_radii];
            nzeros = i - iold;
            temp = (ctx->fap2[iold] - ctx->fap2[i]) / nzeros;

            /* if the radius is too small to be represented as double, set it
            * to the minimum  representable double */
            if (temp < xsmall)      /* if (temp < MAX(xsmall, -xbig)) DARIO Giugno 23 */
                r = DBL_MIN;          /* r = small; */

            /* if the radius is too big to be represented as double, set it
            * to the maximum representable double */
            if (temp > xbig)
                r = DBL_MAX;          /* big;   DARIO Giugno 23 */

            /* if the radius is representable as double, compute it    */
            if ((temp <= xbig) && (temp > xsmall))
                /* if ((temp <= xbig) && (temp > MAX(-xbig, xsmall))) DARIO Giugno 23 */
                r = exp(temp);

            /* if the radius is greater than the radius of the cluster
            * set the radius equal to the radius of the cluster */
            if (clust_rad != 0 && r > clust_rad)
                r = clust_rad;

            cfg.fradii[cfg.n_radii] = r;
            cfg.partitioning[++cfg.n_radii] = i;
        }

    /* Close partitioning */
    cfg.partitioning[cfg.n_radii] = n;

    /* Compact radius that are too near */
    for (i = 0; i < cfg.n_radii; i++)
    {
        /* Scan next radii to see if they are near the
        * i-th that we are considering now  */
        for (j = i + 1; j < cfg.n_radii; j++)
        {
            /* Get an estimate of the distance between approximations that would
            * be obtained collapsing the circles. It's a little understimated
            * but this is due to the fact that changing approximation is only
            * a fallback, and not the preferred action to perform. */
            if (fabs((cfg.fradii[j] - cfg.fradii[i])) >
                MIN(cfg.fradii[j],
                    cfg.fradii[i]) * PI / (cfg.partitioning[j + 1] -
                        cfg.partitioning[i]))
            {
                break;
            }
        }

        /* This is the number of circles that are near */
        j--;
        offset = j - i;

        /* If there is nothing to compact, do not compact it */
        if (offset == 0)
        {
            continue;
        }

        if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            MPS_DEBUG(ctx, "Compacting circles from %d to %d", i, j);

        /* We shall now compact circles between i and j, so
        * we start computing the mean of the radius */
        cfg.fradii[i] *= (cfg.partitioning[i + 1] - cfg.partitioning[i]);
        for (k = i + 1; k <= j; k++)
        {
            cfg.fradii[i] +=
                cfg.fradii[k] * (cfg.partitioning[k + 1] - cfg.partitioning[k]);
        }

        cfg.fradii[i] /= cfg.partitioning[j + 1] - cfg.partitioning[i];

        /* Move other circles backward */
        for (k = j + 1; k < cfg.n_radii; k++)
        {
            cfg.fradii[k - offset] = cfg.fradii[k];
            cfg.partitioning[k - offset] = cfg.partitioning[k];
        }

        /* Set new cfg.n_radii and new partitioning */
        cfg.n_radii = cfg.n_radii - offset;
        cfg.partitioning[cfg.n_radii] = n;
    }

    mps_free(h);

    return cfg;
}

/**
* @brief Compute new starting approximations to the roots
* of the polynomial \f$p(x)\f$ having coefficients of modulus apoly.
*
* Computations is done by
* means of the Rouche'-based criterion of Bini (Numer. Algo. 1996).
* The program can compute all the approximations
* (if \f$n\f$ is the degree of \f$p(x)\f$) or it may compute the
* approximations of the cluster in the <code>cluster_item</code>.
* The status vector is changed into <code>'o'</code> for the components
* that belong to a cluster with relative radius less than <code>eps</code>.
* The status vector is changed into <code>'x'</code> for the components that
* cannot be represented as double.
*
* @param ctx The <code>mps_context</code> associated with the current computation.
* @param n number of roots in the cluster.
* @param cluster_item The element of the <code>mps_clusterization</code> of which
* we are computing the starting points, or NULL if we are computing the starting
* points for all the approximations.
* @param clust_rad radius of cluster.
* @param g gravity center of the cluster.
* @param eps a double that represent the maximum value
* of relative radius (with respect to <code>g</code>) of
* roots whose status must be set to <code>o</code>.
* @param fap array of moduli of the coefficients as double.
*
* @see status
*/
MPS_PRIVATE void
    mps_fstart(mps_context* ctx, int n, mps_cluster_item* cluster_item,
        double clust_rad, double g, rdpe_t eps, double fap[])
{
    MPS_DEBUG_THIS_CALL(ctx);
    int i, j, jj, l, nzeros = 0;
    double sigma, th, ang, r = 0;
    mps_root* it_root = NULL;
    rdpe_t tmp;

    mps_cluster* cluster = NULL;
    mps_root* root = NULL;

    if (cluster_item)
        cluster = cluster_item->cluster;

    if (ctx->random_seed)
        sigma = drand();
    else
    {
        /* If this is the first cluster select sigma = 0. In the other
        * case try to maximize starting points distance. */
        if ((cluster == NULL) || (cluster_item->prev_cluster_item == NULL))
        {
            sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
        }
        else
        {
            sigma = mps_maximize_distance(ctx, ctx->last_sigma, cluster_item, n);
        }
    }

    th = pi2 / n;

    /* In the general case apply the Rouche-based criterion */
    mps_starting_configuration cfg =
        mps_fcompute_starting_radii(ctx, n, cluster_item, clust_rad, g, eps, fap);

    if (g != 0.0 && cluster != NULL)
        it_root = cluster->first_root;

    for (i = 0; i < cfg.n_radii; i++)
    {
        nzeros = cfg.partitioning[i + 1] - cfg.partitioning[i];
        ang = pi2 / nzeros;
        r = cfg.fradii[i];

        for (j = cfg.partitioning[i]; j < cfg.partitioning[i + 1]; j++)
        {
            if (g != 0.0 && it_root)
                l = it_root->k;
            else
                l = j;

            jj = j - cfg.partitioning[i];

            /* if the radius reaches extreme values then set the component
            * of status, corresponding to approximation which fall out the
            * representable range, to 'x' (out)    */
            if ((r == DBL_MIN) || (r == DBL_MAX))
                /* if ((r == small) || (r == big)) DARIO Giugno 23 */
                ctx->approx_root[l]->status = MPS_ROOT_STATUS_NOT_FLOAT;
            cplx_set_d(ctx->approx_root[l]->fvalue,
                r * cos(ang * jj + th * cfg.partitioning[i + 1] +
                    sigma),
                r * sin(ang * jj + th * cfg.partitioning[i + 1] +
                    sigma));
            if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
                MPS_DEBUG_CPLX(ctx, ctx->approx_root[l]->fvalue, "ctx->froot[%d]", l);
            }

            if (it_root)
                it_root = it_root->next_root;
        }


        /* If the new radius of the cluster is relatively smaller, then
        * set the status component equal to 'o' (output) */
        if (g != 0.0 && cluster)
        {
            rdpe_mul_d(tmp, eps, g);
            if (r * nzeros <= rdpe_get_d(tmp))
            {
                for (root = cluster->first_root; root != NULL; root = root->next_root)
                {
                    l = root->k;
                    ctx->approx_root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
                    ctx->approx_root[l]->frad = r * nzeros;
                }
            }
        }
    }

    mps_starting_configuration_clear(ctx, &cfg);
}

/**
    * @brief Compute radii of the circles where the initial approximation
    * will be disposed by mps_dstart()
    *
    * @param ctx mps_context* stuct pointer.
    * @param n number of roots in the cluster.
    * @param cluster_item The element of the <code>mps_clusterization</code> of which
    * we are computing the starting points, or NULL if we are computing the starting
    * points for all the approximations.
    * @param clust_rad radius of the cluster.
    * @param g new gravity center where the polynomial has been shifted.
    * @param eps out epsilon.
    * @param dap[] Array with the moduli of the coefficients.
    *
    * @see mps_dstart()
    */
MPS_PRIVATE mps_starting_configuration
    mps_dcompute_starting_radii(mps_context * ctx, int n,
        mps_cluster_item * cluster_item,
        rdpe_t clust_rad, rdpe_t g, rdpe_t eps,
        rdpe_t dap[])
{
    /* Compute big and small values */
    double xbig, xsmall;

    int i, j, k, nzeros, ni, offset, iold;
    double temp;
    rdpe_t r, tmp;

    xbig = LONG_MAX;
    xsmall = LONG_MIN;

    mps_starting_configuration cfg;
    cfg.n_radii = 0;
    cfg.partitioning = NULL;
    cfg.fradii = NULL;
    cfg.dradii = NULL;

    ni = 0;
    nzeros = 0;
    rdpe_set(r, rdpe_zero);

    if (rdpe_ne(g, rdpe_zero))
    {
        for (i = 0; i <= n; i++)
            if (rdpe_ne(dap[i], rdpe_zero))
            {
                ni = i;
                break;
            }
        if (ni == 0)
            temp = -2.0 * (LONG_MAX * LOG2);
        else
        {
            /* temp = log(dap[ni])+ni*(log(DBL_EPSILON)+log(g*ni*10.d0)) */
            temp = ni * 10.0;
            rdpe_mul_d(tmp, g, temp);
            temp = rdpe_log(tmp);
            temp += ::log(DBL_EPSILON);
            temp *= ni;
            temp += rdpe_log(dap[ni]);
        }
    }
    else
        temp = -2.0 * (LONG_MAX * LOG2);
    for (i = 0; i <= n; i++)
        if (rdpe_ne(dap[i], rdpe_zero))
            ctx->fap2[i] = rdpe_log(dap[i]);
        else
            ctx->fap2[i] = temp;

    /* compute the convex hull */
    int* h = mps_fconvex(ctx, n, ctx->fap2);

    /* Count the number of vertexes in the convex hull. */
    j = 0;
    for (i = 1; i <= n; i++)
        if (h[i])
            j++;

    cfg.dradii = mps_newv(rdpe_t, j + 2);
    cfg.partitioning = mps_newv(int, j + 2);

    /* compute the radii of the circles containing starting approximations  */
    cfg.n_radii = 0;
    cfg.partitioning[0] = 0;
    for (i = 1; i <= n; i++)
        if (h[i])
        {
            iold = cfg.partitioning[cfg.n_radii];
            nzeros = i - iold;
            temp = (ctx->fap2[iold] - ctx->fap2[i]) / nzeros;

            /* if the radius is too small to be represented as double, set it
                * to the minimum  representable double */
            if (temp < xsmall)
                rdpe_set(r, RDPE_MIN);       /* r = small; */

            /* if the radius is too big to be represented as double, set it
                * to the maximum representable double */
            if (temp > xbig)
                rdpe_set(r, RDPE_MAX);

            /* if the radius is representable as double, compute it    */
            if ((temp < xbig) && (temp > xsmall))
                rdpe_set_d(r, temp);
            rdpe_exp_eq(r);

            /* if the radius is greater than the radius of the cluster
                * set the radius equal to the radius of the cluster */
            if (rdpe_ne(clust_rad, rdpe_zero) && rdpe_gt(r, clust_rad))
                rdpe_set(r, clust_rad);

            rdpe_set(cfg.dradii[cfg.n_radii], r);
            cfg.partitioning[++cfg.n_radii] = i;
        }

    /* Close partitioning */
    cfg.partitioning[cfg.n_radii] = n;

    /* Compact radius that are too near */
    for (i = 0; i < cfg.n_radii; i++)
    {
        /* Scan next radii to see if they are near the
            * i-th that we are considering now  */
        for (j = i + 1; j < cfg.n_radii; j++)
        {
            /* Get an estimate of the distance between approximations that would
                * be obtained collapsing the circles. It's a little understimated
                * but this is due to the fact that changing approximation is only
                * a fallback, and not the preferred action to perform. */
            rdpe_sub(tmp, cfg.dradii[j], cfg.dradii[i]);
            rdpe_abs_eq(tmp);
            if (rdpe_lt(cfg.dradii[i], cfg.dradii[j]))
                rdpe_div_eq(tmp, cfg.dradii[i]);
            else
                rdpe_div_eq(tmp, cfg.dradii[j]);
            rdpe_div_eq_d(tmp, PI);
            rdpe_mul_eq_d(tmp, cfg.partitioning[j + 1] - cfg.partitioning[i]);
            if (rdpe_gt(tmp, rdpe_one))
            {
                break;
            }
        }

        /* This is the number of circles that are near */
        j--;
        offset = j - i;


        /* If there is nothing to compact, do not compact it */
        if (offset == 0)
        {
            continue;
        }

        MPS_DEBUG(ctx, "Compacting circles from %d to %d", i, j);

        /* We shall now compact circles between i and j, so
            * we start computing the mean of the radius */
        rdpe_mul_eq_d(cfg.dradii[i],
            cfg.partitioning[i + 1] - cfg.partitioning[i]);
        for (k = i + 1; k <= j; k++)
        {
            rdpe_mul_d(tmp, cfg.dradii[k],
                cfg.partitioning[k + 1] - cfg.partitioning[k]);
            rdpe_add_eq(cfg.dradii[i], tmp);
        }

        rdpe_div_eq_d(cfg.dradii[i],
            cfg.partitioning[j + 1] - cfg.partitioning[i]);

        /* Move other circles backward */
        for (k = j + 1; k < cfg.n_radii; k++)
        {
            rdpe_set(cfg.dradii[k - offset], cfg.dradii[k]);
            cfg.partitioning[k - offset] = cfg.partitioning[k];
        }

        /* Set new cfg.n_radii and new partitioning */
        cfg.n_radii = cfg.n_radii - offset;
        cfg.partitioning[cfg.n_radii] = n;
    }

    mps_free(h);
    return cfg;
}

/**
    * @brief Compute new starting approximations to the roots of the
    * polynomial \f$p(x)\f$ having coefficients of modulus apoly, by
    * means of the Rouche'-based criterion of Bini (Numer. Algo. 1996).
    *
    * The program can compute all the approximations
    * (if \f$n\f$ is the degree of \f$p(x)\f$) or it may compute the
    * approximations of the cluster of the <code>cluster_item</code>.
    * The status vector is changed into <code>'o'</code> for the components
    * that belong to a cluster with relative radius less than <code>eps</code>.
    * The status vector is changed into <code>'f'</code> for the components
    * that cannot be represented as <code>>dpe</code>.
    *
    * @param ctx mps_context struct pointer.
    * @param n number of root in the cluster to consider
    * @param cluster_item The item of the mps_clusterization containing the cluster
    * that we are analyzing.
    * @param clust_rad radius of the cluster.
    * @param g new center in which the the polynomial will be shifted.
    * @param eps maximum radius considered small enough to be negligible.
    * @param dap[] moduli of the coefficients as <code>dpe</code> numbers.
    */
MPS_PRIVATE void
    mps_dstart(mps_context * ctx, int n, mps_cluster_item * cluster_item,
        rdpe_t clust_rad, rdpe_t g, rdpe_t eps, rdpe_t dap[])
{
    int l = 0, i, j, jj, nzeros = 0;
    rdpe_t r, tmp, tmp1;
    double sigma, th, ang;
    mps_boolean flag = false;

    mps_cluster* cluster = NULL;
    mps_root* root = NULL;

    if (cluster_item)
        cluster = cluster_item->cluster;

    if (ctx->random_seed)
        sigma = drand();
    else
    {
        /* If this is the first cluster select sigma = 0. In the other
            * case try to maximize starting points distance. */
        if (!cluster_item || cluster_item->prev_cluster_item == NULL)
        {
            sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
        }
        else
        {
            sigma = mps_maximize_distance(ctx, ctx->last_sigma, cluster_item, n);
        }
    }

    /* In the case of user-defined polynomial choose as starting
        * approximations equispaced points in the unit circle. */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
    {
        ang = pi2 / n;
        for (i = 0; i < n; i++)
        {
            cdpe_set_d(ctx->approx_root[i]->dvalue, cos(ang * i + sigma),
                sin(ang * i + sigma));
        }
        return;
    }


    /* check if it is the case dpe_after_float, in this case set flag=true  */
    for (i = 0; i < n; i++)
    {
        flag = (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT);
        if (flag)
            break;
    }

    /* Compute starting radii with the Rouche based criterion */
    mps_starting_configuration cfg =
        mps_dcompute_starting_radii(ctx, n, cluster_item, clust_rad, g, eps, dap);

    th = pi2 / n;

    for (i = 0; i < cfg.n_radii; i++)
    {
        nzeros = cfg.partitioning[i + 1] - cfg.partitioning[i];
        ang = pi2 / nzeros;
        rdpe_set(r, cfg.dradii[i]);

        if (cluster_item)
            root = cluster->first_root;

        for (j = cfg.partitioning[i]; j < cfg.partitioning[i + 1]; j++)
        {
            if (cluster_item)
            {
                l = root->k;
                root = root->next_root;
            }
            else
                l = j;


            jj = j - cfg.partitioning[i];

            /* If dpe_after_float (i.e., flag is true) recompute the starting
                * values of only the approximations falling out of the range */
            if (flag)
            {
                if (ctx->approx_root[l]->status == MPS_ROOT_STATUS_NOT_FLOAT)
                {
                    cdpe_set_d(ctx->approx_root[l]->dvalue,
                        cos(ang * jj + th * cfg.partitioning[i + 1] +
                            sigma),
                        sin(ang * jj + th * cfg.partitioning[i + 1] +
                            sigma));
                    cdpe_mul_eq_e(ctx->approx_root[l]->dvalue, r);
                    ctx->approx_root[l]->status = MPS_ROOT_STATUS_CLUSTERED;
                    /*#G 27/4/98 if (rdpe_eq(r, big) || rdpe_eq(r, small)) */
                    if (rdpe_eq(r, RDPE_MIN) || rdpe_eq(r, RDPE_MAX))
                        ctx->approx_root[l]->status = MPS_ROOT_STATUS_NOT_DPE;
                }
            }
            else
            {
                /* else compute all the initial approximations */
                cdpe_set_d(ctx->approx_root[l]->dvalue,
                    cos(ang * jj + th * cfg.partitioning[i + 1] +
                        sigma),
                    sin(ang * jj + th * cfg.partitioning[i + 1] +
                        sigma));
                cdpe_mul_eq_e(ctx->approx_root[l]->dvalue, r);
                /*#G 27/4/98 if (rdpe_eq(r, big) || rdpe_eq(r, small)) */
                if (rdpe_eq(r, RDPE_MIN) || rdpe_eq(r, RDPE_MAX))
                    ctx->approx_root[l]->status = MPS_ROOT_STATUS_NOT_DPE;
            }

            if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            {
                MPS_DEBUG_CDPE(ctx, ctx->approx_root[l]->dvalue, "ctx->droot[%d]", l);
            }
        }

        /* If the new radius of the cluster is relatively small, then
            * set the status component equal to 'o' (output) */
        if (cluster_item)
        {
            rdpe_mul(tmp, g, eps);
            rdpe_mul_d(tmp1, r, (double)nzeros);
            if (rdpe_lt(tmp1, tmp))
            {
                for (root = cluster->first_root; root != NULL; root = root->next_root)
                {
                    l = root->k;
                    ctx->approx_root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
                    rdpe_set(ctx->approx_root[l]->drad, tmp1);
                }
            }
        }
    }

    mps_starting_configuration_clear(ctx, &cfg);
}

/**
    * @brief Compute radii of the circles where the initial approximation
    * will be disposed by mps_mstart()
    *
    * @param ctx mps_context* stuct pointer.
    * @param n number of roots in the cluster.
    * @param cluster_item The element of the <code>mps_clusterization</code> of which
    * we are computing the starting points, or NULL if we are computing the starting
    * points for all the approximations.
    * @param clust_rad radius of the cluster.
    * @param g new gravity center where the polynomial has been shifted.
    * @param dap[] Array with the moduli of the coefficients.
    *
    * @see mps_mstart()
    */
MPS_PRIVATE mps_starting_configuration
    mps_mcompute_starting_radii(mps_context * ctx, int n, mps_cluster_item * cluster_item,
        rdpe_t clust_rad, rdpe_t g, rdpe_t dap[])
{
    int i, offset, iold, nzeros, j, k;
    rdpe_t big_rdpe_t, small_rdpe_t, tmp;
    double xbig, xsmall, temp;

    xsmall = rdpe_log(RDPE_MIN);
    xbig = rdpe_log(RDPE_MAX);
    rdpe_set(small_rdpe_t, RDPE_MIN);
    rdpe_set(big_rdpe_t, RDPE_MAX);

    mps_starting_configuration cfg;
    cfg.n_radii = 0;
    cfg.partitioning = NULL;
    cfg.fradii = NULL;
    cfg.dradii = NULL;

    if (rdpe_eq(dap[0], rdpe_zero))
        ctx->fap2[0] = -ctx->mpwp * LOG2;

    /*  check for possible null entries in the trailing coefficients */
    for (i = 0; i <= n; i++)
        if (rdpe_ne(dap[i], rdpe_zero))
            ctx->fap2[i] = rdpe_log(dap[i]);
        else
            ctx->fap2[i] = ctx->fap2[0];

    /* compute the convex hull */
    int* h = mps_fconvex(ctx, n, ctx->fap2);

    /* Count the number of vertexes in the convex hull. */
    j = 0;
    for (i = 1; i <= n; i++)
        if (h[i])
            j++;

    cfg.dradii = mps_newv(rdpe_t, j + 2);
    cfg.partitioning = mps_newv(int, j + 2);

    /* Scan all the vertices of the convex hull */
    cfg.partitioning[0] = 0;
    cfg.n_radii = 0;
    for (i = 1; i <= n; i++)
    {
        if (h[i])
        {
            iold = cfg.partitioning[cfg.n_radii];
            nzeros = i - iold;
            temp = (ctx->fap2[iold] - ctx->fap2[i]) / nzeros;
            /* if the radius is too small or too big to be represented as dpe,
                * output a warning message */
            if (temp < xsmall)
            {
                rdpe_set(cfg.dradii[cfg.n_radii], small_rdpe_t);
                MPS_DEBUG(ctx, "Warning: Some zeros are too small to be\n"
                    "represented as cdpe, they are replaced by\n"
                    "small numbers and the status is set to 'F'.");
            }
            if (temp > xbig)
            {
                rdpe_set(cfg.dradii[cfg.n_radii], big_rdpe_t);
                MPS_DEBUG(ctx, "Warning: Some zeros are too big to be\n"
                    "represented as cdpe, they are replaced by\n"
                    "big numbers and the status is set to 'F'.");
            }

            /* if the radius is representable as dpe, compute it */
            if (temp <= xbig && temp >= xsmall)
            {
                rdpe_set_d(cfg.dradii[cfg.n_radii], temp);
                rdpe_exp_eq(cfg.dradii[cfg.n_radii]);
            }

            /* if the radius is greater than the radius of the cluster
                * set the radius equal to the radius of the cluster */
            if (rdpe_gt(cfg.dradii[cfg.n_radii], clust_rad))
                rdpe_set(cfg.dradii[cfg.n_radii], clust_rad);

            MPS_DEBUG_RDPE(ctx, cfg.dradii[cfg.n_radii], "cfg.dradii[%d]", cfg.n_radii);

            /* Close partitioning and start a new one */
            cfg.partitioning[++cfg.n_radii] = i;
        }
    }

    /* Set last point of the partitioning */
    cfg.partitioning[cfg.n_radii] = n;

    /* Compact radius that are too near, but not if we are shifting. */
    if (rdpe_ne(g, rdpe_zero))
        return cfg;

    for (i = 0; i < cfg.n_radii; i++)
    {
        /* Scan next radii to see if they are near the
            * i-th that we are considering now  */
        for (j = i + 1; j < cfg.n_radii; j++)
        {
            /* Get an estimate of the distance between approximations that would
                * be obtained collapsing the circles. It's a little understimated
                * but this is due to the fact that changing approximation is only
                * a fallback, and not the preferred action to perform. */
            rdpe_sub(tmp, cfg.dradii[j], cfg.dradii[i]);
            rdpe_abs_eq(tmp);
            if (rdpe_lt(cfg.dradii[i], cfg.dradii[j]))
                rdpe_div_eq(tmp, cfg.dradii[i]);
            else
                rdpe_div_eq(tmp, cfg.dradii[j]);
            rdpe_div_eq_d(tmp, PI);
            rdpe_mul_eq_d(tmp, cfg.partitioning[j + 1] - cfg.partitioning[i]);
            if (rdpe_gt(tmp, rdpe_one))
            {
                break;
            }
        }

        /* This is the number of circles that are near */
        j--;
        offset = j - i;

        /* If there is nothing to compact, do not compact it */
        if (offset == 0)
        {
            continue;
        }

        MPS_DEBUG(ctx, "Compacting circles from %d to %d", i, j);

        /* We shall now compact circles between i and j, so
            * we start computing the mean of the radius */
        rdpe_mul_eq_d(cfg.dradii[i],
            cfg.partitioning[i + 1] - cfg.partitioning[i]);
        for (k = i + 1; k <= j; k++)
        {
            rdpe_mul_d(tmp, cfg.dradii[j],
                cfg.partitioning[k + 1] - cfg.partitioning[k]);
            rdpe_add_eq(cfg.dradii[i], tmp);
        }

        rdpe_div_eq_d(cfg.dradii[i],
            cfg.partitioning[j + 1] - cfg.partitioning[i]);

        /* Move other circles backward */
        for (k = j + 1; k < cfg.n_radii; k++)
        {
            rdpe_set(cfg.dradii[k - offset], cfg.dradii[k]);
            cfg.partitioning[k - offset] = cfg.partitioning[k];
        }

        /* Set new cfg.n_radii and new partitioning */
        cfg.n_radii = cfg.n_radii - offset;
        cfg.partitioning[cfg.n_radii] = n;
    }

    mps_free(h);
    return cfg;
}

/**
    * @brief Multiprecision version of mps_fstart()
    * @see mps_fstart()
    */
MPS_PRIVATE void
    mps_mstart(mps_context * ctx, int n, mps_cluster_item * cluster_item,
        rdpe_t clust_rad,
        rdpe_t g, rdpe_t dap[], mpc_t gg)
{
    mps_cluster* cluster = NULL;
    int i, j, jj, iold, l, nzeros;
    double sigma, ang, th;
    rdpe_t big, small_rdpe_t, rtmp1, rtmp2;
    cdpe_t ctmp;
    mpc_t mtmp;
    mps_boolean need_recomputing = true;
    mps_root* root = NULL;

    if (cluster_item)
        cluster = cluster_item->cluster;

    mpc_init2(mtmp, ctx->mpwp);

    rdpe_set(small_rdpe_t, RDPE_MIN);
    rdpe_set(big, RDPE_MAX);

    if (ctx->random_seed)
        sigma = drand();
    else
    {
        /* If this is the first cluster select sigma = 0. In the other
            * case try to maximize starting points distance. */
        if (!cluster_item || cluster_item->prev_cluster_item == NULL)
        {
            sigma = ctx->last_sigma = MPS_STARTING_SIGMA;
        }
        else
        {
            sigma = mps_maximize_distance(ctx, ctx->last_sigma, cluster_item, n);
        }
    }

    nzeros = 0;

    mps_starting_configuration cfg;
    cfg.n_radii = 0;
    cfg.fradii = NULL;
    cfg.dradii = NULL;
    cfg.partitioning = NULL;

    /* Continue to cycle while clusters are _real clusters_,
        * because until now we have detached roots that were not
        * certainly outside of the cluster. */
    while (need_recomputing)
    {
        /* In the general case apply the Rouche-based criterion */
        mps_starting_configuration_clear(ctx, &cfg);
        cfg = mps_mcompute_starting_radii(ctx, n, cluster_item, clust_rad, g, dap);

        /* We need to check that the points that we have kept out of
            * the cluster are really out of the clusters. */
        need_recomputing = false;

        /* Find the maximum radii */
        rdpe_set(rtmp1, rdpe_zero);
        for (i = 0; i < cfg.n_radii; i++)
        {
            if (rdpe_lt(rtmp1, cfg.dradii[i]))
            {
                rdpe_set(rtmp1, cfg.dradii[i]);
            }
        }

        mps_cluster_item* c_item;
        for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
        {
            if (c_item->detached_cluster_item && c_item->detached_cluster_item->cluster == cluster &&
                cluster != NULL)
            {
                mps_root* root = c_item->cluster->first_root;
                i = root->k;

                /* Check if the root touches the cluster */
                mpc_sub(mtmp, ctx->approx_root[i]->mvalue, gg);
                mpc_get_cdpe(ctmp, mtmp);
                cdpe_mod(rtmp2, ctmp);

                /* i is the cluster detached from i_clust, so the unique
                    * root in it is the first. */
                rdpe_sub_eq(rtmp2, ctx->approx_root[i]->drad);
                rdpe_sub_eq(rtmp2, rtmp1);

                /* If they are too near we need to recompact them */
                if (rdpe_lt(rtmp2, rdpe_zero))
                {
                    mps_cluster_item* next_cluster_item = c_item->next_cluster_item;

                    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                        MPS_DEBUG(ctx, "Recompacting cluster with root %d", i);

                    need_recomputing = true;
                    mps_clusterization_remove_cluster(ctx, ctx->clusterization, c_item);
                    mps_cluster_insert_root(ctx, cluster, i);

                    c_item = next_cluster_item;
                }
            }
        }
    }

    th = pi2 / n;

    /* Set initial approximations accordingly to the computed
        * circles  */
    root = (cluster) ? cluster->first_root : ctx->clusterization->first_cluster_item->cluster->first_root;
    for (i = 0; i < cfg.n_radii; i++)
    {
        mps_root* starting_root = root;

        nzeros = cfg.partitioning[i + 1] - cfg.partitioning[i];
        ang = pi2 / nzeros;
        iold = cfg.partitioning[i];

        /* Compute the initial approximations */
        for (j = iold; j < cfg.partitioning[i + 1]; j++)
        {
            jj = j - iold;

            /* Take index relative to the cluster
                * that we are analyzing. */
                /* l = ctx->clust[ctx->punt[i_clust] + j]; */
            l = root->k;

            cdpe_set_d(ctmp,
                cos(ang * jj + th * cfg.partitioning[i + 1] + sigma),
                sin(ang * jj + th * cfg.partitioning[i + 1] + sigma));
            cdpe_mul_eq_e(ctmp, cfg.dradii[i]);
            cdpe_set(ctx->approx_root[l]->dvalue, ctmp);

            if (rdpe_eq(cfg.dradii[i], big) || rdpe_eq(cfg.dradii[i], small_rdpe_t))
            {
                ctx->approx_root[l]->status = MPS_ROOT_STATUS_NOT_DPE;
            }

            root = root->next_root;
        }


        /* If the new radius of the cluster is relatively small, then
            * set the status component equal to 'o' (output)
            * and set the corresponding radius */
        rdpe_set(rtmp1, cfg.dradii[i]);
        rdpe_mul_eq_d(rtmp1, (double)nzeros);
        if (rdpe_lt(rtmp1, clust_rad) && ctx->algorithm == MPS_ALGORITHM_SECULAR_GA)
            rdpe_set(rtmp1, clust_rad);
        rdpe_set(rtmp2, g);
        rdpe_mul_eq(rtmp2, ctx->eps_out);
        MPS_DEBUG(ctx, "Performing relatively small check");
        if (rdpe_le(rtmp1, rtmp2))
        {
            mps_root* root2 = starting_root;
            for (j = cfg.partitioning[i]; j < cfg.partitioning[i + 1]; j++)
            {
                l = root2->k;
                ctx->approx_root[l]->status = MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER;
                rdpe_mul_d(ctx->approx_root[l]->drad, rtmp1, (double)nzeros);
                root2 = root2->next_root;
            }
        }
        rdpe_set(clust_rad, cfg.dradii[i]);
    }

    mpc_clear(mtmp);
    mps_starting_configuration_clear(ctx, &cfg);
}

/**
* @brief This function scans the existing clusters and selects the
* ones where shift in the gravity center must be done.
* Then computes the gravity center g, performs the shift of
* the variable and compute new starting approximations in the
* cluster (floating point version).
*
* @param ctx A pointer to the current mps_context.
*
* Shift in g is perfomed if the approximation is included in the search set
* or its inclusion status has not been determined yet.
*
* To compute g, first compute the weighted mean (super center sc)
* of the approximations in the cluster, where the weight are the
* radii, then compute the radius (super radius sr) of the disk
* centered in the super center containing all the disks of the cluster.
* Apply few steps of Newton's iteration to the (m-1)-st derivative
* of the polynomial starting from the super center and obtain the
* point g where to shift the variable.
* If g is outside the super disk of center sc and radius sr
* output a warning message.
*/
MPS_PRIVATE void
    mps_frestart(mps_context* ctx)
{
    int k, j, l;
    double sr, sum, rtmp, rtmp1;
    cplx_t sc, corr, ctmp;
    mps_boolean tst;
    mps_monomial_poly* p = MPS_MONOMIAL_POLY(ctx->active_poly);
    mps_approximation* g = NULL;

    ctx->operation = MPS_OPERATION_SHIFT;

    /* Variables for cluster analysis */
    mps_cluster* cluster;
    mps_cluster_item* c_item;
    mps_root* root;

    /* For user's polynomials skip the restart stage (not yet implemented) */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        return;

    /* scan the existing clusters and  select the ones where shift in
    * the gravity center must be done. tst=true means do not perform shift */
    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
    {                           /* loop1: */
        cluster = c_item->cluster;

        /* Continue if the root is isolated */
        if (cluster->n == 1)
            continue;

        tst = true;
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {                       /* looptst : */
            l = root->k;
            if (!ctx->approx_root[l]->again)
                goto loop1;
            if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
            {
                if (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                {
                    tst = false;
                    break;
                }
            }
            else if ((ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                || (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_IN))
            {
                tst = false;
                break;
            }
        }                       /* for */
        if (tst)
            goto loop1;

        /* Compute super center sc and super radius sr */
        mps_fsrad(ctx, cluster, sc, &sr);

        /* Check the relative width of the cluster
        * If it is greater than 1 do not shift
        * and set status(:1)='c' that means
        * keep iterating Aberth's step. */
        if (sr > cplx_mod(sc))
        {
            mps_root* r;
            for (r = cluster->first_root; r != NULL; r = r->next_root)
                ctx->approx_root[r->k]->status = MPS_ROOT_STATUS_CLUSTERED;
            MPS_DEBUG(ctx, "Cluster rel. large: skip to the next component");
            goto loop1;
        }

        /* Now check the Newton isolation of the cluster */
        mps_cluster_item* c_item2;
        mps_cluster* cluster2;
        int m;
        for (c_item2 = ctx->clusterization->first_cluster_item; c_item2 != NULL; c_item2 = c_item2->next_cluster_item)
            if (c_item2 != c_item)
            {
                cluster2 = c_item2->cluster;
                for (root = cluster2->first_root; root != NULL; root = root->next_root)
                {
                    m = root->k;
                    cplx_sub(ctmp, sc, ctx->approx_root[m]->fvalue);
                    rtmp = cplx_mod(ctmp);
                    rtmp1 = (sr + ctx->approx_root[m]->frad) * 5 * ctx->n;
                    if (rtmp < rtmp1)
                    {
                        for (root = cluster->first_root; root != NULL; root = root->next_root)
                            ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
                        MPS_DEBUG_WITH_INFO(ctx, "Cluster not Newton isolated: skip "
                            "to the next component.");
                        goto loop1;
                    }
                }
            }
        /* Compute the coefficients of the derivative of p(x) having order
        * equal to the multiplicity of the cluster -1. */
        sum = 0.0;
        for (j = 0; j <= MPS_POLYNOMIAL(p)->degree; j++)
        {
            sum += cplx_mod(p->fpc[j]);
            cplx_set(p->fppc[j], p->fpc[j]);
        }
        for (j = 1; j < cluster->n; j++)
        {
            for (k = 0; k <= MPS_POLYNOMIAL(p)->degree - j; k++)
                cplx_mul_d(p->fppc[k], p->fppc[k + 1], (double)(k + 1));
        }
        for (j = 0; j < MPS_POLYNOMIAL(p)->degree - cluster->n + 2; j++)
            ctx->fap1[j] = cplx_mod(p->fppc[j]);

        /* Apply at most max_newt_it steps of Newton's iterations
        * to the above derivative starting from the super center
        * of the cluster. */

        g = mps_approximation_new(ctx);
        cplx_set(g->fvalue, sc);

        for (j = 0; j < ctx->max_newt_it; j++)
        {
            mps_fnewton(ctx, MPS_POLYNOMIAL(p), g, corr);
            cplx_sub_eq(g->fvalue, corr);
            if (!g->again)
            {
                break;
            }
        }

        if (j == ctx->max_newt_it)
        {
            MPS_DEBUG(ctx, "Exceeded maximum Newton iterations in frestart");
            mps_approximation_free(ctx, g);
            return;
        }
        cplx_sub(ctmp, sc, g->fvalue);
        if (cplx_mod(ctmp) > sr)
        {
            MPS_DEBUG(ctx, "The gravity center falls outside the cluster");
            mps_approximation_free(ctx, g);
            return;
        }
        /* Compute the coefficients of the shifted polynomial p(x+g)
        * and compute new starting approximations
        * First check if shift may cause overflow, in this case skip
        * the shift stage */

        if (ctx->n * ::log(cplx_mod(g->fvalue)) + ::log(sum) > ::log(DBL_MAX))
            goto loop1;
        MPS_DEBUG_CALL(ctx, "mps_fshift");
        mps_fshift(ctx, cluster->n, c_item, sr, g->fvalue, ctx->eps_out);
        rtmp = cplx_mod(g->fvalue);
        rtmp *= DBL_EPSILON * 2;
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;
            /* Choose as new incl. radius 2*multiplicity*(radius of the circle) */
            ctx->approx_root[l]->frad = 2 * cluster->n * cplx_mod(ctx->approx_root[l]->fvalue);
            cplx_add_eq(ctx->approx_root[l]->fvalue, g->fvalue);
            if (ctx->approx_root[l]->frad < rtmp)        /* DARIO* aggiunto 1/5/97 */
                ctx->approx_root[l]->frad = rtmp;
        }

    loop1:
        if (g)
        {
            mps_approximation_free(ctx, g);
            g = NULL;
        }
        ;
    }
}

/**
* @brief This function scans the existing clusters and selects the
* ones where shift in the gravity center must be done.
* Then computes the gravity center g, performs the shift of
* the variable and compute new starting approximations in the
* cluster (DPE version).
*
* @param ctx A pointer to the current mps_context.
*
* Shift in g is perfomed if the approximation is included in the search set
* or its inclusion status has not been determined yet.
*
* To compute g, first compute the weighted mean (super center sc)
* of the approximations in the cluster, where the weight are the
* radii, then compute the radius (super radius sr) of the disk
* centered in the super center containing all the disks of the cluster.
* Apply few steps of Newton's iteration to the (m-1)-st derivative
* of the polynomial starting from the super center and obtain the
* point g where to shift the variable.
* If g is outside the super disk of center sc and radius sr
* output a warning message.
*/
MPS_PRIVATE void
    mps_drestart(mps_context* ctx)
{
    int k, j, l;
    rdpe_t sr, rtmp, rtmp1;
    cdpe_t sc, corr, ctmp;
    mps_boolean tst;
    mps_monomial_poly* p = MPS_MONOMIAL_POLY(ctx->active_poly);

    /* Cluster related variables */
    mps_cluster_item* c_item;
    mps_cluster* cluster;
    mps_root* root;
    mps_approximation* g = NULL;

    ctx->operation = MPS_OPERATION_SHIFT;

    /*  For user's polynomials skip the restart stage (not yet implemented) */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        return;

    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
    {                           /* loop1: */
        cluster = c_item->cluster;

        if (cluster->n == 1)
            continue;

        tst = true;
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {                       /* looptst: */
            l = root->k;
            if (!ctx->approx_root[l]->again)
                goto loop1;
            if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
            {
                if (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                {
                    tst = false;
                    break;
                }
            }
            else if ((ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN) ||
                (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_IN))
            {
                tst = false;
                break;
            }
        }                       /* for */
        if (tst)
            goto loop1;

        /* Compute super center sc and super radius sr */
        mps_dsrad(ctx, cluster, sc, sr);

        /* Check the relative width of the cluster
        * If it is greater than 1 do not shift
        * and set statu(:1)='c' that means
        * keep iterating Aberth's step. */
        cdpe_mod(rtmp, sc);
        if (rdpe_gt(sr, rtmp))
        {
            for (root = cluster->first_root; root != NULL; root = root->next_root)
            {
                ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
                /* err(clust[j])=true  */
            }
            MPS_DEBUG(ctx, "cluster rel. large: skip to the next component");
            goto loop1;
        }

        /* Now check the Newton isolation of the cluster */
        mps_cluster_item* c_item2;
        mps_cluster* cluster2;

        for (c_item2 = ctx->clusterization->first_cluster_item; c_item2 != NULL; c_item2 = c_item2->next_cluster_item)
            if (c_item2 != c_item)
            {
                cluster2 = c_item2->cluster;
                for (root = cluster2->first_root; root != NULL; root = root->next_root)
                {
                    cdpe_sub(ctmp, sc, ctx->approx_root[root->k]->dvalue);
                    cdpe_mod(rtmp, ctmp);
                    rdpe_add(rtmp1, sr, ctx->approx_root[root->k]->drad);
                    rdpe_mul_eq_d(rtmp1, 2.0 * ctx->n);
                    if (rdpe_lt(rtmp, rtmp1))
                    {
                        for (root = cluster->first_root; root != NULL; root = root->next_root)
                            ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
                        MPS_DEBUG(ctx, "Cluster not Newton isolated: skip to the next component.");
                        goto loop1;
                    }
                }
            }

        /* Compute the coefficients of the derivative of p(x) having order
        * equal to the multiplicity of the cluster -1. */
        for (j = 0; j <= ctx->n; j++)
            cdpe_set(ctx->dpc2[j], p->dpc[j]);
        for (j = 1; j < cluster->n; j++)
        {
            for (k = 0; k <= ctx->n - j; k++)
                cdpe_mul_d(ctx->dpc2[k], ctx->dpc2[k + 1], (double)(k + 1));
        }
        for (j = 0; j < ctx->n - cluster->n + 2; j++)
            cdpe_mod(ctx->dap1[j], ctx->dpc2[j]);

        /* Apply at most max_newt_it steps of Newton's iterations
        * to the above derivative starting from the super center
        * of the cluster. */
        g = mps_approximation_new(ctx);
        cdpe_set(g->dvalue, sc);
        for (j = 0; j < ctx->max_newt_it; j++)
        {                       /* loop_newt: */
            mps_dnewton(ctx, MPS_POLYNOMIAL(p), g, corr);
            cdpe_sub_eq(g->dvalue, corr);
            if (!g->again)
            {
                break;
            }
        }
        if (j == ctx->max_newt_it)
        {
            MPS_DEBUG(ctx, "Exceeded maximum Newton iterations in frestart");
            mps_approximation_free(ctx, g);
            return;
        }
        cdpe_sub(ctmp, sc, g->dvalue);
        cdpe_mod(rtmp, ctmp);
        if (rdpe_gt(rtmp, sr))
        {
            MPS_DEBUG(ctx, "The gravity center falls outside the cluster");
            mps_approximation_free(ctx, g);
            return;
        }
        /* Shift the variable and compute new approximations */
        MPS_DEBUG_CALL(ctx, "mps_dshift");
        mps_dshift(ctx, cluster->n, c_item, sr, g->dvalue, ctx->eps_out);
        cdpe_mod(rtmp, g->dvalue);
        rdpe_mul_eq_d(rtmp, DBL_EPSILON * 2);
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;

            /* Choose as new incl. radius 2*multiplicity*(radius of the circle) */
            cdpe_mod(ctx->approx_root[l]->drad, ctx->approx_root[l]->dvalue);
            rdpe_mul_eq_d(ctx->approx_root[l]->drad,
                (double)(2 * cluster->n));
            cdpe_add_eq(ctx->approx_root[l]->dvalue, g->dvalue);
            if (rdpe_lt(ctx->approx_root[l]->drad, rtmp))
                rdpe_set(ctx->approx_root[l]->drad, rtmp);
        }

    loop1:
        if (g)
        {
            mps_approximation_free(ctx, g);
            g = NULL;
        }
        ;
    }
}

/**
* @brief Walk through the clusterization to find clusters that were detached
* from the given one. Check if they were detached correctly according to the
* new inclusion radii computed by mps_mrestart().
*
* @param ctx The context in which MPSolve is running.
* @param clusterization The clusterization where clusters should be looked up
* @param cluster_item The item containing the cluster from which the other clusters
* may have been detached.
*/
static void
    mps_cluster_check_detachment(mps_context* ctx, mps_clusterization* clusterization,
        mps_cluster_item* cluster_item)
{
    MPS_DEBUG_CALL(ctx, "mps_cluster_check_detachment");

    mps_cluster_item* item = NULL;
    mps_cluster* cluster = cluster_item->cluster;
    rdpe_t radius, distance_abs;
    mpc_t center, distance;

    mpc_init2(center, ctx->mpwp);
    mpc_init2(distance, ctx->mpwp);

    mps_msrad(ctx, cluster, center, radius);
    rdpe_mul_eq_d(radius, 2.0 * ctx->n);

    for (item = clusterization->first_cluster_item; item != NULL; item = item->next_cluster_item)
    {
        mps_cluster* detached_cluster = item->cluster;
        if (item->detached_cluster_item == cluster_item)
        {
            int k = detached_cluster->first_root->k;
            mpc_sub(distance, center, ctx->approx_root[k]->mvalue);
            mpc_rmod(distance_abs, distance);

            if (rdpe_lt(distance_abs, radius))
            {
                mps_cluster_item* next_cluster_item = item->next_cluster_item;
                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                    MPS_DEBUG(ctx,
                        "Cluster containing root %d has not been correctly detached, reattaching.", k);

                mps_cluster_insert_root(ctx, cluster, k);
                mps_clusterization_remove_cluster(ctx, clusterization, item);

                item = next_cluster_item;
            }
            else
            {
                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                    MPS_DEBUG(ctx,
                        "Cluster containing root %d was successfully detached.", k);

                /* We need to stop marking this cluster as detached, that means
                * "experimental" in this context. */
                item->detached_cluster_item = NULL;
            }
        }
    }

    mpc_clear(center);
    mpc_clear(distance);
}

static void
    mps_cluster_reattach_all_detached_clusters(mps_context* ctx, mps_clusterization* clusterization,
        mps_cluster_item* cluster_item)
{
    MPS_DEBUG_CALL(ctx, "mps_cluster_reattach_all_detached_clusters");

    mps_cluster_item* item = NULL;
    mps_cluster* cluster = cluster_item->cluster;

    /* Find cluster that have been detached from this and
    * attach them. */
    for (item = clusterization->first_cluster_item; item != NULL; )
    {
        /* Pre-cache the next_cluster_item cluster item se we don't ruin the
        * loop if we need to free the current one. */
        mps_cluster* d_cluster = item->cluster;
        mps_cluster_item* next_cluster_item = item->next_cluster_item;

        if (item->detached_cluster_item == cluster_item)
        {
            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                MPS_DEBUG(ctx,
                    "Reattaching root %ld to its original cluster", d_cluster->first_root->k);

            mps_cluster_insert_root(ctx, cluster, d_cluster->first_root->k);
            mps_clusterization_remove_cluster(ctx, clusterization, item);
        }

        item = next_cluster_item;
    }
}

/**
* @brief This function scans the existing clusters and selects the
* ones where shift in the gravity center must be done.
* Then computes the gravity center g, performs the shift of
* the variable and compute new starting approximations in the
* cluster (MP version).
*
* @param ctx A pointer to the current mps_context.
*
* Shift in g is perfomed if the approximation is included in the search set
* or its inclusion status has not been determined yet.
*
* To compute g, first compute the weighted mean (super center sc)
* of the approximations in the cluster, where the weight are the
* radii, then compute the radius (super radius sr) of the disk
* centered in the super center containing all the disks of the cluster.
* Apply few steps of Newton's iteration to the (m-1)-st derivative
* of the polynomial starting from the super center and obtain the
* point g where to shift the variable.
* If g is outside the super disk of center sc and radius sr
* output a warning message.
*/
MPS_PRIVATE void
    mps_mrestart(mps_context* ctx)
{
    mps_boolean tst;
    int j, k, l;
    rdpe_t sr, rtmp, rtmp1, rtmp2, ag;
    cdpe_t tmp;
    mpf_t rea, srmp;
    mpc_t sc, corr, temp;
    mps_monomial_poly* p = MPS_MONOMIAL_POLY(ctx->active_poly);

    /* Variables for cluster iteration */
    mps_cluster_item* c_item;
    mps_cluster* cluster;
    mps_root* root;
    mps_approximation* g = NULL;

    long int starting_wp = ctx->mpwp;

    MPS_DEBUG_THIS_CALL(ctx);

    ctx->operation = MPS_OPERATION_SHIFT;

    /* For user's polynomials skip the restart stage (not yet implemented) */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
    {
        MPS_DEBUG_WITH_INFO(ctx, "Skipping restart on user polynomial");
        return;
    }

    mpf_init2(rea, ctx->mpwp);
    mpf_init2(srmp, ctx->mpwp);
    mpc_init2(sc, ctx->mpwp);
    mpc_init2(corr, ctx->mpwp);
    mpc_init2(temp, ctx->mpwp);

    /* Try to detach quasi-convergent elements from the clusters */
    mps_clusterization_detach_clusters(ctx, ctx->clusterization);

    k = 0;
    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
        k = MAX(k, c_item->cluster->n);

    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
    {
        /* loop1: */
        cluster = c_item->cluster;

        if (cluster->n == 1)
            continue;

        tst = true;
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {                       /* looptst: */
            l = root->k;

            if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
            {
                if (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                {
                    tst = false;
                    break;
                }
            }
            else if ((ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN) ||
                (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_IN))
            {
                tst = false;
                break;
            }
        }                       /* for */

        if (tst)
            goto clean_detached_cluster;

        /* Compute super center sc and super radius sr */
        mps_msrad(ctx, cluster, sc, sr);

        /* MPS_DEBUG (ctx, "Clust = %d", i); */
        MPS_DEBUG_MPC(ctx, 10, sc, "Super center");
        MPS_DEBUG_RDPE(ctx, sr, "Super radius");

        /* Check the relative width of the cluster
        * If it is greater than 1 do not shift
        * and set status[:1)='c' that means
        * keep iterating Aberth's step.
        * Check also the Newton-isolation of the cluster */

        mpc_get_cdpe(tmp, sc);
        cdpe_mod(rtmp, tmp);

        if (ctx->DOLOG)
        {
            rdpe_div(rtmp2, sr, rtmp);
        }
        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
        {
            MPS_DEBUG_RDPE(ctx, rtmp2, "Relative width");
        }

        if (rdpe_gt(sr, rtmp))
        {
            for (root = cluster->first_root; root != NULL; root = root->next_root)
                ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
            MPS_DEBUG(ctx, "Cluster relat. large: skip to the next component");

            goto clean_detached_cluster;
        }

        /* Now check the Newton isolation of the cluster */
        mps_cluster_item* c_item2;
        mps_cluster* cluster2;

        rdpe_set(rtmp2, rdpe_zero);
        for (c_item2 = ctx->clusterization->first_cluster_item; c_item2 != NULL; c_item2 = c_item2->next_cluster_item)
        {
            if (c_item2 != c_item)
            {
                cluster2 = c_item2->cluster;
                for (root = cluster2->first_root; root != NULL; root = root->next_root)
                {
                    mpc_sub(temp, sc, ctx->approx_root[root->k]->mvalue);
                    mpc_get_cdpe(tmp, temp);
                    cdpe_mod(rtmp, tmp);
                    rdpe_sub_eq(rtmp, ctx->approx_root[root->k]->drad);
                    rdpe_sub_eq(rtmp, sr);
                    rdpe_inv_eq(rtmp);
                    rdpe_add_eq(rtmp2, rtmp);
                }
            }
        }
        rdpe_mul_eq(rtmp2, sr);
        rdpe_set_d(rtmp1, 0.3);

        if (rdpe_gt(rtmp2, rtmp1))
        {
            for (root = cluster->first_root; root != NULL; root = root->next_root)
                ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
            MPS_DEBUG(ctx, "Cluster not Newton isolated: skip to the next component");
            goto clean_detached_cluster;
        }

        mps_debug_cluster_structure(ctx);

        /* Create the approximation used to find the derivative zero */
        g = mps_approximation_new(ctx);
        mpc_set(g->mvalue, sc);

        /* Compute the coefficients of the derivative of p(x) having order
        * equal to the multiplicity of the cluster -1. */
        mps_monomial_poly* der = mps_monomial_poly_derive(ctx, p, cluster->n - 1, ctx->mpwp);
        rdpe_t epsilon;
        rdpe_t error;

        rdpe_set(epsilon, rdpe_zero);

        /* Evaluate the necessary precision to perform the Newton iterations. */
        do
        {
            long int prec = mps_monomial_poly_get_precision(ctx, der);
            mps_mhorner_with_error2(ctx, der, g->mvalue, corr, error, prec);

            mpc_rmod(rtmp, corr);
            rdpe_set_2dl(epsilon, 1.0, -ctx->mpwp);
            rdpe_mul_eq(epsilon, rtmp);

            if (rdpe_gt(error, epsilon))
            {
                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                    MPS_DEBUG(ctx, "Raising precision of the derivative to %ld bits", prec * 2);
                mps_polynomial_raise_data(ctx, MPS_POLYNOMIAL(der), prec * 2);
                mpc_set_prec(g->mvalue, prec * 2);
            }
        } while (rdpe_gt(error, epsilon) && mps_monomial_poly_get_precision(ctx, der) <= cluster->n * ctx->mpwp);

        /*  Apply at most max_newt_it steps of Newton's iterations */
        /*  to the above derivative starting from the super center */
        /*  of the cluster. */

        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
        {
            MPS_DEBUG_MPC(ctx, 30, g->mvalue, "g before newton");
        }

        for (j = 0; j < ctx->max_newt_it; j++)
        {                       /* loop_newt: */
            mps_mnewton(ctx, MPS_POLYNOMIAL(der), g, corr,
                mpc_get_prec(g->mvalue));
            mpc_sub_eq(g->mvalue, corr);

            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            {
                MPS_DEBUG_RDPE(ctx, g->drad, "Radius of the cluster");
                MPS_DEBUG_MPC(ctx, 100, g->mvalue, "Iteration %d on the derivative", j);
            }

            mpc_rmod(ag, g->mvalue);

            if (!g->again || rdpe_Esp(g->drad) < rdpe_Esp(ag) - ctx->mpwp)
                break;
        }

        mps_polynomial_free(ctx, MPS_POLYNOMIAL(der));

        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            MPS_DEBUG(ctx, "Performed %d Newton iterations", j + 1);
        if (j == ctx->max_newt_it)
        {
            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                MPS_DEBUG(ctx, "Exceeded maximum number of Newton iterations.");

            goto clean_detached_cluster;
        }

        mpc_sub(temp, sc, g->mvalue);
        mpc_get_cdpe(tmp, temp);
        cdpe_mod(rtmp, tmp);
        if (rdpe_gt(rtmp, sr))
        {
            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                MPS_DEBUG(ctx, "The gravity center falls outside the cluster");
            goto clean_detached_cluster;
        }

        /* shift the variable and compute new approximations */
        MPS_DEBUG_CALL(ctx, "mps_mshift");
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;
            mpc_get_cdpe(ctx->approx_root[l]->dvalue, ctx->approx_root[l]->mvalue);
        }

        /*#D perform shift only if the new computed sr is smaller than old*0.25 */
        /* if (ctx->algorithm == MPS_ALGORITHM_SECULAR_GA) */
        /*        rdpe_set (rtmp, sr); */
        /* else */
        rdpe_mul_d(rtmp, sr, 0.25);

        /*#D AGO99 Factors: 0.1 (MPS2.0), 0.5 (GIUGN98) */
        mps_mshift(ctx, c_item->cluster->n, c_item, sr, g->mvalue);
        if (rdpe_le(sr, rtmp))
        {                       /* Perform shift only if the new clust is smaller */
            mpc_get_cdpe(tmp, g->mvalue);
            cdpe_mod(rtmp, tmp);
            if (ctx->lastphase == mp_phase)
                rdpe_set_2dl(rtmp1, 1.0, -starting_wp);
            else
                rdpe_set_d(rtmp1, DBL_EPSILON);
            rdpe_mul_eq(rtmp, rtmp1);
            rdpe_mul_eq_d(rtmp, 2);

            MPS_DEBUG_RDPE(ctx, rtmp, "rtmp");
            MPS_DEBUG_RDPE(ctx, rtmp1, "rtmp1");

            for (root = cluster->first_root; root != NULL; root = root->next_root)
            {
                l = root->k;
                mpc_set_cdpe(ctx->approx_root[l]->mvalue, ctx->approx_root[l]->dvalue);
                mpc_add_eq(ctx->approx_root[l]->mvalue, g->mvalue);
                cdpe_mod(rtmp1, ctx->approx_root[l]->dvalue);
                rdpe_mul_d(ctx->approx_root[l]->drad, rtmp1,
                    2.0 * cluster->n);
                if (rdpe_lt(ctx->approx_root[l]->drad, rtmp))
                    rdpe_set(ctx->approx_root[l]->drad, rtmp);

                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                {
                    MPS_DEBUG_MPC(ctx, ctx->mpwp, ctx->approx_root[l]->mvalue, "Repositioned root %4d", l);
                    MPS_DEBUG_RDPE(ctx, ctx->approx_root[l]->drad, "Radius for root %4d", l);
                }
            }

            /* Check if the clusters that have been detached from this have been
            * detached for a good reason. */
            mps_cluster_check_detachment(ctx, ctx->clusterization, c_item);
        }
        else
        {
            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                MPS_DEBUG(ctx, "DO NOT PERFORM RESTART, "
                    "new radius of the cluster is larger");

            goto clean_detached_cluster;
        }

    clean_detached_cluster:
        mps_cluster_reattach_all_detached_clusters(ctx, ctx->clusterization, c_item);

        if (g != NULL)
        {
            mps_approximation_free(ctx, g);
            g = NULL;
        }
    }

    mpc_clear(temp);
    mpc_clear(corr);
    mpc_clear(sc);
    mpf_clear(srmp);
    mpf_clear(rea);
}

/**
* @brief Check if the clusters are Newton isolated, in a way that we can
* apply the shift in the gravity center without problems deriving from
* other roots.
*
* @param ctx A pointer to the current mps_context.
*/
MPS_PRIVATE void
    mps_mnewtis(mps_context* ctx)
{
    mps_boolean tst;
    int k, l;
    rdpe_t sr, rtmp, rtmp1;
    cdpe_t tmp;
    mpf_t rea, srmp;
    mpc_t sc, temp;
    rdpe_t rtmp2;

    /* Cluster variables */
    mps_cluster_item* c_item;
    mps_cluster* cluster;
    mps_root* root;

    /* For user's polynomials skip the restart stage (not yet implemented) */
    if (!MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        return;
    mpf_init2(rea, ctx->mpwp);
    mpf_init2(srmp, ctx->mpwp);
    mpc_init2(sc, ctx->mpwp);
    mpc_init2(temp, ctx->mpwp);

    k = 0;
    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
        k = MAX(k, c_item->cluster->n);

    for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
    {                           /* loop1: */
        cluster = c_item->cluster;

        if (cluster->n == 1)
            continue;

        tst = true;

        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {                       /* looptst: */
            l = root->k;
            if (!ctx->approx_root[l]->again)
                goto loop1;
            if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
            {
                if (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                {
                    tst = false;
                    break;
                }
            }
            else if ((ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN) ||
                (ctx->approx_root[l]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[l]->inclusion == MPS_ROOT_INCLUSION_IN))
            {
                tst = false;
                break;
            }
        }                       /* for */
        if (tst)
            goto loop1;

        /* Compute super center sc and super radius sr */
        mpf_set_ui(srmp, 0);
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;
            mpf_set_rdpe(rea, ctx->approx_root[l]->drad);
            mpf_add(srmp, srmp, rea);
        }
        mpc_set_ui(sc, 0, 0);
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;
            mpf_set_rdpe(rea, ctx->approx_root[l]->drad);
            mpc_mul_f(temp, ctx->approx_root[l]->mvalue, rea);
            mpc_add_eq(sc, temp);
        }
        mpc_div_eq_f(sc, srmp);
        rdpe_set(sr, rdpe_zero);
        for (root = cluster->first_root; root != NULL; root = root->next_root)
        {
            l = root->k;
            mpc_sub(temp, sc, ctx->approx_root[l]->mvalue);
            mpc_get_cdpe(tmp, temp);
            cdpe_mod(rtmp, tmp);
            rdpe_add_eq(rtmp, ctx->approx_root[l]->drad);
            if (rdpe_lt(sr, rtmp))
                rdpe_set(sr, rtmp);
        }

        /* Check the relative width of the cluster
        * If it is greater than 1 do not shift
        * and set status[:1)='c' that means
        * keep iterating Aberth's step.
        * Check also the Newton-isolation of the cluster */
        mpc_get_cdpe(tmp, sc);
        cdpe_mod(rtmp, tmp);
        rdpe_div(rtmp2, sr, rtmp);
        if (rdpe_gt(sr, rtmp))
        {
            for (root = cluster->first_root; root != NULL; root = root->next_root)
                ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
            MPS_DEBUG(ctx, "Custer relatively large: "
                "skip to the next compoent");
            goto loop1;
        }

        /* Now check the Newton isolation of the cluster */
        mps_cluster_item* c_item2;
        mps_cluster* cluster2;
        rdpe_set(rtmp2, rdpe_zero);

        for (c_item2 = ctx->clusterization->first_cluster_item; c_item2 != NULL; c_item2 = c_item2->next_cluster_item)
        {
            if (c_item2 != c_item)
            {
                cluster2 = c_item2->cluster;
                for (root = cluster2->first_root; root != NULL; root = root->next_root)
                {
                    mpc_sub(temp, sc, ctx->approx_root[root->k]->mvalue);
                    mpc_get_cdpe(tmp, temp);
                    cdpe_mod(rtmp, tmp);
                    rdpe_sub_eq(rtmp, ctx->approx_root[root->k]->drad);
                    rdpe_sub_eq(rtmp, sr);
                    rdpe_inv_eq(rtmp);
                    rdpe_add_eq(rtmp2, rtmp);
                }
            }
        }
        rdpe_mul_eq(rtmp2, sr);
        rdpe_set_d(rtmp1, 0.3);

        if (rdpe_gt(rtmp2, rtmp1))
        {
            for (root = cluster->first_root; root != NULL; root = root->next_root)
                ctx->approx_root[root->k]->status = MPS_ROOT_STATUS_CLUSTERED;
            MPS_DEBUG(ctx, "Cluster not Newton isolated: "
                "skip to the next component");
            goto loop1;
        }
        ctx->newtis = 1;

    loop1:
        ;
    }

    mpc_clear(temp);
    mpc_clear(sc);
    mpf_clear(srmp);
    mpf_clear(rea);
}
