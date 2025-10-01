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
#include <mps/mps.h>
#include <float.h>

/**
* @brief Set <code>again[i]</code> to <code>true</code> or to <code>false</code>
* according to the values of <code>status</code> and <code>inclusion</code>
*  in the mps_approximation and the current <code>goal</code>.
*
* More precisely:
*
* - If goal is <code>MPS_OUTPUT_GOAL_COUNT</code>: <code>true</code> only
*   if inclusion is <code>MPS_ROOT_INCLUSION_UNKNOWN</code> and status is not
*   <code>MPS_ROOT_STATUS_APPROXIMATED</code>, <code>MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER</code>
*   or <code>MPS_ROOT_STATUS_NOT_FLOAT</code>.
*
* - If goal is <code>MPS_OUTPUT_GOAL_ISOLATE</code> or <code>MPS_OUTPUT_GOAL_APPROXIMATE</code>:
*   <code>true</code> if <code>status</code> is <code>MPS_ROOT_STATUS_CLUSTERED</code> or <code>inclusion</code>
*   is <code>MPS_ROOT_INCLUSION_UNKNOWN</code>.
*/
MPS_PRIVATE void
    mps_update(mps_context* ctx)
{
    int i;

    for (i = 0; i < ctx->n; i++)
        ctx->approx_root[i]->again = false;
    switch (ctx->output_config->goal)
    {
    case MPS_OUTPUT_GOAL_COUNT:                  /*  count */
        for (i = 0; i < ctx->n; i++)
        {
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                if (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                    ctx->approx_root[i]->status != MPS_ROOT_STATUS_NOT_DPE)
                    ctx->approx_root[i]->again = true;

            if (ctx->output_config->multiplicity &&
                ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
                ctx->approx_root[i]->again = true;

            if (ctx->output_config->root_properties &&
                (ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                    (ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_UNKNOWN ||
                        (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                            ctx->approx_root[i]->status != MPS_ROOT_STATUS_NOT_DPE))))
                ctx->approx_root[i]->again = true;
        }
        break;

    case MPS_OUTPUT_GOAL_ISOLATE:                  /* isolate */
        for (i = 0; i < ctx->n; i++)
        {
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN ||
                (ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_IN))
                if (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                    (ctx->approx_root[i]->status != MPS_ROOT_STATUS_ISOLATED ||
                        ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_IN))
                    ctx->approx_root[i]->again = true;

            if (ctx->output_config->multiplicity &&
                ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
                ctx->approx_root[i]->again = true;

            if (ctx->output_config->root_properties)
            {
                if (ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                    (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                        ctx->approx_root[i]->status != MPS_ROOT_STATUS_NOT_DPE))
                    ctx->approx_root[i]->again = true;
            }
        }

        break;

    case MPS_OUTPUT_GOAL_APPROXIMATE:                  /* approximate (the same as isolate) */
        for (i = 0; i < ctx->n; i++)
        {
            if (ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN ||
                (ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                    ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_IN))
                if (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status))
                    ctx->approx_root[i]->again = true;

            if (ctx->output_config->multiplicity &&
                ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
                ctx->approx_root[i]->again = true;

            if (ctx->output_config->root_properties)
            {
                if (ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                    MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status))
                    ctx->approx_root[i]->again = true;
            }
        }
        break;
    }
}

/**
* @brief Compute super center and super radius
*
* This routines the super radius of the <code>i</code>-th cluster,
* i.e. the radius of the inclusion disc for the whole cluster
*
* @param ctx The <code>mps_context</code> associated with the current computation.
* @param cluster cluster of which the super center and radius should be computed.
* @param sc Center of the cluster;
* @param sr Double that will be set to the super radius of the cluster;
*/
MPS_PRIVATE void
    mps_fsrad(mps_context* ctx, mps_cluster* cluster, cplx_t sc, double* sr)
{
    cplx_t ctmp;
    double sum;
    int l;

    mps_root* root;

    sum = 0.0;
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        sum += ctx->approx_root[l]->frad;
    }
    cplx_set(sc, cplx_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        cplx_mul_d(ctmp, ctx->approx_root[l]->fvalue, ctx->approx_root[l]->frad);
        cplx_add_eq(sc, ctmp);
    }
    cplx_div_eq_d(sc, sum);
    *sr = 0.0;
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        cplx_sub(ctmp, sc, ctx->approx_root[l]->fvalue);
        *sr = MAX(*sr, ctx->approx_root[l]->frad + cplx_mod(ctmp));
    }
}

/**
* @brief <code>dpe</code> version of <code>fsrad()</code>
*/
MPS_PRIVATE void
    mps_dsrad(mps_context* ctx, mps_cluster* cluster, cdpe_t sc, rdpe_t sr)
{
    cdpe_t ctmp;
    rdpe_t sum, rtmp;
    int l;
    mps_root* root;

    rdpe_set(sum, rdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        rdpe_add_eq(sum, ctx->approx_root[l]->drad);
    }
    cdpe_set(sc, cdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        cdpe_mul_e(ctmp, ctx->approx_root[l]->dvalue, ctx->approx_root[l]->drad);
        cdpe_add_eq(sc, ctmp);
    }
    cdpe_div_eq_e(sc, sum);
    rdpe_set(sr, rdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        cdpe_sub(ctmp, sc, ctx->approx_root[l]->dvalue);
        cdpe_mod(rtmp, ctmp);
        rdpe_add_eq(rtmp, ctx->approx_root[l]->drad);
        if (rdpe_lt(sr, rtmp))
            rdpe_set(sr, rtmp);
    }
}

/**
* @brief Multiprecision versione of <code>fsrad()</code>
*/
MPS_PRIVATE void
    mps_msrad(mps_context* ctx, mps_cluster* cluster, mpc_t sc, rdpe_t sr)
{
    int l;
    rdpe_t rtmp;
    cdpe_t cdtmp;
    mpf_t ftmp, sum;
    mpc_t ctmp;
    mps_root* root;

    mpc_init2(ctmp, ctx->mpwp);
    mpf_init2(ftmp, ctx->mpwp);
    mpf_init2(sum, ctx->mpwp);

    mpf_set_ui(sum, 0);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        mpf_set_rdpe(ftmp, ctx->approx_root[l]->drad);
        mpf_add(sum, sum, ftmp);
    }

    mpc_set_ui(sc, 0, 0);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        mpf_set_rdpe(ftmp, ctx->approx_root[l]->drad);
        mpc_mul_f(ctmp, ctx->approx_root[l]->mvalue, ftmp);
        mpc_add_eq(sc, ctmp);
    }

    mpc_div_eq_f(sc, sum);
    rdpe_set(sr, rdpe_zero);
    for (root = cluster->first_root; root != NULL; root = root->next_root)
    {
        l = root->k;
        mpc_sub(ctmp, sc, ctx->approx_root[l]->mvalue);
        mpc_get_cdpe(cdtmp, ctmp);
        cdpe_mod(rtmp, cdtmp);
        rdpe_add_eq(rtmp, ctx->approx_root[l]->drad);
        if (rdpe_lt(sr, rtmp))
            rdpe_set(sr, rtmp);
        else
        {
            if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            {
                MPS_DEBUG_RDPE(ctx, sr, "Super radius of the cluster");
            }
        }
    }

    mpf_clear(sum);
    mpf_clear(ftmp);
    mpc_clear(ctmp);
}



/**
* @brief Check if the roots are computed with the required
* precision.
*
* Set <code>computed</code> to <code>true</code>
* if the stop condition is satisfied,
* otherwise set <code>computed</code> to <code>false</code>.
*
* The stop condition is obtained from the vector <code>status</code>
* as follows:
*
* If the <code>goal</code> is count stop if
*  -  <code>**u</code> does not exist, except for
*     <code>a*u</code>, <code>o*u</code>, <code>f*u</code>
*  - Mult. and does not exist <code>c**</code>
*  - Real. and does not exist <code>*u*</code>, except for <code>au*</code>,
*    <code>ou*</code>;
*  - Imag  and does not exist <code>*v*</code>, except for
*    <code>av*</code>, <code>ov*</code>;
*
* If the <code>goal</code> is isolate or approximate stop if:
* - <code>**u</code> does not exist, except for  <code>a*u</code>,
*   <code>o*u</code>, <code>f*u</code>
*   and if <code>c*i</code> does not exist;
* - Mult. and does not exist <code>c*i</code>, <code>o*i</code>
* - Real. and does not exist <code>*ui</code>, except for <code>aui</code>,
*   <code>oui</code>;
* - Imag  and does not exist <code>*vi</code>,
*   except for <code>avi</code>, <code>ovi</code>;
*
* @see status
*/
MPS_PRIVATE mps_boolean
    mps_check_stop(mps_context* ctx)
{
    MPS_DEBUG_THIS_CALL(ctx);

    int i;
    mps_boolean computed;

    computed = false;
    /* count */
    if (ctx->output_config->goal == MPS_OUTPUT_GOAL_COUNT)
    {
        for (i = 0; i < ctx->n; i++)
        {
            if (!MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN)
                return computed;

            if (ctx->output_config->multiplicity &&
                ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
                return computed;

            if (ctx->output_config->root_properties &&
                ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
                !MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                ctx->approx_root[i]->status != MPS_ROOT_STATUS_MULTIPLE)
                return computed;
        }
        computed = true;
    }

    /* isolate or approximate */
    if (ctx->output_config->goal == MPS_OUTPUT_GOAL_ISOLATE ||
        ctx->output_config->goal == MPS_OUTPUT_GOAL_APPROXIMATE)
    {
        for (i = 0; i < ctx->n; i++)
        {
            if ((ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_UNKNOWN || ctx->approx_root[i]->inclusion == MPS_ROOT_INCLUSION_IN) &&
                !MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
            {
                if (ctx->debug_level & MPS_DEBUG_PACKETS)
                    MPS_DEBUG(ctx, "Cannot stop because root %d is not approximated, and its inclusion is not certain", i);
                return computed;
            }

            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT)
            {
                if (ctx->debug_level & MPS_DEBUG_PACKETS)
                    MPS_DEBUG(ctx, "Cannot stop because root %d is clustered and not certainly out of the target set", i);
                return computed;
            }

            if (ctx->output_config->multiplicity &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
                ctx->approx_root[i]->status == MPS_ROOT_STATUS_CLUSTERED)
            {
                if (ctx->debug_level & MPS_DEBUG_PACKETS)
                    MPS_DEBUG(ctx, "Cannot stop because root %d is not out of the target set and is clustered", i);
                return computed;
            }

            if (ctx->output_config->root_properties &&
                ctx->approx_root[i]->attrs == MPS_ROOT_ATTRS_NONE &&
                ctx->approx_root[i]->inclusion != MPS_ROOT_INCLUSION_OUT &&
                !MPS_ROOT_STATUS_IS_APPROXIMATED(ctx->approx_root[i]->status) &&
                ctx->approx_root[i]->status != MPS_ROOT_STATUS_MULTIPLE)
            {
                if (ctx->debug_level & MPS_DEBUG_PACKETS)
                    MPS_DEBUG(ctx, "Cannot stop because properties of root %d have not been detected, it's not out of the target set, nor approximated or multiple", i);
                return computed;
            }
        }

        MPS_DEBUG(ctx, "All roots are computed, stopping Aberth");
        computed = true;
    }

    return computed;
}

/**
* @brief Actually solve the polynomial
*
* This routine performs the following computations:
* -# Select starting approximations and check if some of them need dpe.
*    Initialize the vector again which is true if the corresponding
*    approximation is out of the root neighbourhood.
* -# Performs max_pack packets of Aberth iterations on all the
*    components out of the root neighbourhood belonging to the set S
*    and on which it is possible to iterate with float.
*    More precisely, each packet performs max_it iterations on all the
*    components where again is true. At each iteration check if the
*    current approximation is in the root neighbourhood; in this case
*    set 'again' to false.
* -# If at the end of the general packet all the approximations
*    are inside the root neighbourhood, i.e., 'again' is false in all
*    the components then return.
*    else, perform cluster analysis, select new starting approximations
*    update the vector 'statu', update the vector 'again' that selects
*    the components on which to iterate, according to the goal, and
*    repeat until the max number of allowed packets is reached.
*    In the latter case output FAILURE.
*
* The local variable <code>again</code> controls the iteration: i.e.,
*   <code>again[i]=true</code> means iterate on the <code>i</code>-th
*   component
*
* @param ctx The <code>mps_context</code> associated with the current computation.
* @param d_after_f this variable is <code>true</code> if dpe
* are needed after the floating point pass.
*/
MPS_PRIVATE void
    mps_fsolve(mps_context* ctx, mps_boolean* d_after_f)
{
    mps_boolean excep;
    int it_pack, iter, nit, oldnclust, i, j, required_zeros = ctx->n;
    mps_polynomial* p = ctx->active_poly;
    double* frad = double_valloc(ctx->n);

    /* == 1 ==  Initialize variables */
    it_pack = 0;
    mps_cluster_reset(ctx);
    for (i = 0; i < ctx->n; i++)
    {
        ctx->approx_root[i]->again = true;
        cplx_set(ctx->approx_root[i]->fvalue, cplx_zero);
        ctx->approx_root[i]->frad = DBL_MAX;
    }

    /* choose starting approximations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "FSOLVE: call fstart");

    mps_polynomial_fstart(ctx, p, ctx->approx_root);

    /***************
        this part of code performs shift in the gravity center of the roots
        In order to use it, uncomment the part below and comment the
        instruction above. Dangerous for overflow.
    ************/
    /*
        {
        cplx_t ft;
        cplx_mul_d(ft, fpc[n], -n);
        cplx_div(ft, fpc[n-1], ft);
        fshift(n, 0, 100, ft, eps_out);
        }
    *//* till here */

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* Check if there are too large or too small approximations */
    *d_after_f = false;
    for (i = 0; i < ctx->n; i++)
        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
        {
            ctx->approx_root[i]->again = false;
            *d_after_f = true;
        }

    /* == 2 ==  Perform max_pack packets of Aberth's iterations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   FSOLVE:  call fpolzer\n");
    for (iter = 0; iter < ctx->max_pack; iter++)
    {                           /* floop: */
    /* mps_fpolzer(ctx, &nit, &excep);   */
        mps_thread_fpolzer(ctx, &nit, &excep, required_zeros--);
        it_pack += nit;

        /* This flag will be set in case any floating point
        * exception happens during this phase. In this case, we are
        * required to switch to the DPE phase. */
        if (ctx->skip_float)
        {
            *d_after_f = true;
            ctx->skip_float = false;

            /* Fake that all the roots are not representable in floating point.
            * In this way the following routines will take care of reinitializing the
            * itearation completely. */
            for (i = 0; i < ctx->n; i++)
            {
                ctx->approx_root[i]->status = MPS_ROOT_STATUS_NOT_FLOAT;
            }

            goto fsolve_final_cleanup;
        }

        if (ctx->DOLOG)
            fprintf(ctx->logstr, "Packet %d  iterations= %d\n", iter, nit);

        /* perform cluster analysis, shift, restart, update 'statu', and
        * update 'again'
        */
        if (excep)
        {
            oldnclust = ctx->clusterization->n;

            if (ctx->DOLOG)
                fprintf(ctx->logstr, "   FSOLVE: call fcluster\n");
            /* cluster analysis */

            /* Compute the inclusion radii with Gerschgorin so we can compute
            * clusterizations for the roots. */
            mps_fradii(ctx, ctx->active_poly, frad);
            mps_fcluster(ctx, frad, 2 * ctx->n);   /* Isolation factor */

            MPS_DEBUG(ctx, "oldncluster = %d, ncluster = %ld", oldnclust, ctx->clusterization->n);

            if (oldnclust == ctx->clusterization->n)
            {
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   FSOLVE: cycle\n");
                continue;
            }
            else
            {
                /* modify the vector status and mark also the old
                * clusters with 'C'
                */
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   FSOLVE: call modify\n");
                mps_fmodify(ctx, true);

                if (iter == 0)
                    for (i = 0; i < ctx->n; i++)
                        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                            ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

                /* If the polynomial is not given in terms of its coeff. then
                * skip the restart stage */
                if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
                {
                    /* choose new starting approximations only for new clusters */
                    if (ctx->DOLOG)
                        fprintf(ctx->logstr, "   FSOLVE: call frestart\n");
                    mps_frestart(ctx);
                }
                /* reset the status vector */
                for (j = 0; j < ctx->n; j++)
                {
                    if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                        ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                    ctx->again_old[j] = ctx->approx_root[j]->again;
                }

                /* update 'again' */
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   FSOLVE: call update\n");
                mps_update(ctx);

                /* adjust 'again' This is needed since we are
                * between two packets */
                for (i = 0; i < ctx->n; i++)
                    if (!ctx->again_old[i])
                        ctx->approx_root[i]->again = false;
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   FSOLVE: call checkstop\n");
                /* Check the stop condition */
                if (mps_check_stop(ctx))
                    goto fsolve_final_cleanup;
            }
        }
        else
            break;
    }

    /* The 'floop' has been completed:
    * If the max number of iteration has been reached then output FAILURE */
    if (iter == ctx->max_pack)
    {
        mps_dump(ctx);
        mps_error(ctx,
            "Float: reached the maximum number of packet iterations");
    }
    /* Otherwise exit since all the approximations are
    * in the root neighbourhood, except for the ones that cannot be
    * represented as double. */

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "FLOAT: nit= %d\n", it_pack);

    /* Update */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   FSOLVE: call fcluster\n");
    oldnclust = ctx->clusterization->n;

    /* Compute the inclusion radii with Gerschgorin so we can compute
    * clusterizations for the roots. */
    mps_fradii(ctx, ctx->active_poly, frad);
    mps_fcluster(ctx, frad, 2 * ctx->n);   /* Isolation factor */

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   FSOLVE: call modify\n");
    mps_fmodify(ctx, true);

    /* reset the status vector */
    for (j = 0; j < ctx->n; j++)
    {
        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NOT_FLOAT)
            *d_after_f = true;
    }


fsolve_final_cleanup:
    double_vfree(frad);
}

/**
    * @brief This routine applies <code>nit</code> iterations of
    * Aberth's method.
    *
    * The method is applied to
    * the <code>i</code>-th component of the approximations for which
    * <code>again[i]</code>
    * is <code>true</code>. Set <code>again[i]=false</code> if the <code>i</code>-th
    * approximation is in
    * the root neighbourhood. Stop if <code>again[i]=false</code> for any <code>i</code>.
    *
    * @param ctx The <code>mps_context</code> associated with the current computation.
    * @param it Index of the component on which the iteration is needed.
    * @param excep This variable is set to <code>true</code> if after <code>nit</code>
    * iterations some approximation is still out of the root neighbourhood.
    */
MPS_PRIVATE void
    mps_fpolzer(mps_context* ctx, int* it, mps_boolean* excep)
{
    int i, iter, nzeros;
    cplx_t corr, abcorr;
    double rad1_double, modcorr;
    mps_polynomial* p = ctx->active_poly;

    /* initialize the iteration counter */
    *it = 0;
    *excep = false;

    /* count the number of approximations in the root neighbourhood */
    nzeros = 0;
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
        return;

    /* Start Aberth's iterations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "FPOLZER: starts aberth it\n");

    for (iter = 0; iter < ctx->max_it; iter++)
    {                           /* do_iter : DO iter=1,nit */
        if (ctx->DOLOG)
        {
            fprintf(ctx->logstr, "FPOLZER: iteration %d\n", iter);
            mps_dump(ctx);
        }

        for (i = 0; i < ctx->n; i++)
        {                       /* do_index */
            if (ctx->approx_root[i]->again)
            {
                (*it)++;
                rad1_double = ctx->approx_root[i]->frad;
                mps_polynomial_fnewton(ctx, p, ctx->approx_root[i], corr);
                if (iter == 0 && !ctx->approx_root[i]->again && ctx->approx_root[i]->frad > rad1_double && rad1_double
                    != 0)
                    ctx->approx_root[i]->frad = rad1_double;
                /***************************************
                    The above condition is needed to cope with the case
                    where at the first iteration the starting point
                    is already in the root neighbourhood and the actually
                    computed radius is too big since the value of the first
                    derivative is too small.
                    In this case the previous radius bound, obtained by
                    means of Rouche' is more reliable and strict
                    **************************************/

                if (ctx->approx_root[i]->again ||
                    /* the correction is performed only if iter!=1 or rad(i)!=rad1_double */
                    iter != 0 || ctx->approx_root[i]->frad != rad1_double)
                {
                    mps_faberth(ctx, ctx->approx_root[i], abcorr);
                    cplx_mul_eq(abcorr, corr);
                    cplx_sub(abcorr, cplx_one, abcorr);
                    cplx_div(abcorr, corr, abcorr);
                    cplx_sub_eq(ctx->approx_root[i]->fvalue, abcorr);
                    modcorr = cplx_mod(abcorr);
                    ctx->approx_root[i]->frad += modcorr;
                }

                /* check for new approximated roots */
                if (!ctx->approx_root[i]->again)
                {
                    nzeros++;
                    if (nzeros == ctx->n)
                        return;
                }
            }
        }
    }
    *excep = true;
}

/**
    * @brief <code>dpe</code> version of <code>fpolzer()</code>.
    */
MPS_PRIVATE void
    mps_dpolzer(mps_context* ctx, int* it, mps_boolean* excep)
{
    int iter, i, nzeros;
    rdpe_t rad1_rdpe_t, rtmp;
    cdpe_t corr, abcorr;
    mps_polynomial* p = ctx->active_poly;

    /* initialize the iteration counter */
    *it = 0;
    *excep = false;

    /* count the number of approximations in the root neighbourhood */
    nzeros = 0;
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
        return;

    /* Start Aberth's iterations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "DPOLZER: starts aberth\n");
    for (iter = 0; iter < ctx->max_it; iter++)
    {                           /* do_iter: */
        for (i = 0; i < ctx->n; i++)
        {                       /* do_index: */
            if (ctx->approx_root[i]->again)
            {
                (*it)++;
                rdpe_set(rad1_rdpe_t, ctx->approx_root[i]->drad);
                mps_polynomial_dnewton(ctx, p, ctx->approx_root[i], corr);
                if (iter == 0 && !ctx->approx_root[i]->again && rdpe_gt(ctx->approx_root[i]->drad, rad1_rdpe_t)
                    && rdpe_ne(rad1_rdpe_t, rdpe_zero))
                    rdpe_set(ctx->approx_root[i]->drad, rad1_rdpe_t);
                /************************************************
                    The above condition is needed to manage with the case where
                    at the first iteration the starting point is already in the
                    root neighbourhood and the actually computed radius is too
                    big since the value of the first derivative is too small.
                    In this case the previous radius bound, obtained by means of
                    Rouche' is more reliable and strict
                    **********************************************/

                if (ctx->approx_root[i]->again ||
                    /* the correction is performed only if iter!=1 or rad(i)!=rad1_rdpe_t */
                    iter != 0
                    || rdpe_ne(ctx->approx_root[i]->drad, rad1_rdpe_t))
                {
                    mps_daberth(ctx, ctx->approx_root[i], abcorr);
                    cdpe_mul_eq(abcorr, corr);
                    cdpe_sub(abcorr, cdpe_one, abcorr);
                    cdpe_div(abcorr, corr, abcorr);
                    cdpe_sub_eq(ctx->approx_root[i]->dvalue, abcorr);
                    cdpe_mod(rtmp, abcorr);
                    rdpe_add_eq(ctx->approx_root[i]->drad, rtmp);
                }

                /* check for new approximated roots */
                if (!ctx->approx_root[i]->again)
                {
                    nzeros++;
                    if (nzeros == ctx->n)
                        return;
                }
            }
        }
    }
    *excep = true;
}

/**
    * @brief <code>dpe</code> version of <code>fsolve()</code>.
    */
MPS_PRIVATE void
    mps_dsolve(mps_context* ctx, mps_boolean d_after_f)
{
    int it_pack, iter, nit, oldnclust, i, j, required_zeros = ctx->n;
    mps_boolean excep;
    mps_polynomial* p = ctx->active_poly;
    rdpe_t* drad = rdpe_valloc(ctx->n);

    if (ctx->DOLOG)
    {
        fprintf(ctx->logstr, "   DSOLVE: d_after_f= %d\n", d_after_f);
    }

    /* == 1 == Initialize variables */
    it_pack = 0;

    if (d_after_f)
        for (i = 0; i < ctx->n; i++)
            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
            {
                ctx->approx_root[i]->again = true;
                rdpe_set_d(ctx->approx_root[i]->drad, DBL_MAX);
            }
            else
                ctx->approx_root[i]->again = false;
    else
    {
        mps_cluster_reset(ctx);
        for (i = 0; i < ctx->n; i++)
        {
            ctx->approx_root[i]->again = true;
            rdpe_set_d(ctx->approx_root[i]->drad, DBL_MAX);
            cdpe_set(ctx->approx_root[i]->dvalue, cdpe_zero);
        }
    }

    /* Choose starting approximations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   DSOLVE: call dstart with again=\n");

    mps_polynomial_dstart(ctx, p, ctx->approx_root);

    /* Now adjust the status vector */
    if (d_after_f)
        for (i = 0; i < ctx->n; i++)
            if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NOT_FLOAT)
                ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
        mps_dump(ctx);

    /* == 2 == Perform ctx->max_pack  packets of Aberth's iterations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   DSOLVE: call dpolzero\n");

    for (iter = 0; iter < ctx->max_pack; iter++)
    {                           /* dloop : DO iter=1,ctx->max_pack */
        /* mps_dpolzer(ctx, &nit, &excep);  */
        mps_thread_dpolzer(ctx, &nit, &excep, required_zeros--);
        it_pack += nit;

        MPS_DEBUG(ctx, "DPE packet completed in %d iterations", nit);

        if (ctx->DOLOG)
            fprintf(ctx->logstr, "Packet %d iterations= %d\n", iter, nit);

        if (excep)
        {
            oldnclust = ctx->clusterization->n;

            /* cluster analysis */
            if (ctx->DOLOG)
                fprintf(ctx->logstr, "   DSOLVE: call dcluster\n");

            mps_dradii(ctx, ctx->active_poly, drad);
            mps_dcluster(ctx, drad, 2 * ctx->n);   /* Isolation factor */
            if (oldnclust == ctx->clusterization->n)
            {
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   DSOLVE:  CYCLE\n");
                continue;
            }
            else
            {
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   DSOLVE: call dmodify\n");
                mps_dmodify(ctx, true);

                if (iter == 0 && !d_after_f)
                    for (i = 0; i < ctx->n; i++)
                        if (ctx->approx_root[i]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                            ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;

                /* If the polynomial is not given in terms of its
                    * coeff. then skip the restart stage */
                if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
                {
                    /* choose new starting approximations only for new clusters */
                    if (ctx->DOLOG)
                        fprintf(ctx->logstr, "   DSOLVE: call drestart\n");
                    mps_drestart(ctx);
                }
                /* reset the status vector */
                for (j = 0; j < ctx->n; j++)
                    if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                        ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                for (j = 0; j < ctx->n; j++)
                    ctx->again_old[j] = ctx->approx_root[j]->again;

                /* update 'again' */
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   DSOLVE: call update\n");
                mps_update(ctx);
                /* adjust 'again'
                    * This is needed since we are between two packets
                    */
                for (i = 0; i < ctx->n; i++)
                    if (!ctx->again_old[i])
                        ctx->approx_root[i]->again = false;
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "   DSOLVE: call checkstop\n");
                if (mps_check_stop(ctx))
                    goto dsolve_final_cleanup;
            }
        }
        else
            break;
    }
    if (iter == ctx->max_pack)
    {
        mps_dump(ctx);
        mps_error(ctx,
            "DPE: reached the maximum number of packet iterations");
    }
    /* Otherwise exit since all the approximations are
        * in the root neighbourhood, except for the ones that cannot be
        * represented as double.
        */

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "DPE: nit=%d\n", nit);

    /* Update */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   DSOLVE: now update: call dcluster\n");
    oldnclust = ctx->clusterization->n;

    mps_dradii(ctx, ctx->active_poly, drad);
    mps_dcluster(ctx, drad, 2 * ctx->n);   /* Isolation factor */

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "   DSOLVE: now call dmodify\n");
    mps_dmodify(ctx, true);

    /* reset the status vector */
    for (j = 0; j < ctx->n; j++)
        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

dsolve_final_cleanup:
    rdpe_vfree(drad);
}

/**
    * @brief Multiprecision version of <code>fsolve()</code>.
    */
MPS_PRIVATE void
    mps_msolve(mps_context* ctx)
{
    int iter, nit, oldnclust, i, j, it_pack, required_zeros = ctx->n;
    mps_boolean excep;
    int nzc;
    rdpe_t* drad = rdpe_valloc(ctx->n);

    /* == 1 == Initialize variables */
    it_pack = 0;

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE: call restart\n");
    if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
        mps_mrestart(ctx);
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE: call update1\n");
    mps_update(ctx);
    if (ctx->DOLOG)
    {
        fprintf(ctx->logstr, "  MSOLVE: again after update = ");
        for (i = 0; i < ctx->n; i++)
            fprintf(ctx->logstr, "%d", ctx->approx_root[i]->again);
        fprintf(ctx->logstr, "\n");
    }

    for (i = 0; i < ctx->n; i++)
        if (ctx->approx_root[i]->again)
            ctx->approx_root[i]->wp = ctx->mpwp;

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE: call checkstop\n");
    if (mps_check_stop(ctx))
    {
        oldnclust = ctx->clusterization->n;

        mps_mmodify(ctx, true);

        /* reset the status vector */
        for (j = 0; j < ctx->n; j++)
            if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

        goto msolve_final_cleanup;
    }
    nzc = 0;
    if (ctx->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE &&
        !ctx->output_config->multiplicity &&
        !ctx->output_config->root_properties)
        for (i = 0; i < ctx->n; i++)
            if (MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
                nzc++;
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE: nzc=%d\n", nzc);

    if (nzc == ctx->n)
    {
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "  MSOLVE: call mmodify and return\n");
        mps_mmodify(ctx, true);

        /* reset the status vector */
        for (j = 0; j < ctx->n; j++)
            if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

        goto msolve_final_cleanup;
    }

    /* Perform ctx->max_pack  packets of Aberth's iterations */
    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE: Perform packets of Aberth\n");

    for (iter = 0; iter < ctx->max_pack; iter++)
    {                           /* mloop : DO iter=1,ctx->max_pack */
        if (ctx->DOLOG)
        {
            fprintf(ctx->logstr, "  MSOLVE: packet= %d\n", iter);
            fprintf(ctx->logstr, "  MSOLVE: again before mpolzer =");
            for (i = 0; i < ctx->n; i++)
                fprintf(ctx->logstr, "%d", ctx->approx_root[i]->again);
            fprintf(ctx->logstr, "\n");
            fprintf(ctx->logstr, "  MSOLVE: call mpolzer\n");
        }
        /* mps_mpolzer(ctx, &nit, &excep); */
        mps_thread_mpolzer(ctx, &nit, &excep, required_zeros--);

        if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
            mps_dump(ctx);

        if (ctx->DOLOG)
            fprintf(ctx->logstr, "  MSOLVE: Packet %d: iterations= %d\n", iter, nit);

        it_pack += nit;
        nzc = 0;
        if (ctx->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE &&
            !ctx->output_config->multiplicity &&
            !ctx->output_config->root_properties)  /* DARIO APRILE 98 */
            for (i = 0; i < ctx->n; i++)
                if (MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
                    nzc++;
        if (ctx->DOLOG)
            fprintf(ctx->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
        if (nzc == ctx->n)
        {
            if (ctx->DOLOG)
                fprintf(ctx->logstr, "  MSOLVE: call mmodify and return\n");
            mps_mmodify(ctx, true);

            /* reset the status vector */
            for (j = 0; j < ctx->n; j++)
                if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                    ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
            if (ctx->DOLOG)
            {
                MPS_DEBUG(ctx, "Dumping status");
                for (j = 0; j < ctx->n; j++)
                    MPS_DEBUG(ctx, " %4d: %s", j,
                        MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[j]->status));
            }

            goto msolve_final_cleanup;
        }

        if (ctx->DOLOG)
            fprintf(ctx->logstr, "  MSOLVE: isolated %d roots excep=%d\n", nzc,
                excep);

        if (excep)
        {
            oldnclust = ctx->clusterization->n;

            /* cluster analysis */
            if (ctx->DOLOG)
                fprintf(ctx->logstr, "  MSOLVE: call mcluster\n");

            mps_mradii(ctx, ctx->active_poly, drad);
            mps_mcluster(ctx, drad, 2 * ctx->n);   /* Isolation factor */

            ctx->newtis_old = ctx->newtis;
            if (ctx->newtis == 0)
                mps_mnewtis(ctx);
            if (ctx->DOLOG)
                fprintf(ctx->logstr,
                    "  MSOLVE: newtis_old=%d, newtis=%d, oldncl=%d, ctx->nclust=%ld\n",
                    ctx->newtis_old, ctx->newtis, oldnclust, ctx->clusterization->n);

            if (oldnclust == ctx->clusterization->n && !(ctx->newtis == 1 && ctx->newtis_old == 0))
                /*#             if(&& iter != 0) AGO99 */
            {
                /*#D !newtis */
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "  MSOLVE: CYCLE\n");
                continue;
            }
            else
            {
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "  MSOLVE: call modify\n");
                mps_mmodify(ctx, true);
                if (iter == 0)
                    /* if first packet: reset the status vector */
                    for (j = 0; j < ctx->n; j++)
                        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                if (ctx->DOLOG)
                {
                    MPS_DEBUG(ctx, "Dumping status");
                    for (j = 0; j < ctx->n; j++)
                        MPS_DEBUG(ctx, " %4d: %s", j,
                            MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[j]->status));
                }
                /* If the polynomial is not given in terms of its coeff. then
                    * skip the restart stage */
                if (MPS_IS_MONOMIAL_POLY(ctx->active_poly))
                {
                    /* choose new starting approximations only for new clusters */
                    if (ctx->DOLOG)
                        fprintf(ctx->logstr,
                            "  MSOLVE: call mrestart for new clusters\n");
                    mps_mrestart(ctx);
                }
                /* reset the ctx->status vector */
                for (j = 0; j < ctx->n; j++)
                {
                    if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                        ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;
                    ctx->again_old[j] = ctx->approx_root[j]->again;
                }
                /* update 'again' */
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "  MSOLVE: call update2 : ");
                mps_update(ctx);
                if (ctx->DOLOG)
                {
                    fprintf(ctx->logstr, "  MSOLVE: again = ");
                    for (j = 0; j < ctx->n; j++)
                        fprintf(ctx->logstr, "%d", ctx->approx_root[j]->again);
                    fprintf(ctx->logstr, "\n");
                }
                /* adjust 'again'  This is needed since we are between two packets */
                for (i = 0; i < ctx->n; i++)
                    if (!ctx->again_old[i])
                        ctx->approx_root[i]->again = false;
                if (ctx->DOLOG)
                {
                    fprintf(ctx->logstr, "  MSOLVE: adjusted again = ");
                    for (j = 0; j < ctx->n; j++)
                        fprintf(ctx->logstr, "%d", ctx->approx_root[j]->again);
                    fprintf(ctx->logstr, "\n");
                }
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "  MSOLVE: call checkstop\n");
                if (mps_check_stop(ctx))
                {
                    mps_mmodify(ctx, true);

                    /* reset the ctx->status vector */
                    for (j = 0; j < ctx->n; j++)
                        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

                    goto msolve_final_cleanup;
                }

                nzc = 0;
                if (ctx->output_config->search_set == MPS_SEARCH_SET_COMPLEX_PLANE &&
                    ctx->output_config->multiplicity &&
                    !ctx->output_config->root_properties)
                    for (i = 0; i < ctx->n; i++)
                        if (MPS_ROOT_STATUS_IS_COMPUTED(ctx->approx_root[i]->status))
                            nzc++;
                if (ctx->DOLOG)
                    fprintf(ctx->logstr, "  MSOLVE: check again nzc=%d\n", nzc);
                if (nzc == ctx->n)
                {
                    if (ctx->DOLOG)
                        fprintf(ctx->logstr, "  MSOLVE: call mmodify and return");
                    mps_mmodify(ctx, true);

                    /* reset the ctx->status vector */
                    for (j = 0; j < ctx->n; j++)
                        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
                            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

                    goto msolve_final_cleanup;
                }
            }
        }
        else
            break;
    }

    if (iter == ctx->max_pack)
    {
        mps_error(ctx, "MP: reached the maximum number of packet iteration");
    }

    if (ctx->DOLOG)
    {
        fprintf(ctx->logstr, "  MSOLVE: MP: nit= %d\n", nit);
        fprintf(ctx->logstr, "  MSOLVE: call mcluster\n");
    }

    mps_mradii(ctx, ctx->active_poly, drad);
    mps_mcluster(ctx, drad, 2 * ctx->n);   /* Isolation factor */

    if (ctx->DOLOG)
        fprintf(ctx->logstr, "  MSOLVE:  call mmodify\n");

    oldnclust = ctx->clusterization->n;

    mps_mmodify(ctx, true);

    for (j = 0; j < ctx->n; j++)
        if (ctx->approx_root[j]->status == MPS_ROOT_STATUS_NEW_CLUSTERED)
            ctx->approx_root[j]->status = MPS_ROOT_STATUS_CLUSTERED;

    if (ctx->DOLOG)
    {
        MPS_DEBUG(ctx, "Dumping status: ");
        for (j = 0; j < ctx->n; j++)
            MPS_DEBUG(ctx, " %4d: %s", j,
                MPS_ROOT_STATUS_TO_STRING(ctx->approx_root[j]->status));
    }
    mps_update(ctx);

    if (ctx->DOLOG)
    {
        for (j = 0; j < ctx->n; j++)
            fprintf(ctx->logstr, "%d", ctx->approx_root[j]->again);
        fprintf(ctx->logstr, "\n");
    }

msolve_final_cleanup:
    rdpe_vfree(drad);
}

/**
    * @brief Multiprecision versione of <code>fpolzer()</code>.
    */
MPS_PRIVATE void
    mps_mpolzer(mps_context* ctx, int* it, mps_boolean* excep)
{
    int nzeros, i, iter, l;
    mpc_t corr, abcorr;
    rdpe_t eps, rad1_rdpe_t, rtmp;
    cdpe_t ctmp;
    mps_polynomial* p = ctx->active_poly;

    mpc_init2(abcorr, ctx->mpwp);
    mpc_init2(corr, ctx->mpwp);

    rdpe_mul_d(eps, ctx->mp_epsilon, (double)4 * ctx->n);

    /* initialize the iteration counter */
    *it = 0;
    *excep = false;

    /* count the number of approximations in the root neighbourhood */
    nzeros = 0;
    for (i = 0; i < ctx->n; i++)
        if (!ctx->approx_root[i]->again)
            nzeros++;
    if (nzeros == ctx->n)
        goto endfun;

    mps_cluster_item* c_item;
    mps_cluster* cluster;
    mps_root* root;

    /* Start Aberth's iterations */
    for (iter = 0; iter < ctx->max_it; iter++)
    {                           /* do_iter: */
        for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item)
        {                       /* do_clust: */
            cluster = c_item->cluster;
            for (root = cluster->first_root; root != NULL; root = root->next_root)
            {                   /* do_indice: */
                l = root->k;

                MPS_DEBUG(ctx, "Iterating on root %d, iter %d", l, iter);

                if (ctx->approx_root[l]->again)
                {
                    (*it)++;
                    /* sparse/dense polynomial */
                    rdpe_set(rad1_rdpe_t, ctx->approx_root[l]->drad);
                    mps_polynomial_mnewton(ctx, p, ctx->approx_root[l], corr, ctx->mpwp);
                    if (iter == 0 && !ctx->approx_root[l]->again && rdpe_gt(ctx->approx_root[l]->drad,
                        rad1_rdpe_t)
                        && rdpe_ne(rad1_rdpe_t, rdpe_zero))
                        rdpe_set(ctx->approx_root[l]->drad, rad1_rdpe_t);

                    /************************************************
                        The above condition is needed to cope with the case
                        where at the first iteration the starting point is
                        already in the root neighbourhood and the actually
                        computed radius is too big since the value of the
                        first derivative is too small.
                        In this case the previous radius bound, obtained by
                        means of Rouche' is more reliable and strict
                        ***********************************************/

                    if (ctx->approx_root[l]->again ||
                        /* the correction is performed only if iter!=1 or rad[l]!=rad1_rdpe_t */
                        iter != 0
                        || rdpe_ne(ctx->approx_root[l]->drad, rad1_rdpe_t))
                    {
                        mps_maberth_s(ctx, ctx->approx_root[l], cluster, abcorr);
                        mpc_mul_eq(abcorr, corr);
                        mpc_neg_eq(abcorr);
                        mpc_add_eq_ui(abcorr, 1, 0);
                        mpc_div(abcorr, corr, abcorr);
                        mpc_sub_eq(ctx->approx_root[l]->mvalue, abcorr);
                        mpc_get_cdpe(ctmp, abcorr);
                        cdpe_mod(rtmp, ctmp);
                        rdpe_add_eq(ctx->approx_root[l]->drad, rtmp);
                    }

                    /* check for new approximated roots */
                    if (!ctx->approx_root[l]->again)
                    {
                        nzeros++;
                        if (nzeros == ctx->n)
                            goto endfun;
                    }

                    MPS_DEBUG_MPC(ctx, 15, ctx->approx_root[l]->mvalue, "ctx->mroot[%d]", l);
                    MPS_DEBUG_RDPE(ctx, ctx->approx_root[l]->drad, "ctx->drad[%d]", l);
                }
            }
        }
        if (nzeros == ctx->n)
            goto endfun;
    }
    *excep = true;

endfun:                        /* free local MP variables */
    mpc_clear(corr);
    mpc_clear(abcorr);
}
