/*
* This file is part of MPSolve 3.2.2
*
* Copyright (C) 2001-2020, Dipartimento di Matematica "L. Tonelli", Pisa.
* License: http://www.gnu.org/licenses/gpl.html GPL version 3 or higher
*
* Authors:
*   Leonardo Robol <leonardo.robol@unipi.it>
*/

#include <mps/mps.h>
#include <string.h>

MPS_PRIVATE void
    mps_cluster_analysis(mps_context* ctx, mps_polynomial* p)
{
    switch (ctx->lastphase)
    {
    case float_phase:
    {
        double* radii = double_valloc(ctx->n);

        mps_fradii(ctx, p, radii);
        mps_fcluster(ctx, radii, 2 * ctx->n);
        mps_fmodify(ctx, false);

        mps_free(radii);
        break;
    }

    case dpe_phase:
    {
        rdpe_t* radii = rdpe_valloc(ctx->n);

        mps_dradii(ctx, p, radii);
        mps_dcluster(ctx, radii, 2 * ctx->n);
        mps_dmodify(ctx, false);

        mps_free(radii);
        break;
    }

    case mp_phase:
    {
        rdpe_t* radii = rdpe_valloc(ctx->n);

        mps_mradii(ctx, p, radii);
        mps_mcluster(ctx, radii, 2 * ctx->n);
        mps_mmodify(ctx, false);

        mps_free(radii);
        break;
    }

    case no_phase:
        break;
    }
}

/**
* This subroutine makes cluster analysis, i.e., detects
* overlapping disks, where two disks overlap if the distances
* of their centers is less than the sum of their radii
* multiplied by <code>nf</code>.
*
* Observe that \f$nf=1\f$ then this concept corresponds to overlapping,
* if \f$nf =2 \cdot n\f$, this concept corresponds to Newton isolation.
*
* This routine set the vector <code>clust</code> so that it
* contains the indices of the
* disks in each overlapping group, while  <code>punt[i]</code>
* points to the
* index of <code>clust</code> where the i-th group starts. Moreover
* <code>m_clust[i]</code>
* contains the  multiplicity of the i-th cluster.
* <code>nclust</code> is the
* number of clusters.
*
* @param ctx  The <code>mps_context</code> associated with the current
*           computaion.
* @param frad The vector of radii to use for cluster analysis.
* @param nf see above for a detailed description.
*/
MPS_PRIVATE void
    mps_fcluster(mps_context* ctx, double* frad, int nf)
{
    ctx->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

    /* We need to scan every cluster and make it in pieces, if possible */
    mps_clusterization* new_clusterization = mps_clusterization_empty(ctx);
    mps_cluster_item* item;
    int analyzed_roots = 0;
    int i, j;

    /* This value is set to false if the radius are not newton isolated
    * by means of the newton radii. */
    mps_boolean newton_isolation = true;

    /* Debug clusterization status if debugging was required */
    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        int i;
        MPS_DEBUG(ctx, "Debugging the radius and approximations obtained for the roots before cluster analysis");
        for (i = 0; i < ctx->n; i++)
        {
            MPS_DEBUG_CPLX(ctx, ctx->approx_root[i]->fvalue, "Root %d", i);
            MPS_DEBUG(ctx, "radius for root %4d: %e", i, frad[i]);
        }

        MPS_DEBUG(ctx, "Debugging cluster structure before cluster analysis");
        mps_debug_cluster_structure(ctx);
    }

    /*
    * Mark newton isolated roots as newton isolated.
    */
    double* newton_radii = mps_newv(double, ctx->n);
    for (i = 0; i < ctx->n; i++)
        newton_radii[i] = ctx->approx_root[i]->frad;

    for (i = 0; i < ctx->n; i++)
    {
        for (j = 0; j < ctx->n; j++)
        {
            if ((i != j) && mps_ftouchnwt(ctx, newton_radii, nf, i, j))
            {
                newton_isolation = false;
                break;
            }
            /* ctx->root_status[i] = MPS_ROOT_STATUS_NEWTON_ISOLATED; */
        }
    }

    mps_free(newton_radii);

    item = ctx->clusterization->first_cluster_item;
    while (item)
    {
        mps_cluster* cluster = item->cluster;
        mps_cluster_item* next_item = item->next_cluster_item;

        /* Keep isolated cluster isolated, moving them in the new clusterization. */
        if (cluster->n == 1)
        {
            mps_clusterization_insert_cluster(ctx, new_clusterization,
                mps_cluster_with_root(ctx, cluster->first_root->k));
            mps_clusterization_remove_cluster(ctx, ctx->clusterization, item);
            analyzed_roots++;
        }

        item = next_item;
    }

    /* Now do cluster analysis with the rest of the clusters, but only if
    * newton isolation is not guaranteed by means of the newton radii. */
    if (!newton_isolation)
    {
        /* if (MPS_INPUT_CONFIG_IS_USER (ctx->input_config))  */
        /*        {  */
        /*          mps_clusterization_free (ctx, new_clusterization);  */
        /*          return; */
        /*        } */

        item = ctx->clusterization->first_cluster_item;
        while (analyzed_roots < ctx->n)
        {
            /* Create a new cluster to be inserted in the cluster analysis */
            mps_root* base_root;
            mps_cluster* cluster = item->cluster;
            mps_cluster* new_cluster = mps_cluster_empty(ctx);

            while (cluster->n == 0)
            {
                item = item->next_cluster_item;
                cluster = item->cluster;
            }


            base_root = mps_cluster_insert_root(ctx, new_cluster, cluster->first_root->k);
            analyzed_roots++;
            mps_cluster_remove_root(ctx, cluster, cluster->first_root);

            /* Check if this root touches others root, and if new roots were added to the
            * cluster repeat the checks. */
            while (base_root)
            {
                /* mps_cluster_item * c_item; */
                mps_root* iter_root;

                /* for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_item) */
                {
                    /* mps_cluster * iter_cluster = c_item->cluster; */
                    mps_cluster* iter_cluster = cluster;

                    iter_root = iter_cluster->first_root;
                    while (iter_root)
                    {
                        if (mps_ftouchnwt(ctx, frad, nf, base_root->k, iter_root->k))
                        {
                            mps_root* next_root = iter_root->next_root;
                            mps_cluster_insert_root(ctx, new_cluster, iter_root->k);
                            mps_cluster_remove_root(ctx, iter_cluster, iter_root);
                            analyzed_roots++;
                            iter_root = next_root;
                        }
                        else
                            iter_root = iter_root->next_root;
                    }
                }

                base_root = base_root->prev_root;
            }

            /* Now insert the cluster in the new clusterization */
            mps_clusterization_insert_cluster(ctx, new_clusterization, new_cluster);

            /* Check if the new cluster is isolated and, in that case, set the gerschgorin
            * radius as inclusion radius if it's more conveniente than the old one.
            * In general the Gerschgorin radius cannot be used as inclusion radius, because
            * it may touch another radius and so it may be empty. */
            if (new_cluster->n == 1)
            {
                int k = new_cluster->first_root->k;
                double new_rad;

                new_rad = cplx_mod(ctx->approx_root[k]->fvalue) * 4.0f * DBL_EPSILON + frad[k];

                /* Check if the computed radius is more convenient than the old one.  */
                /*      If that's the case, apply it as inclusion radius   */
                if (new_rad < ctx->approx_root[k]->frad)
                    ctx->approx_root[k]->frad = new_rad;
            }
        }
    }

    if (newton_isolation)
    {
        mps_clusterization_free(ctx, new_clusterization);
        new_clusterization = mps_clusterization_empty(ctx);

        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            MPS_DEBUG(ctx, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

        for (i = 0; i < ctx->n; i++)
        {
            mps_clusterization_insert_cluster(ctx, new_clusterization,
                mps_cluster_with_root(ctx, i));
        }
    }

    /* Set the new clusterizaition in the mps_context */
    mps_clusterization_free(ctx, ctx->clusterization);
    ctx->clusterization = new_clusterization;

    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        MPS_DEBUG(ctx, "Debugging cluster structure after cluster analysis");
        mps_debug_cluster_structure(ctx);
    }
}

/**
* @brief Perform cluster analysis to each existing cluster by
* applying <code>mps_xcluster</code> to each existing cluster.
*
* Rebuild the vectors <code>ctx->clust</code>,
* <code>ctx->punt</code>, and the integer <code>ctx->nclust</code>.
*
*/
MPS_PRIVATE void
    mps_dcluster(mps_context* ctx, rdpe_t* drad, int nf)
{
    ctx->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

    /* We need to scan every cluster and make it in pieces, if possible */
    mps_clusterization* new_clusterization = mps_clusterization_empty(ctx);
    mps_cluster_item* item;
    int i, j;

    /* This value is set to false if the radius are not newton isolated
    * by means of the newton radii. */
    mps_boolean newton_isolation = true;

    /* Debug clusterization status if debugging was required */
    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        int i;
        MPS_DEBUG(ctx, "Debugging the radius and approximations obtained for the roots before cluster analysis");
        for (i = 0; i < ctx->n; i++)
        {
            MPS_DEBUG_CDPE(ctx, ctx->approx_root[i]->dvalue, "Root %d", i);
            MPS_DEBUG_RDPE(ctx, drad[i], "radius for root %4d", i);
        }

        MPS_DEBUG(ctx, "Debugging cluster structure before cluster analysis");
        mps_debug_cluster_structure(ctx);
    }

    int analyzed_roots = 0;

    /* Do a first check of clusterization using the newton
    * radii. These are not valid to perform cluster analysis in
    * general, but can be used if they provide *COMPLETE* Newton
    * isolation. */
    rdpe_t* newton_radii = rdpe_valloc(ctx->n);
    for (i = 0; i < ctx->n; i++)
        rdpe_set(newton_radii[i], ctx->approx_root[i]->drad);

    for (i = 0; i < ctx->n; i++)
    {
        for (j = 0; j < ctx->n; j++)
        {
            if ((i != j) && mps_dtouchnwt(ctx, newton_radii, nf, i, j))
            {
                newton_isolation = false;
                break;
            }
        }
    }

    rdpe_vfree(newton_radii);

    /* If newton isolation has not been reached check with Gerschgorin */
    {
        /* if (MPS_INPUT_CONFIG_IS_USER (ctx->input_config))  */
        /*        {  */
        /*          mps_clusterization_free (ctx, new_clusterization);  */
        /*          return;  */
        /*        } */

        item = ctx->clusterization->first_cluster_item;
        while (item)
        {
            mps_cluster* cluster = item->cluster;
            mps_cluster_item* next_item = item->next_cluster_item;

            /* Keep isolated cluster isolated, moving them in the new clusterization. */
            if (cluster->n == 1)
            {
                mps_clusterization_insert_cluster(ctx, new_clusterization,
                    mps_cluster_with_root(ctx, cluster->first_root->k));
                mps_clusterization_remove_cluster(ctx, ctx->clusterization, item);
                analyzed_roots++;
            }

            item = next_item;
        }

        /* Now do cluster analysis with the rest of the clusters. */
        item = ctx->clusterization->first_cluster_item;
        while (analyzed_roots < ctx->n)
        {
            /* Create a new cluster to be inserted in the cluster analysis */
            mps_root* base_root;
            mps_cluster* cluster = item->cluster;
            mps_cluster* new_cluster = mps_cluster_empty(ctx);

            while (cluster->n == 0)
            {
                item = item->next_cluster_item;
                cluster = item->cluster;
            }

            base_root = mps_cluster_insert_root(ctx, new_cluster, cluster->first_root->k);
            analyzed_roots++;
            mps_cluster_remove_root(ctx, cluster, cluster->first_root);

            /* Check if this root touches others root, and if new roots were added to the
            * cluster repeat the checks. */
            while (base_root)
            {
                /* mps_cluster_item * c_item; */
                mps_root* iter_root;

                /* for (c_item = ctx->clusterization->first_cluster_item; c_item != NULL; c_item = c_item->next_cluster_item) */
                {
                    /* mps_cluster * iter_cluster = c_item->cluster; */
                    mps_cluster* iter_cluster = cluster;

                    iter_root = iter_cluster->first_root;
                    while (iter_root)
                    {
                        if (mps_dtouchnwt(ctx, drad, nf, base_root->k, iter_root->k))
                        {
                            mps_root* next_root = iter_root->next_root;
                            mps_cluster_insert_root(ctx, new_cluster, iter_root->k);
                            mps_cluster_remove_root(ctx, iter_cluster, iter_root);
                            analyzed_roots++;
                            iter_root = next_root;
                        }
                        else
                            iter_root = iter_root->next_root;
                    }
                }

                base_root = base_root->prev_root;
            }

            /* Now insert the cluster in the new clusterization */
            mps_clusterization_insert_cluster(ctx, new_clusterization, new_cluster);

            /* Check if the new cluster is isolated and, in that case, set the gerschgorin
            * radius as inclusion radius if it's more conveniente than the old one.
            * In general the Gerschgorin radius cannot be used as inclusion radius, because
            * it may touch another radius and so it may be empty. */
            if (new_cluster->n == 1)
            {
                int k = new_cluster->first_root->k;
                rdpe_t new_rad;

                cdpe_mod(new_rad, ctx->approx_root[k]->dvalue);
                rdpe_mul_eq_d(new_rad, 4 * DBL_EPSILON);
                rdpe_add_eq(new_rad, drad[k]);


                /* Check if the computed radius is more convenient than the old one.
                If that's the case, apply it as inclusion radius */
                if (rdpe_lt(new_rad, ctx->approx_root[k]->drad))
                    rdpe_set(ctx->approx_root[k]->drad, new_rad);
            }
        }
    }


    if (newton_isolation)
    {
        mps_clusterization_free(ctx, new_clusterization);
        new_clusterization = mps_clusterization_empty(ctx);

        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            MPS_DEBUG(ctx, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

        for (i = 0; i < ctx->n; i++)
        {
            mps_clusterization_insert_cluster(ctx, new_clusterization,
                mps_cluster_with_root(ctx, i));
        }
    }

    /* Set the new clusterizaition in the mps_context */
    mps_clusterization_free(ctx, ctx->clusterization);
    ctx->clusterization = new_clusterization;

    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        MPS_DEBUG(ctx, "Debugging cluster structure after cluster analysis");
        mps_debug_cluster_structure(ctx);
    }
}

struct _mps_cluster_worker_data {
    mps_context* ctx;
    mps_cluster* cluster;
    int* analyzed_roots;
    int base_root;
    int start_root;
    int end_root;
    rdpe_t* drad;
    int nf;
    mps_mutex_t* block_mutex;  // this will hold a ptr to a mutex not a 'new' mutex
    mps_cluster** original_clusters;
};

void*
    _mps_mcluster_worker(void* data_ptr)
{
    _mps_cluster_worker_data* data = (_mps_cluster_worker_data*)data_ptr;
    int i, n = 0;
    mps_root* first_root = NULL;
    mps_root* last_root = NULL;
    mps_cluster* cluster = data->original_clusters[data->base_root];

    for (i = data->start_root; i < data->end_root; i++)
    {
        if (!data->analyzed_roots[i] && (data->original_clusters[i] == cluster))
        {
            if (mps_mtouchnwt(data->ctx, data->drad, data->nf, data->base_root, i))
            {
                if (!data->analyzed_roots[i])
                {
                    data->analyzed_roots[i] = true;

                    if (first_root == NULL)
                    {
                        mps_new_obj(mps_root, first_root, sizeof(mps_root));
                        last_root = first_root;
                        last_root->next_root = first_root->next_root = last_root->prev_root = first_root->prev_root = NULL;
                        last_root->k = i;
                    }
                    else
                    {
                        mps_new_obj(mps_root, new_root, sizeof(mps_root));
                        new_root->next_root = first_root;
                        first_root->prev_root = new_root;
                        first_root = new_root;
                        new_root->prev_root = NULL;
                        new_root->k = i;
                    }

                    n++;
                }
            }
        }
    }

    if (n > 0)
    {
        mps_mutex_lock(data->cluster->list_lock);

        last_root->next_root = data->cluster->first_root;
        data->cluster->first_root->prev_root = last_root;
        data->cluster->first_root = first_root;
        data->cluster->n += n;

        mps_mutex_unlock(data->cluster->list_lock);
    }

    // unlock it only if I own it
    mps_mutex_guarded_unlock(*data->block_mutex);

    mps_delete_obj(data);

    return NULL;
}

/**
* @brief Perform cluster analysis to each existing cluster by
* applying <code>mps_xcluster</code> to each existing cluster.
*
*
* Rebuild the vectors <code>ctx->clust</code>,
* <code>ctx->punt</code>, and the integer <code>ctx->nclust</code>.
*
* @see mps_xcluster
*/
MPS_PRIVATE void
    mps_mcluster(mps_context* ctx, rdpe_t* drad, int nf)
{
    MPS_DEBUG_THIS_CALL(ctx);

    ctx->operation = MPS_OPERATION_CLUSTER_ANALYSIS;

    /* We need to scan every cluster and make it in pieces, if possible */
    mps_clusterization* new_clusterization = mps_clusterization_empty(ctx);
    mps_cluster_item* item;
    int i, j;

    /* This value is set to false if the radius are not newton isolated
    * by means of the newton radii. */
    mps_boolean newton_isolation = true;

    /* Debug clusterization status if debugging was required */
    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        int i;
        MPS_DEBUG(ctx, "Debugging the radius and approximations obtained for the roots before cluster analysis");
        for (i = 0; i < ctx->n; i++)
        {
            MPS_DEBUG_MPC(ctx, mpc_get_prec(ctx->approx_root[i]->mvalue), ctx->approx_root[i]->mvalue, "Root %d", i);
            MPS_DEBUG_RDPE(ctx, drad[i], "radius for root %4d", i);
        }

        MPS_DEBUG(ctx, "Debugging cluster structure before cluster analysis");
        mps_debug_cluster_structure(ctx);
    }

    /* Do a first check of clusterization using the newton
    * radii. These are not valid to perform cluster analysis in
    * general, but can be used if they provide *COMPLETE* Newton
    * isolation. */
    rdpe_t* newton_radii = rdpe_valloc(ctx->n);
    for (i = 0; i < ctx->n; i++)
        rdpe_set(newton_radii[i], ctx->approx_root[i]->drad);

    for (i = 0; i < ctx->n; i++)
    {
        for (j = 0; j < ctx->n; j++)
        {
            if ((i != j) && mps_mtouchnwt(ctx, newton_radii, nf, i, j))
            {
                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                    MPS_DEBUG(ctx, "Failing newton isolation on root %d and %d", i, j);

                newton_isolation = false;
                break;
            }
            /* ctx->root_status[i] = MPS_ROOT_STATUS_NEWTON_ISOLATED; */
        }

        if (!newton_isolation)
            break;
    }

    rdpe_vfree(newton_radii);

    /* Perform parallel analysis of the Gerschgorin disks. */
    int analyzed_roots = 0;
    int* already_analyzed_roots = mps_newv(int, ctx->n);
    mps_cluster** original_clusters = mps_newv(mps_cluster*, ctx->n);
    mps_root* root = NULL;

    memset(already_analyzed_roots, 0, sizeof(int) * ctx->n);

    int block_size = 128;
    int block_number = (ctx->n - 1) / block_size + 1;
    mps_new_array_obj(mps_mutex_t, block_mutexes, sizeof(mps_mutex_t), ctx->n);

    for (j = 0; j < block_number; j++)
        mps_mutex_init(block_mutexes[j]);

    item = ctx->clusterization->first_cluster_item;
    while (item)
    {
        mps_cluster* cluster = item->cluster;
        mps_root* my_root = cluster->first_root;
        while (my_root)
        {
            original_clusters[my_root->k] = cluster;
            my_root = my_root->next_root;
        }
        item = item->next_cluster_item;
    }

    while (analyzed_roots < ctx->n)
    {
        int j;

        /* Find the first not analyzed root */
        j = 0;
        while (already_analyzed_roots[j])
            j++;
        // fprintf (stderr, "New base root = %d\n", j);

        if (j > ctx->n)
            break;

        item = mps_clusterization_insert_cluster(ctx, new_clusterization, mps_cluster_with_root(ctx, j));
        root = item->cluster->first_root;

        already_analyzed_roots[j] = true;

        do
        {
            /* We need to check which other approximation touches our current base root
            * and add it to our cluster. */
            for (j = 0; j < block_number; j++)
            {
                mps_new_obj(_mps_cluster_worker_data, data, sizeof(_mps_cluster_worker_data));

                data->ctx = ctx;
                data->cluster = item->cluster;
                data->base_root = root->k;
                data->start_root = j * block_size;
                data->end_root = MIN((j + 1) * block_size, ctx->n);
                data->analyzed_roots = already_analyzed_roots;
                data->drad = drad;
                data->nf = nf;
                data->block_mutex = &block_mutexes[j];
                data->original_clusters = original_clusters;
                mps_mutex_guarded_lock(*data->block_mutex);
                mps_thread_pool_assign(ctx, ctx->pool, _mps_mcluster_worker, data);
            }

            mps_thread_yield();

            /* This is a cheap way of getting the next root in the cluster.
            * If it's already present (root->prev_root != NULL) just use it, otherwise
            * lock the mutexes that the clusters hold while analyzing a single
            * block. Getting lock implies that the cluster will be already
            * extended, thus root->prev_root != NULL unless there were no other
            * approximations inside the cluster in that block. In that case, try the
            * next until they are all analyzed. */
            for (j = 0; root->prev_root == NULL && j < block_number; j++)
            {
                // this may find we already own the lock and do nothing
                mps_mutex_guarded_lock(block_mutexes[j]);
                // release it only if I still own it (e.g.,  a worker doing it in this thread has not released it for me)
                mps_mutex_guarded_unlock(block_mutexes[j]);
            }

            analyzed_roots++;

        } while ((root = root->prev_root) != NULL);
    }

    for (j = 0; j < block_number; j++)
    {
        mps_mutex_destroy(block_mutexes[j]);
    }
    mps_delete_array_obj(block_mutexes);

    mps_free(already_analyzed_roots);
    mps_free(original_clusters);
    mps_clusterization_free(ctx, ctx->clusterization);
    ctx->clusterization = new_clusterization;

    for (item = ctx->clusterization->first_cluster_item; item != NULL; item = item->next_cluster_item)
    {
        mps_cluster* new_cluster = item->cluster;

        /* Check if the new cluster is isolated and, in that case, set the gerschgorin
        * radius as inclusion radius if it's more conveniente than the old one.
        * In general the Gerschgorin radius cannot be used as inclusion radius, because
        * it may touch another radius and so it may be empty. */
        if (new_cluster->n == 1)
        {
            int k = new_cluster->first_root->k;
            cdpe_t c;
            rdpe_t new_rad;

            /* Check if the computed radius is more convenient than the old one.
            If that's the case, apply it as inclusion radius */
            mpc_get_cdpe(c, ctx->approx_root[k]->mvalue);
            cdpe_mod(new_rad, c);
            rdpe_mul_eq(new_rad, ctx->mp_epsilon);
            rdpe_mul_eq_d(new_rad, 4.0f);
            rdpe_add_eq(new_rad, drad[k]);

            if (rdpe_lt(new_rad, ctx->approx_root[k]->drad))
            {
                rdpe_set(ctx->approx_root[k]->drad, new_rad);
                ctx->approx_root[k]->frad = rdpe_get_d(ctx->approx_root[k]->drad);
            }
        }
    }

    if (newton_isolation)
    {
        mps_clusterization_free(ctx, new_clusterization);
        new_clusterization = mps_clusterization_empty(ctx);

        if (ctx->debug_level & MPS_DEBUG_CLUSTER)
            MPS_DEBUG(ctx, "Reached isolation using Newton radii, so skipping every other check with Gerschgorin");

        for (i = 0; i < ctx->n; i++)
        {
            mps_clusterization_insert_cluster(ctx, new_clusterization,
                mps_cluster_with_root(ctx, i));
        }

        ctx->clusterization = new_clusterization;
    }

    if (ctx->debug_level & MPS_DEBUG_CLUSTER)
    {
        MPS_DEBUG(ctx, "Debugging cluster structure after cluster analysis");
        mps_debug_cluster_structure(ctx);
    }

    if (ctx->debug_level & MPS_DEBUG_APPROXIMATIONS)
    {
        mps_dump(ctx);
    }
}
