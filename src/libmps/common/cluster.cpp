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


#include <string.h>
#include <float.h>
#include <assert.h>
#include <mps/mps.h>
#include <math.h>

/**
* @brief Get an empty mps_cluster, with no roots.
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_cluster*
    mps_cluster_empty(mps_context* ctx)
{
    mps_new_obj(mps_cluster, cluster, sizeof(mps_cluster));

    cluster->first_root = NULL;
    cluster->n = 0;

    mps_mutex_init(cluster->list_lock);

    return cluster;
}

/**
* @brief Create a cluster containing only the selected root.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param root_index The root that must be in the cluster.
*/
mps_cluster*
    mps_cluster_with_root(mps_context* ctx, long int root_index)
{
    mps_new_obj(mps_cluster, cluster, sizeof(mps_cluster));

    mps_new_defined_obj(mps_root, cluster->first_root, sizeof(mps_root));
    cluster->n = 1;

    cluster->first_root->k = root_index;
    cluster->first_root->next_root = NULL;
    cluster->first_root->prev_root = NULL;

    mps_mutex_init(cluster->list_lock);

    return cluster;
}

/**
* @brief Free a previously allocated cluster with all the roots in
* it.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param cluster The cluster to free.
*/
void
    mps_cluster_free(mps_context* ctx, mps_cluster* cluster)
{
    mps_root* root = cluster->first_root;
    mps_root* old_root;

    /* Free all the roots in the cluster */
    while (root)
    {
        old_root = root;
        root = root->next_root;
        mps_delete_obj(old_root);
    }

    mps_mutex_destroy(cluster->list_lock);
    mps_delete_obj(cluster);
}

/**
* @brief Insert a root in a cluster.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param cluster The cluster in which the root must be inserted.
* @param root_index The index of the root to insert.
*/
mps_root*
    mps_cluster_insert_root(mps_context* ctx,
        mps_cluster* cluster,
        long int root_index)
{
    mps_new_obj(mps_root, root, sizeof(mps_root));

    /* Inserting root in the starting of the cluster */
    root->k = root_index;
    root->prev_root = NULL;

    root->next_root = cluster->first_root;
    cluster->n++;

    if (cluster->first_root)
        cluster->first_root->prev_root = root;

    cluster->first_root = root;

    return root;
}

/**
* @brief Remove a root from a cluster.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param cluster The cluster from which the root must be removed.
* @param root The root to remove.
*
* Please note the the root specified must be in the cluster, otherwise
* an assertion error or segmentation fault will be triggered.
*/
void
    mps_cluster_remove_root(mps_context* ctx, mps_cluster* cluster, mps_root* root)
{
    mps_root* prev_root = root->prev_root;
    mps_root* next_root = root->next_root;

    if (prev_root)
        prev_root->next_root = next_root;
    if (next_root)
        next_root->prev_root = prev_root;

    if (cluster->first_root == root)
        cluster->first_root = root->next_root;

    /* Decrease root count */
    cluster->n--;

    /* Free the root */
    mps_delete_obj(root);
}

/**
* @brief Join two cluster in one big cluster containing the roots of
* both. Please note that the cluster must not overlap.
*
* @param ctx The <code>mps_context</code> of the current computation.
* @param cluster_a The first cluster
* @param cluster_b The second cluster
* @return A new cluster containing the roots of both.
*/
mps_cluster*
    mps_cluster_join(mps_context* ctx, mps_cluster* cluster_a, mps_cluster* cluster_b)
{
    mps_root* root;
    mps_cluster* small_cluster;
    mps_cluster* big_cluster;

    mps_cluster* new_cluster = mps_cluster_empty(ctx);

    if (cluster_a->n < cluster_b->n)
    {
        small_cluster = cluster_a;
        big_cluster = cluster_b;
    }
    else
    {
        small_cluster = cluster_b;
        big_cluster = cluster_a;
    }

    root = small_cluster->first_root;
    while (root->next_root != NULL)
        root = root->next_root;
    root->next_root = big_cluster->first_root;

    new_cluster->first_root = small_cluster->first_root;

    new_cluster->n = small_cluster->n + big_cluster->n;

    return new_cluster;
}

/**
* @brief Create a new empty clusterization.
*
* @param ctx The <code>mps_context</code> of the current computation.
*/
mps_clusterization*
    mps_clusterization_empty(mps_context* ctx)
{
    mps_new_obj(mps_clusterization, clusterization, sizeof(mps_clusterization));

    clusterization->n = 0;
    clusterization->first_cluster_item = NULL;
    return clusterization;
}

/**
* @brief Insert a new cluster into a root clusterization.
* @param ctx The <code>mps_context</code> of the current computation.
* @param c The clusterization in which the cluster should be inserted.
* @param cluster The cluster that should be inserted.
*/
mps_cluster_item*
    mps_clusterization_insert_cluster(mps_context* ctx, mps_clusterization* clusterization, mps_cluster* cluster)
{
    mps_new_obj(mps_cluster_item, item, sizeof(mps_cluster_item));

    /* Set previous item as NULL and next item as the first now */
    item->prev_cluster_item = NULL;
    item->next_cluster_item = clusterization->first_cluster_item;
    item->detached_cluster_item = NULL;

    /* Store the pointer to the cluster */
    item->cluster = cluster;

    /* The first cluster is the one we're inserting */
    clusterization->first_cluster_item = item;

    /* The previous item of before was NULL, now should be the cluster that
    * we are inserting. */
    if (item->next_cluster_item)
        item->next_cluster_item->prev_cluster_item = item;

    clusterization->n++;

    return item;
}

/**
* @brief Pop out a cluster from a clusterization.
* @param ctx The <code>mps_context</code> of the current computation.
* @param c The clusterization from which the cluster_item should be popped.
* @param cluster_item The cluster item to remove.
*/
void
    mps_clusterization_pop_cluster(mps_context* ctx, mps_clusterization* clusterization, mps_cluster_item* cluster_item)
{
    mps_cluster_item* next_cluster_item = cluster_item->next_cluster_item;
    mps_cluster_item* prev_cluster_item = cluster_item->prev_cluster_item;

    if (prev_cluster_item)
        prev_cluster_item->next_cluster_item = next_cluster_item;
    if (next_cluster_item)
        next_cluster_item->prev_cluster_item = prev_cluster_item;

    if (clusterization->first_cluster_item == cluster_item)
    {
        clusterization->first_cluster_item = next_cluster_item;
    }
    clusterization->n--;
}

/**
* @brief Remove a cluster item from a clusterization, freeing it.
* @param ctx The <code>mps_context</code> of the current computa
tion.
* @param c The clusterization from where the cluster_item should be removed.
* @param cluster_item The cluster item to remove.
*/
void
    mps_clusterization_remove_cluster(mps_context* ctx, mps_clusterization* clusterization, mps_cluster_item* cluster_item)
{
    mps_clusterization_pop_cluster(ctx, clusterization, cluster_item);
    mps_cluster_free(ctx, cluster_item->cluster);
    mps_delete_obj(cluster_item);
}

/**
* @brief Free a clusterization and all the cluster in it.
* @param ctx The <code>mps_context</code> of the current computation.
* @param c The clusterization to free.
*/
void
    mps_clusterization_free(mps_context* ctx, mps_clusterization* clusterization)
{
    mps_cluster_item* cluster_item = clusterization->first_cluster_item;
    mps_cluster_item* next_cluster_item;

    while (cluster_item != NULL)
    {
        mps_cluster_free(ctx, cluster_item->cluster);
        next_cluster_item = cluster_item->next_cluster_item;
        mps_delete_obj(cluster_item);
        cluster_item = next_cluster_item;
    }

    mps_delete_obj(clusterization);
}

/**
* @brief Reset cluster structure information contained in <code>ctx</code>. After
* the call to this routine the roots will be considered as a unique big cluster,
* discarding every information present before.
*
* @param ctx the mps_context pointer.
*/
void
    mps_cluster_reset(mps_context* ctx)
{
    /* Reset cluster status of the roots */
    int i;
    mps_cluster* cluster;

    for (i = 0; i < ctx->n; i++)
    {
        ctx->approx_root[i]->status = MPS_ROOT_STATUS_CLUSTERED;
        ctx->approx_root[i]->attrs = MPS_ROOT_ATTRS_NONE;
        ctx->approx_root[i]->inclusion = MPS_ROOT_INCLUSION_UNKNOWN;
    }

    if (ctx->clusterization != NULL)
        mps_clusterization_free(ctx, ctx->clusterization);
    ctx->clusterization = mps_clusterization_empty(ctx);

    /* Fill in the roots */
    cluster = mps_cluster_empty(ctx);
    for (i = 0; i < ctx->n; i++)
        mps_cluster_insert_root(ctx, cluster, i);

    mps_clusterization_insert_cluster(ctx, ctx->clusterization, cluster);
}

void
    mps_clusterization_detach_clusters(mps_context* ctx, mps_clusterization* clusterization)
{
    MPS_DEBUG_THIS_CALL(ctx);

    /* Disable this function since it is not working as it should. */
    /* The problem, now, is that more than a root could be removed from
    * a cluster and will be then marked as isolated even if it is only isolated
    * from the base cluster and from the other detached roots. */
    return;

    mps_cluster_item* cluster_item;
    rdpe_t rtmp;
    int k;

    for (cluster_item = clusterization->first_cluster_item; cluster_item != NULL; cluster_item = cluster_item->next_cluster_item)
    {
        mps_root* root;

        /* Skip isolated clusters */
        if (cluster_item->cluster->n == 1)
            continue;

        /* Scan the cluster for quasi approximated roots */
        root = cluster_item->cluster->first_root;
        while (root != NULL)
        {
            k = root->k;
            mpc_rmod(rtmp, ctx->approx_root[k]->mvalue);

            /* We need a complex condition here since the heuristic used to determine if a root
            * is a simple root in a cluster is based on Newton radii.
            * These have different behavious based on the algorithm that has been selected, so
            * we introduce here two different guesses that work in each one. */
            if (((ctx->algorithm == MPS_ALGORITHM_STANDARD_MPSOLVE) &&
                ((rdpe_Esp(rtmp) - rdpe_Esp(ctx->approx_root[k]->drad) > ctx->mpwp / ::sqrt(cluster_item->cluster->n) + 1) ||
                    (ctx->approx_root[k]->status == MPS_ROOT_STATUS_APPROXIMATED_IN_CLUSTER))) ||
                ((ctx->algorithm == MPS_ALGORITHM_SECULAR_GA) &&
                    (rdpe_Esp(rtmp) - rdpe_Esp(ctx->approx_root[k]->drad) > ctx->mpwp - 4)))
            {
                if (ctx->debug_level & MPS_DEBUG_CLUSTER)
                {
                    MPS_DEBUG(ctx, "Temporary removing root %d from its cluster since it is quasi approximated", k);
                }

                mps_cluster* detached_cluster = mps_cluster_with_root(ctx, k);
                mps_root* next_root = root->next_root;
                mps_cluster_remove_root(ctx, cluster_item->cluster, root);
                root = next_root;

                /* Insert the cluster in the clusterization */
                mps_cluster_item* new_cluster_item = mps_clusterization_insert_cluster(ctx,
                    ctx->clusterization,
                    detached_cluster);
                /* Set the new item as detached from the old one */
                new_cluster_item->detached_cluster_item = cluster_item;
            }
            else
                root = root->next_root;

            /* If we have left only an isolated roots stop checking this cluster */
            if (cluster_item->cluster->n == 1)
                break;
        }
    }
}


void
    mps_clusterization_reassemble_clusters(mps_context* ctx, mps_clusterization* clusterization)
{
    MPS_DEBUG_THIS_CALL(ctx);

    mps_cluster_item* cluster_item;

    cluster_item = ctx->clusterization->first_cluster_item;
    while (cluster_item != NULL)
    {
        mps_cluster_item* next_cluster_item = cluster_item->next_cluster_item;

        if (cluster_item->detached_cluster_item)
        {
            mps_cluster_insert_root(ctx, cluster_item->detached_cluster_item->cluster, cluster_item->cluster->first_root->k);
            mps_clusterization_remove_cluster(ctx, ctx->clusterization, cluster_item);
        }

        cluster_item = next_cluster_item;
    }
}

void
    mps_cluster_detachment_reset(mps_context* ctx)
{
    mps_clusterization_reassemble_clusters(ctx, ctx->clusterization);
}
