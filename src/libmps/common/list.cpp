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

/**
* @brief Create a new empty list.
*/
mps_list*
    mps_list_new(void)
{
    mps_new_obj(mps_list, l, sizeof(mps_list));

    l->first_element = l->last_element = NULL;
    l->size = 0;

    return l;
}

/**
* @brief Free a list and all the elements inside it.
*/
void
    mps_list_free(mps_list* list)
{
    /* First delete all the elements in the list */
    mps_list_element* list_element = NULL;
    for (list_element = mps_list_first(list); list_element != NULL; list_element = mps_list_element_next(list_element))
    {
        mps_list_element_free(list_element);
    }

    mps_del_obj(list);
}

/**
* @brief Return the number of elements in a list.
*/
int
    mps_list_size(mps_list* list)
{
    return list->size;
}

void
    mps_list_append(mps_list* list, mps_list_element* list_element)
{
    if (list->last_element)
    {
        list->last_element->next_element = list_element;
        list_element->previous_element = list->last_element;
        list_element->next_element = NULL;
        list->last_element = list_element;
    }
    else
    {
        list->first_element = list->last_element = list_element;
        list_element->previous_element = list_element->next_element = NULL;
    }

    list->size++;
}

mps_list_element*
    mps_list_first(mps_list* list)
{
    return list->first_element;
}

mps_list_element*
    mps_list_last_element(mps_list* list)
{
    return list->last_element;
}
