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

mps_list_element*
	mps_list_element_new(void* value)
{
	mps_new_obj(mps_list_element, list_element, sizeof(mps_list_element));

	list_element->value = value;
	list_element->previous_element = list_element->next_element = NULL;

	return list_element;
}

void
	mps_list_element_free(mps_list_element* list_element)
{
	mps_del_obj(list_element);
}

mps_list_element*
	mps_list_element_previous(mps_list_element* list_element)
{
	return list_element->previous_element;
}

mps_list_element*
	mps_list_element_next(mps_list_element* list_element)
{
	return list_element->next_element;
}

void*
	mps_list_element_value(mps_list_element* list_element)
{
	return (list_element != NULL) ? list_element->value : NULL;
}
