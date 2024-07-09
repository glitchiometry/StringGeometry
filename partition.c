#include "partition.h"

void partition_init(partition *p, int size)
{
	aarray_int_init(&((*p).parts), 0);
	array_int_init(&((*p).addr_0), size);
	array_int_init(&((*p).addr_1), size);
	(*p).addr_0.len = size;
	(*p).addr_1.len = size;
	for (int i = 0; i < size; i++)
	{
		(*p).addr_0.e[i] = -1;
		(*p).addr_1.e[i] = -1;
	}
}

void partition_union(partition *p, int pi, int pii)
{
	// OPT: consider estimating the total work in advance for each 
	// 	approach (e.g. copying parts.e[pi] to parts.e[pii] or vice versa)
	// 	to improve performance.
	if (pi < pii) partition_union_a2b(p, pii, pi);
	else partition_union_a2b(p, pi, pii);
}

void extend_partition(partition *p)
{
	add2array_int(&((*p).addr_0), -1);
	add2array_int(&((*p).addr_1), -1);
}

void extend_partition_n(partition *p, int n)
{
	// OPT: resize arrays and set elements manually
	for (int i = 0; i < n; i++) extend_partition(p);
}

void remove_from_partition(partition *p, int i)
{
	int i0 = (*p).addr_0.e[i];
	if (i0 > -1)
	{
		int i1 = (*p).addr_1.e[i];
		remove_array_int(&((*p).parts.e[i0]), i1);
		(*p).addr_0.e[i] = -1;
		(*p).addr_1.e[i] = -1; // OPT: one of these is unnecessary
	}
}

void partition_set(partition *p, int i, int pi)
{
	int pi0 = (*p).addr_0.e[i];
	if (pi0 != pi && pi < (*p).parts.len)
	{
		int pi0_a = (*p).addr_1.e[i];
		(*p).addr_0.e[i] = pi;
		(*p).addr_1.e[i] = (*p).parts.e[pi].len;
		add2array_int(&((*p).parts.e[pi]), i);
		if (pi0 > -1) remove_array_int(&((*p).parts.e[pi0]), pi0_a);
	}
}

void free_partition(partition *p)
{
	free_aarray_int(&((*p).parts));
	free_array_int(&((*p).addr_0));
	free_array_int(&((*p).addr_1));
}

void partition_add_partition(partition *p)
{
	extend_aarray_int(&((*p).parts));
}

int partition_rep(partition *p, int i)
{
  return (*p).addr_0.e[i] == -1 ? i : (*p).parts.e[(*p).addr_0.e[i]].e[0];
}

int partition_get(partition *p, int i)
{
	return (*p).addr_0.e[i];
}

void partition_union_a2b(partition *p, int p1, int p0)
{
  for (int ei = 0; ei < (*p).parts.e[p1].len; ei++)
	{
		int i = (*p).parts.e[p1].e[ei];
		(*p).addr_0.e[i] = p0;
		(*p).addr_1.e[i] = (*p).parts.e[p0].len;
		add2array_int(&((*p).parts.e[p0]), i);
	}
	remove_aarray_int(&((*p).parts), p1);
	if ((*p).parts.len > p1) 
	{
		for (int ei = 0; ei < (*p).parts.e[p1].len; ei++)
		{
			int i = (*p).parts.e[p1].e[ei];
			(*p).addr_0.e[i] = p1;
		}
	}
}

void partition_join_a2b(partition *p, int i, int j)
{
  	int ci = (*p).addr_0.e[i];
	int cj = (*p).addr_0.e[j];
	if (ci != -1 && cj != -1)
	{
		if (ci != cj) partition_union_a2b(p, ci, cj);
	}
	else
	{
		if (cj != -1) partition_set(p, i, cj);
		else if (ci != -1)
		  {
		    int ri = (*p).parts.e[ci].e[0];
		    (*p).parts.e[ci].e[0] = j;
		    partition_set(p, ri, ci);
		  }
		else
		{
			int pj = (*p).parts.len;
			partition_add_partition(p);
			partition_set(p, j, pj);
			partition_set(p, i, pj);
		}
	}
}

void partition_join(partition *p, int i, int j)
{
	int ci = (*p).addr_0.e[i];
	int cj = (*p).addr_0.e[j];
	if (ci != -1 && cj != -1)
	{
		if (ci != cj) partition_union(p, ci, cj);
	}
	else
	{
		if (cj != -1) partition_set(p, i, cj);
		else if (ci != -1) partition_set(p, j, ci);
		else
		{
			int pi = (*p).parts.len;
			partition_add_partition(p);
			partition_set(p, i, pi);
			partition_set(p, j, pi);
		}
	}
}

void fprintf_partition(partition *p, FILE *ofile)
{
  if (ofile != NULL) {}
  else ofile = stdout;
  fprintf_array_int(&((*p).addr_0), ofile);
}
