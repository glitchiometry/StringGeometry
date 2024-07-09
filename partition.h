#include "basics.h"

typedef struct
{
	aarray_int parts;
	array_int addr_0;
	array_int addr_1;
} partition;

void partition_init(partition *p, int size);
void partition_union(partition *p, int pi, int pii);
void partition_union_a2b(partition *p, int p1, int p0);
void partition_join(partition *p, int i, int j);
void partition_join_a2b(partition *p, int i, int j);
void extend_partition(partition *p);
void extend_partition_n(partition *p, int n);
void remove_from_partition(partition *p, int i);
void partition_set(partition *p, int i, int pi);
int partition_get(partition *p, int i);
void free_partition(partition *p);
void partition_add_partition(partition *p);
int partition_rep(partition *p, int i);
void fprintf_partition(partition *p, FILE *ofile);

