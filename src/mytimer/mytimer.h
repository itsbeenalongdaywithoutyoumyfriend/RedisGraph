#include "../execution_plan/ops/op.h"
typedef struct Timer_mql
{
	OpBase *p;
	double time_sum;
    int record_sum;
	uint64_t * distinct_records;
}Timer_mql;
void simpletimer_start_mql();
double simpletimer_end_mql(); 
Timer_mql** get_timers_mql();
void timers_init_mql();
void timers_append_mql(OpBase *p);
void add_to_timer_mql(OpBase *p,double addtime,int addrecord,int adddistinctrecord);
void timers_output_CondTraverse_mql();