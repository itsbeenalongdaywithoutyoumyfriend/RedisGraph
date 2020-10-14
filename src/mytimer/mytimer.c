#include "mytimer.h"
#include "../execution_plan/ops/op.h"
#include "../execution_plan/ops/op_conditional_traverse.h"
#include<time.h>

Timer_mql** get_timers_mql()
{
	static Timer_mql* timers;
	return &timers;
}
void timers_init_mql()
{
	Timer_mql **timers_p=get_timers_mql();
	*timers_p=array_new(Timer_mql,64);
}
void timers_append_mql(OpBase *p)
{
    Timer_mql* timers = *get_timers_mql();
    int len=array_len(timers);
	int ifexist=0;
    for(int i=0;i<len;++i)
	{
		if(timers[i].p==p)
		{
            ifexist=1;
        }
    }
    if(ifexist)return;
    Timer_mql* newrecord=rm_malloc(sizeof(Timer_mql));
    newrecord->time_sum=0;
    newrecord->record_sum=0;
    newrecord->p=p;
    timers=array_append(timers,*newrecord);
}

void add_to_timer_mql(OpBase *p,double addtime,int addrecord)
{
	Timer_mql* timers = *get_timers_mql();
	int len=array_len(timers);
	int ifexist=0;
	for(int i=0;i<len;++i)
	{
		if(timers[i].p==p)
		{
			timers[i].time_sum+=addtime;
            timers[i].record_sum+=addrecord;
			ifexist=1;
			break;
		}
	}
    assert(ifexist);
}

void timers_output_CondTraverse_mql()
{
    FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
    Timer_mql* timers = *get_timers_mql();
    int len=array_len(timers);
    for(int i=0;i<len;++i)
	{
        CondTraverse *op = (CondTraverse *)timers[i].p;
        fprintf(fp,"%s %s CondTraverse time_used:%lfms,record_scanned:%d\n",
            op->ae->operand.src,
            op->ae->operand.dest,
            timers[i].time_sum,
            timers[i].record_sum
            );
    }
    fclose(fp);
        
}

clock_t simpletimer_cal_mql()
{
	static clock_t oldclock=0;
	clock_t newclock=clock();
	clock_t d=newclock-oldclock;
	oldclock=newclock;
	return d;
}
void simpletimer_start_mql()
{
	simpletimer_cal_mql();
}
double simpletimer_end_mql() // ms
{
	clock_t d=simpletimer_cal_mql();
	return d*1.0/CLOCKS_PER_SEC*1000;
}