#include "mytimer.h"
#include "../execution_plan/ops/op.h"
#include "../execution_plan/ops/op_conditional_traverse.h"
#include <sys/time.h>
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
        fprintf(fp,"%d -> %d CondTraverse time_used:%lfms,record_scanned:%d\n",
            op->srcNodeIdx,
            op->destNodeIdx,
            timers[i].time_sum,
            timers[i].record_sum
            );
    }
    fclose(fp);
        
}

double simpletimer_cal_mql()
{
	static double oldclock=0;
	struct timeval nowtime;
    gettimeofday(&nowtime,0);
    double newclock=1000000*nowtime.tv_sec+nowtime.tv_usec;
	double d=newclock-oldclock;
	oldclock=newclock;
	return d/1000.0;
}
void simpletimer_start_mql()
{
	simpletimer_cal_mql();
}
double simpletimer_end_mql() // ms
{
	double d=simpletimer_cal_mql();
	return d;
}