#include "mytimer.h"
#include<time.h>
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
int simpletimer_end_mql() // ms
{
	clock_t d=simpletimer_cal_mql();
	return d;
}