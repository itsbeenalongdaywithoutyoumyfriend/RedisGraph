#include "path_recorder.h"
#include "../util/arr.h"
#include<string.h>
PathRecorder* get_pathrecorder_mql()
{
    static PathRecorder pr;
    return &pr;
}
void pathrecorder_init_mql()
{
    PathRecorder* p = get_pathrecorder_mql();
    p->paths=array_new(char**,0);
}
bool add_to_pathrecorder_mql(QGEdge **path)
{
    int plen=array_len(path);
    char** p=array_new(char*,plen);
    for(int i=0;i<plen;++i)p=array_append(p,path[i]->alias);
    bool ifexist=0;
    PathRecorder* pr=get_pathrecorder_mql();
    int prlen=array_len(pr->paths);
    for(int i=0;i<prlen;++i)
    {
        int thislen=array_len(pr->paths[i]);
        if(thislen==plen)
        {
            int cnt=0;
            for(int j=0;j<plen;++j)
            {
                for(int k=0;k<plen;++k)
                {
                    if(strcmp(p[j],pr->paths[i][k])==0)
                    {
                        ++cnt;
                        break;
                    }
                }
            }
            if(cnt==plen)
            {
                ifexist=1;
                break;
            }
        }
    }
    if(ifexist)return true;
    pr->paths=array_append(pr->paths,p);
}