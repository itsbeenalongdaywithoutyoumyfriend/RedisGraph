#include "../graph/entities/qg_edge.h"
typedef struct PathRecorder
{
    char*** paths;
}PathRecorder;
PathRecorder* get_pathrecorder_mql();
void pathrecorder_init_mql();
bool add_to_pathrecorder_mql(QGEdge **path);//if exist return 1 else return 0 then add