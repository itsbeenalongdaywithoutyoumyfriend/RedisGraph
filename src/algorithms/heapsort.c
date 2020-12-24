#include "../graph/entities/graph_entity.h"
void swap_mql(NodeID* a, NodeID* b)
{
    NodeID temp = *b;
    *b = *a;
    *a = temp;
}

void max_heapify_mql(NodeID arr[], int start, int end) 
{
    //建立父节点指标和子节点指标
    int dad = start;
    int son = dad * 2 + 1;
    while (son <= end)  //若子节点指标在范围内才做比较
        {
            if (son + 1 <= end && arr[son] < arr[son + 1]) 
            //先比较两个子节点大小，选择最大的
            son++;
        if (arr[dad] > arr[son]) //如果父节点大於子节点代表调整完毕，直接跳出函数
            return;
        else  //否则交换父子内容再继续子节点和孙节点比较
        {
            swap_mql(&arr[dad], &arr[son]);
            dad = son;
            son = dad * 2 + 1;
        }
    }
}
 
void heap_sort_mql(NodeID arr[], int len) 
{
    int i;
    //初始化，i从最後一个父节点开始调整
    for (i = len / 2 - 1; i >= 0; i--)
        max_heapify_mql(arr, i, len - 1);
    //先将第一个元素和已排好元素前一位做交换，再重新调整，直到排序完毕
    for (i = len - 1; i > 0; i--) 
    {
        swap_mql(&arr[0], &arr[i]);
        max_heapify_mql(arr, 0, i - 1);
    }
}