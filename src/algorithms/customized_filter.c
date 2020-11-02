
#include "../arithmetic/algebraic_expression.h"
#include "dfs.h"

void customized_filter_on_cycle_mql(QueryGraph *qg)
{
	DFS_mql(qg);
}


void build_customized_filter_on_cycle_mql(QGNode *n, int path_len, QGEdge ***path, bool *transpositions, QueryGraph *qg)
{
	assert(array_len(*path)==path_len);
	uint firstPathIndex;
	for(firstPathIndex=0;firstPathIndex<path_len;++firstPathIndex)
	{
		if(transpositions[firstPathIndex]==false&&(*path)[firstPathIndex]->src==n)break;
		if(transpositions[firstPathIndex]==true&&(*path)[firstPathIndex]->dest==n)break;
	}
	assert(firstPathIndex<path_len);
	uint part_path_len=path_len-firstPathIndex;
	QGEdge **part_path = array_new(QGEdge *, part_path_len);
	bool part_transpositions[part_path_len];
	uint transposeCount = 0;
	for(uint i=0,j=firstPathIndex;j<path_len;++j,++i)
	{
		part_path=array_append(part_path,(*path)[j]);
		part_transpositions[i]=transpositions[j];
		if(part_transpositions[i])++transposeCount;
	}
	if(transposeCount > (part_path_len - transposeCount)) {
		_reversePath(part_path, part_path_len, part_transpositions);
	}
	// Apply transpose.
	for(uint i = 0; i < part_path_len; i++) {
		QGEdge *e = part_path[i];
		if(part_transpositions[i]) QGEdge_Reverse(e);
	}

	for(uint i=0;i<part_path_len;++i)
	{
		QGEdge **rotated_path = array_new(QGEdge *, part_path_len);
		bool rotated_transpositions[part_path_len];
		for(uint j=0;j<part_path_len;++j)
		{
			rotated_path=array_append(rotated_path,part_path[(j+i)%part_path_len]);
			rotated_transpositions[j]=part_transpositions[(j+i)%part_path_len];
		}
		NodeID *filters = get_filter_mql_on_cycle_mql(rotated_path,rotated_transpositions);
		fill_customized_filter_mql(part_path[i]->src->customized_filter,filters);
	}

	//undo  transpose.
	for(uint i = 0; i < part_path_len; i++) {
		QGEdge *e = part_path[i];
		if(part_transpositions[i]) QGEdge_Reverse(e);
	}


}
