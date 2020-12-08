/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/
#include "../algebraic_expression.h"
#include "utils.h"
#include "../../util/arr.h"
#include "../../query_ctx.h"
#include "../../algorithms/algorithms.h"
#include "../../mytimer/mytimer.h"

/* Node with (income + outcome degree) > 2
 * is considered a highly connected node. */
static bool _highly_connected_node(const QueryGraph *qg, const char *alias) {
	// Look up node in qg.
	QGNode *n = QueryGraph_GetNodeByAlias(qg, alias);
	return n->highly_connected;
}

static inline bool _referred_entity(const char *alias) {
	AST *ast = QueryCtx_GetAST();
	return AST_AliasIsReferenced(ast, alias);
}

/* If the edge is referenced or of a variable length, it should populate the AlgebraicExpression. */
static inline bool _should_populate_edge(QGEdge *e) {
	return (_referred_entity(e->alias) || QGEdge_VariableLength(e));
}

/* Checks if given expression contains a variable length edge. */
static bool _AlgebraicExpression_ContainsVariableLengthEdge
(
	const QueryGraph *qg,
	const AlgebraicExpression *root
) {
	uint child_count = 0;
	switch(root->type) {
	case AL_OPERATION:
		child_count = AlgebraicExpression_ChildCount(root);
		for(uint i = 0; i < child_count; i++) {
			if(_AlgebraicExpression_ContainsVariableLengthEdge(qg, root->operation.children[i])) return true;
		}
		return false;
	case AL_OPERAND:
		if(root->operand.edge) {
			QGEdge *e = QueryGraph_GetEdgeByAlias(qg, root->operand.edge);
			return QGEdge_VariableLength(e);
		}
		break;
	default:
		assert("Unknow algebraic expression node type" && false);
	}
	return false;
}

static void _RemovePathFromGraph(QueryGraph *g, QGEdge **path) {
	uint edge_count = array_len(path);
	for(uint i = 0; i < edge_count; i++) {
		QGEdge *e = path[i];
		QGNode *src = e->src;
		QGNode *dest = e->dest;

		QueryGraph_RemoveEdge(g, e);
		QGEdge_Free(e);
		if(QGNode_EdgeCount(src) == 0) {
			QueryGraph_RemoveNode(g, src);
			QGNode_Free(src);
		}
		if(src != dest && QGNode_EdgeCount(dest) == 0) {
			QueryGraph_RemoveNode(g, dest);
			QGNode_Free(dest);
		}
	}
}

static inline bool _should_divide_expression(QGEdge **path, int idx, const QueryGraph *qg) {
	QGEdge *e = path[idx];

	return (_should_populate_edge(e)                    ||  // This edge is populated.
			_should_populate_edge(path[idx + 1])        ||  // The next edge is populated.
			_highly_connected_node(qg, e->dest->alias)  ||  // Destination node in+out degree > 2.
			_referred_entity(e->dest->alias));              // Destination node is referenced.
}

/* Variable length expression must contain only a single operand, the edge being
 * traversed multiple times, in cases such as (:labelA)-[e*]->(:labelB) both label A and B
 * are applied via a label matrix operand, this function migrates A and B from a
 * variable length expression to other expressions. */
static AlgebraicExpression **_AlgebraicExpression_IsolateVariableLenExps(
	const QueryGraph *qg, AlgebraicExpression **expressions) {
	/* Return value is a new set of expressions, where each variable length expression
	 * is guaranteed to have a single operand, as such in the worst case the number of
	 * expressions doubles + 1. */
	size_t expCount = array_len(expressions);
	AlgebraicExpression **res = array_new(AlgebraicExpression *, expCount * 2 + 1);

	/* Scan through each expression, locate expression which
	 * have a variable length edge in them. */
	for(size_t expIdx = 0; expIdx < expCount; expIdx++) {
		AlgebraicExpression *exp = expressions[expIdx];
		if(!_AlgebraicExpression_ContainsVariableLengthEdge(qg, exp)) {
			res = array_append(res, exp);
			continue;
		}

		// Expression contains a variable length edge.
		QGNode *src = QueryGraph_GetNodeByAlias(qg, AlgebraicExpression_Source(exp));

		// A variable length expression with a labeled source node
		// We only care about the source label matrix, when it comes to
		// the first expression, as in the following expressions
		// src is the destination of the previous expression.
		if(expIdx == 0 && src->label) {
			// Remove src node matrix from expression.
			AlgebraicExpression *op = AlgebraicExpression_RemoveSource(&exp);
			res = array_append(res, op);
		}

		res = array_append(res, exp);

		// If the expression has a labeled destination, separate it into its own expression.
		QGNode *dest = QueryGraph_GetNodeByAlias(qg, AlgebraicExpression_Destination(exp));
		if(dest->label) {
			// Remove dest node matrix from expression.
			AlgebraicExpression *op = AlgebraicExpression_RemoveDest(&exp);

			/* See if dest mat can be prepended to the following expression.
			 * If not create a new expression. */
			if(expIdx < expCount - 1 &&
			   !_AlgebraicExpression_ContainsVariableLengthEdge(qg, expressions[expIdx + 1])) {
				expressions[expIdx + 1] = _AlgebraicExpression_MultiplyToTheLeft(op, expressions[expIdx + 1]);
			} else {
				res = array_append(res, op);
			}
		}
	}

	array_free(expressions);
	return res;
}

/* Break down path into sub paths.
 * Considering referenced intermidate nodes and edges. */
static QGEdge ***_Intermediate_Paths
(
	QGEdge **path,
	const QueryGraph *qg
) {
	QGEdge *e = NULL;
	int pathLen = array_len(path);

	/* Allocating maximum number of expression possible. */
	QGEdge ***paths = array_new(QGEdge **, pathLen);
	QGEdge **intermediate_path = array_new(QGEdge *, pathLen);
	paths = array_append(paths, intermediate_path);

	/* Scan path left to right,
	 * construct intermidate paths by "breaking" on referenced entities. */
	for(int i = 0; i < pathLen - 1; i++) {
		e = path[i];
		intermediate_path = array_append(intermediate_path, e);
		if(_should_divide_expression(path, i, qg)) {
			// Break! add current path to paths and create a new path.
			intermediate_path = array_new(QGEdge *, pathLen);
			paths = array_append(paths, intermediate_path);
		}
	}

	// Handle last hop.
	e = path[pathLen - 1];
	intermediate_path = array_append(intermediate_path, e);

	return paths;
}

static AlgebraicExpression *_AlgebraicExpression_OperandFromNode
(
	QGNode *n
) {
	bool diagonal = true;
	bool transpose = false;
	return AlgebraicExpression_NewOperand(GrB_NULL, diagonal, n->alias, n->alias, NULL, n->label);
}

static AlgebraicExpression *_AlgebraicExpression_OperandFromNodeFilter_mql
(
	QGNode *n,
	NodeID **p
) {
	bool diagonal = true;
	bool transpose = false;
	AlgebraicExpression *ret=AlgebraicExpression_NewOperand(GrB_NULL, diagonal, n->alias, n->alias, NULL, n->label);
	ret->operand.customized_filter_pointer=p;
	return ret;
}


static AlgebraicExpression *_AlgebraicExpression_OperandFromEdge
(
	QGEdge *e,
	bool transpose
) {
	GrB_Matrix mat;
	uint reltype_id;
	Graph *g = QueryCtx_GetGraph();
	AlgebraicExpression *add = NULL;
	AlgebraicExpression *root = NULL;
	AlgebraicExpression *src_filter = NULL;

	QGNode *src_node = e->src;
	QGNode *dest_node = e->dest;

	// Use original `src` and `dest` for algebraic operands.
	const char *src = (transpose) ? dest_node->alias : src_node->alias;
	const char *dest = (transpose) ? src_node->alias : dest_node->alias;
	const char *edge = _should_populate_edge(e) ? e->alias : NULL;
	bool var_len_traversal = QGEdge_VariableLength(e);

	// If src node has a label, multiply to the left by label matrix.
	if(src_node->label) {
		src_filter = _AlgebraicExpression_OperandFromNode(src_node);
	}

	/* No hops: (a)-[:R*0]->(b)
	 * in this case we want to use the identity matrix
	 * f * I  = f */
	if(!var_len_traversal && e->minHops == 0) {
		root = AlgebraicExpression_NewOperand(IDENTITY_MATRIX, true, src, dest, edge, "I");
	} else {
		uint reltype_count = array_len(e->reltypeIDs);
		switch(reltype_count) {
		case 0: // No relationship types specified; use the full adjacency matrix
			root = AlgebraicExpression_NewOperand(GrB_NULL, false, src, dest, edge, NULL);
			break;
		case 1: // One relationship type
			root = AlgebraicExpression_NewOperand(GrB_NULL, false, src, dest, edge, e->reltypes[0]);
			break;
		default: // Multiple edge type: -[:A|:B]->
			add = AlgebraicExpression_NewOperation(AL_EXP_ADD);
			for(uint i = 0; i < reltype_count; i++) {
				AlgebraicExpression *operand = AlgebraicExpression_NewOperand(GrB_NULL, false, src, dest,
																			  edge, e->reltypes[i]);
				AlgebraicExpression_AddChild(add, operand);
			}
			root = add;
			break;
		}

		if(e->bidirectional) {
			/* ()-[]-()
			 * The Adj + Transpose(The Adj)
			 *
			 * ()-[:R]-()
			 * R + Transpose(R)
			 *
			 * ()-[:R0|R1]-()
			 * (R0 + R1) + Transpose(R0 + R1) */
			add = AlgebraicExpression_NewOperation(AL_EXP_ADD);
			AlgebraicExpression_AddChild(add, root);

			AlgebraicExpression *op_transpose = AlgebraicExpression_NewOperation(AL_EXP_TRANSPOSE);
			AlgebraicExpression_AddChild(op_transpose, AlgebraicExpression_Clone(root));
			AlgebraicExpression_AddChild(add, op_transpose);
			root = add;
		}

		/* Expand fixed variable length edge.
		 * -[A*2..2]->
		 * A*A
		 * -[A|B*2..2]->
		 * (A+B) * (A+B) */
		if(!var_len_traversal && e->minHops > 1) {
			AlgebraicExpression *mul = AlgebraicExpression_NewOperation(AL_EXP_MUL);
			AlgebraicExpression_AddChild(mul, root);
			for(int i = 1; i < e->minHops; i++) {
				// Clone to avoid double free.
				AlgebraicExpression_AddChild(mul, AlgebraicExpression_Clone(root));
			}
			root = mul;
		}
	}

	// Transpose entire expression.
	if(transpose) {
		AlgebraicExpression *op_transpose = AlgebraicExpression_NewOperation(AL_EXP_TRANSPOSE);
		AlgebraicExpression_AddChild(op_transpose, root);
		root = op_transpose;
	}

	// Apply source filter by multiplying to the left.
	if(src_filter) {
		root = _AlgebraicExpression_MultiplyToTheLeft(src_filter, root);
	}

	return root;
}

/* In case edges `a` and `b` share a node:
 * (a)-[E0]->(b)<-[E1]-(c)
 * than the shared entity is returned
 * if edges are disjoint, NULL is returned. */
static QGNode *_SharedNode
(
	const QGEdge *a,
	const QGEdge *b
) {
	assert(a && b);
	if(a->dest == b->src) return a->dest;   // (a)-[E0]->(b)-[E1]->(c)
	if(a->src == b->dest) return a->src;    // (a)<-[E0]-(b)<-[E1]-(c)
	if(a->src == b->src) return a->src;     // (a)<-[E0]-(b)-[E1]->(c)
	if(a->dest == b->dest) return a->dest;  // (a)-[E0]->(b)<-[E1]-(c)
	return NULL;
}

static void _reversePath
(
	QGEdge **path,
	uint path_len,
	bool *transpositions
) {
	for(uint i = 0; i < path_len; i++) {
		/* A reversed path should have its transpositions flipped
		 * transpose(transpose(A)) = A */
		transpositions[i] = !transpositions[i];
	}

	/* Transpose(A*B) = Transpose(B) * Transpose(A)
	 * (a)<-[A]-(b)<-[B]-(c)-[C]->(d)
	 * At * Bt * C
	 * Transpose(At * Bt * C) =
	 * = Ct * B * A
	 * (d)-[Ct]->(c)-[B]->(b)-[A]->(a) */

	// Reverse the path array as well as the transposition array
	for(uint i = 0; i < path_len / 2; i++) {
		uint opposite = path_len - i - 1;
		QGEdge *tmp = path[opposite];
		path[opposite] = path[i];
		path[i] = tmp;

		bool transpose_tmp = transpositions[opposite];
		transpositions[opposite] = transpositions[i];
		transpositions[i] = transpose_tmp;
	}
}

static void _normalizePath
(
	QGEdge **path,          // Path to normalize.
	uint path_len,          // Path length.
	bool *transpositions    // Specifies which edges need to be transposed.
) {
	// Initialize `transpositions` array.
	for(uint i = 0; i < path_len; i++) transpositions[i] = false;

	// A single leg path.
	if(path_len <= 1) return;

	uint transposeCount = 0;
	// For every edge except the last:
	for(uint i = 0; i < path_len - 1; i++) {
		QGEdge *e = path[i];
		QGEdge *follow = path[i + 1];
		QGNode *shared = _SharedNode(e, follow);
		assert(shared);

		/* The edge should be transposed if its destination is not shared.
		 * (dest)<-[e]-(shared)-[follow]->()
		 * (dest)<-[e]-(shared)<-[follow]-() */
		if(e->dest != shared) {
			transpositions[i] = true;
			transposeCount++;
		}
	}

	// For the last edge, transpose if its source is not shared.
	QGEdge *e = path[path_len - 1];
	QGNode *shared = _SharedNode(path[path_len - 2], e);
	if(e->src != shared) {
		transposeCount++;
		transpositions[path_len - 1] = true;
	}

	// Reverse entire path if the majority of edges must be transposed.
	if(transposeCount > (path_len - transposeCount)) {
		_reversePath(path, path_len, transpositions);
	}

	// Apply transpose.
	for(uint i = 0; i < path_len; i++) {
		QGEdge *e = path[i];
		if(transpositions[i]) QGEdge_Reverse(e);
	}
}

static AlgebraicExpression *_AlgebraicExpression_FromPath
(
	QGEdge **path,
	bool *transpositions
) {
	assert(path);

	QGEdge *e = NULL;
	uint path_len = array_len(path);
	assert(path_len > 0);
	AlgebraicExpression *root = NULL;

	/* Treating path as a chain
	 * we're aligning all edges to "point right"
	 * (A)-[E0]->(B)-[E0]->(C)-[E0]->(D).
	 * e.g.
	 * (A)-[E0]->(B)<-[E1]-(C)-[E2]->(D)
	 * E1 will be transposed:
	 * (A)-[E0]->(B)-[E1]->(C)-[E2]->(D) */

	// Construct expression.
	for(int i = 0; i < path_len; i++) {
		e = path[i];
		// Add Edge matrix.
		AlgebraicExpression *op = _AlgebraicExpression_OperandFromEdge(e, transpositions[i]);

		if(!root) {
			root = op;
		} else {
			// Connect via a multiplication node.
			root = _AlgebraicExpression_MultiplyToTheRight(root, op);
		}
	}   // End of path traversal.

	// If last node on path has a label, multiply by label matrix.
	if(e->dest->label) {
		root = _AlgebraicExpression_MultiplyToTheRight(root, _AlgebraicExpression_OperandFromNode(e->dest));
	}

	if((e->dest->customized_filter)!=NULL){
		root = _AlgebraicExpression_MultiplyToTheRight(root, _AlgebraicExpression_OperandFromNodeFilter_mql(e->dest,&(e->dest->customized_filter)));
	}
	
	return root;
}

//------------------------------------------------------------------------------
// AlgebraicExpression construction.
//------------------------------------------------------------------------------

static AlgebraicExpression *_AlgebraicExpression_FromPath_mql
(
	QGEdge **path,
	bool *transpositions
) {
	assert(path);

	QGEdge *e = NULL;
	uint path_len = array_len(path);
	assert(path_len > 0);
	AlgebraicExpression *root = NULL;

	/* Treating path as a chain
	 * we're aligning all edges to "point right"
	 * (A)-[E0]->(B)-[E0]->(C)-[E0]->(D).
	 * e.g.
	 * (A)-[E0]->(B)<-[E1]-(C)-[E2]->(D)
	 * E1 will be transposed:
	 * (A)-[E0]->(B)-[E1]->(C)-[E2]->(D) */

	// Construct expression.
	for(int i = 0; i < path_len; i++) {
		e = path[i];
		// Add Edge matrix.
		AlgebraicExpression *op = _AlgebraicExpression_OperandFromEdge(e, transpositions[i]);

		if(!root) {
			root = op;
		} else {
			// Connect via a multiplication node.
			root = _AlgebraicExpression_MultiplyToTheRight(root, op);
		}
	}   // End of path traversal.

	// If last node on path has a label, multiply by label matrix.
	if(e->dest->label) {
		root = _AlgebraicExpression_MultiplyToTheRight(root, _AlgebraicExpression_OperandFromNode(e->dest));
	}
	return root;
}

NodeID * get_filter_on_cycle_mql
(
	QGEdge **path,
	bool *transpositions
)
{
	const int recordsCap=16;
	GraphContext *gc = QueryCtx_GetGraphCtx();
	AlgebraicExpression *exp=_AlgebraicExpression_FromPath_mql(path,transpositions);
	size_t required_dim = Graph_RequiredMatrixDim(gc->g);
	NodeID *filters = array_new(NodeID,required_dim);
	GrB_Matrix res= GrB_NULL;
	GrB_Matrix recordsBulk=GrB_NULL;
	GrB_Matrix_new(&recordsBulk, GrB_BOOL, recordsCap, required_dim);
	GrB_Matrix_new(&res, GrB_BOOL, recordsCap, required_dim);
	AlgebraicExpression_MultiplyToTheLeft(&exp, recordsBulk);
	AlgebraicExpression_Optimize(&exp);
	assert(exp);
	assert(exp->type == AL_OPERATION);
	NodeID* original_filter=path[0]->src->customized_filter;
	if(original_filter==NULL)
	{
		original_filter=array_new(NodeID,required_dim);
		for(int i=0;i<required_dim;++i)
			original_filter=array_append(original_filter,i);
	}
	int original_filter_len=array_len(original_filter);
	for(int i=0;i<original_filter_len;i+=recordsCap)
	{
		for(int j=0;j<recordsCap&&i+j<original_filter_len;++j)
		{
			GrB_Matrix_setElement_BOOL(recordsBulk, true, j, original_filter[i+j]);
		}
		AlgebraicExpression_Eval(exp, res);
		GrB_Matrix_clear(recordsBulk);
		bool v;
		for(int j=0;j<recordsCap&&i+j<original_filter_len;++j)
		{
			if(GrB_Matrix_extractElement_BOOL(&v,res,j,original_filter[i+j])==GrB_SUCCESS)
				filters= array_append(filters,original_filter[i+j]);
		}
	}
	GrB_Matrix_free(&res);
	GrB_Matrix_free(&recordsBulk);
	FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	fprintf(fp,"getfiltercycle %d\n",array_len(filters));
	fclose(fp);
	return filters;
}

NodeID * get_filter_from_starters_mql
(
	QGEdge **path,
	bool *transpositions,
	NodeID *starters
)
{
	const int recordsCap=16;
	GraphContext *gc = QueryCtx_GetGraphCtx();
	AlgebraicExpression *exp=_AlgebraicExpression_FromPath_mql(path,transpositions);
	size_t required_dim = Graph_RequiredMatrixDim(gc->g);
	NodeID *filters = array_new(NodeID,required_dim);
	GrB_Matrix res= GrB_NULL;
	GrB_Matrix recordsBulk=GrB_NULL;
	GrB_Matrix_new(&recordsBulk, GrB_BOOL, recordsCap, required_dim);
	GrB_Matrix_new(&res, GrB_BOOL, recordsCap, required_dim);
	AlgebraicExpression_MultiplyToTheLeft(&exp, recordsBulk);
	AlgebraicExpression_Optimize(&exp);
	assert(exp);
	assert(exp->type == AL_OPERATION);
	int starters_len=array_len(starters);
	for(int i=0;i<starters_len;i+=recordsCap)
	{
		for(int j=0;j<recordsCap&&i+j<starters_len;++j)
		{
			GrB_Matrix_setElement_BOOL(recordsBulk, true, j, starters[i+j]);
		}
		AlgebraicExpression_Eval(exp, res);
		GrB_Matrix_clear(recordsBulk);
		bool v;
		GxB_MatrixTupleIter *iter=NULL;
		GxB_MatrixTupleIter_new(&iter, res);
		NodeID src_id = INVALID_ENTITY_ID;
		NodeID dest_id = INVALID_ENTITY_ID;
		bool depleted = false;
		while(true)
		{
			if(iter) GxB_MatrixTupleIter_next(iter, &src_id, &dest_id, &depleted);
			if(depleted) break;
			filters= array_append(filters,dest_id);
		}
	}
	GrB_Matrix_free(&res);
	GrB_Matrix_free(&recordsBulk);
	int filters_len=array_len(filters);
	heap_sort_mql(filters,filters_len);
	NodeID *return_filters = array_new(NodeID,required_dim);
	for(int i=0;i<filters_len;++i)
	{
		if(i==filters_len-1||filters[i]!=filters[i+1])
			return_filters=array_append(return_filters,filters[i]);
	}
	array_free(filters);
	return return_filters;
}

NodeID * get_filter_mql
(
	QGEdge **path,
	bool *transpositions,
	bool src_or_dest // src 1 dest 0
)
{
	GraphContext *gc = QueryCtx_GetGraphCtx();
	AlgebraicExpression *exp=_AlgebraicExpression_FromPath_mql(path,transpositions);
	GrB_Matrix res= GrB_NULL;
	size_t required_dim = Graph_RequiredMatrixDim(gc->g);
	NodeID *filters = array_new(NodeID,required_dim);
	GrB_Matrix_new(&res, GrB_BOOL, required_dim, required_dim);
	AlgebraicExpression_Optimize(&exp);
	assert(exp);
	if(exp->type == AL_OPERATION){
		AlgebraicExpression_Eval(exp, res);
	}
	else
	{
		assert(exp->type == AL_OPERAND);
		res = exp->operand.matrix;
	}
	
	GxB_MatrixTupleIter *iter=NULL;
	GxB_MatrixTupleIter_new(&iter, res);
	NodeID src_id = INVALID_ENTITY_ID;
	NodeID dest_id = INVALID_ENTITY_ID;
	bool depleted = false;
	while(true)
	{
		if(iter) GxB_MatrixTupleIter_next(iter, &src_id, &dest_id, &depleted);
		if(depleted) break;
		if(src_or_dest)filters= array_append(filters,src_id);
		else filters= array_append(filters,dest_id);
	}
	FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	fprintf(fp,"get_filter_mql %d\n",array_len(filters));
	fclose(fp);
	return filters;
}


void fill_customized_filter_mql
(
	NodeID **to_be_filled,
	NodeID *filter_array // must be sorted
){//will free used filter_array and original to_be_filled
	GraphContext *gc = QueryCtx_GetGraphCtx();
	size_t required_dim = Graph_RequiredMatrixDim(gc->g);
	uint filter_len=array_len(filter_array);
	int cnt=0;
	if(*to_be_filled==NULL)
	{
		NodeID *new_filter_array = array_new(NodeID, required_dim);
		for(int i=0;i<filter_len;++i)
		{
			if(i>0&&filter_array[i]==filter_array[i-1])continue;
			new_filter_array=array_append(new_filter_array,filter_array[i]);
		}
		*to_be_filled = new_filter_array;
		cnt=array_len(new_filter_array);
		array_free(filter_array);
	}
	else
	{
		NodeID *original_filter_array=*to_be_filled;
		uint original_filter_len = array_len(original_filter_array);
		NodeID *new_filter_array = array_new(NodeID, required_dim);
		for(uint i=0,j=0;j<original_filter_len;++j)
		{
			if(j>0&&original_filter_array[j]==original_filter_array[j-1])continue;
			while(i<filter_len&&filter_array[i]<original_filter_array[j])++i;
			if(i<filter_len&&filter_array[i]==original_filter_array[j])
			{
				new_filter_array = array_append(new_filter_array,filter_array[i]);
			}
		}
		*to_be_filled = new_filter_array;
		array_free(original_filter_array);
		array_free(filter_array);
		cnt=array_len(new_filter_array);
	}
	FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	fprintf(fp,"%d\n",cnt);
	fclose(fp);
}




void customized_filter_mql
(
	QGEdge **path,
	bool *transpositions,
	const QueryGraph *qg
){
	simpletimer_start_mql();
	QGEdge *e = NULL;
	GrB_Index nvals;
	int pathLen = array_len(path);
	
	/* Scan path left to right,
	 * construct intermidate paths by "breaking" on referenced entities. */
	for(int i = 0; i < pathLen - 1; i++) {
		// e = path[i];
		// intermediate_path = array_append(intermediate_path, e);
		if(_should_divide_expression(path, i, qg)) {
			// Break! add current path to paths and create a new path.
			// intermediate_path = array_new(QGEdge *, pathLen);
			// paths = array_append(paths, intermediate_path);
			QGEdge **path1 = array_new(QGEdge *, pathLen);
			QGEdge **path2 = array_new(QGEdge *, pathLen);
			for(int j=0;j<=i;++j)
			{
				e=path[j];
				path1=array_append(path1, e);
			}
			for(int j=i+1;j<pathLen;++j)
			{
				e=path[j];
				path2=array_append(path2, e);
			}
			e = path[i];

			
			

			NodeID *filters1 = get_filter_mql(path1,transpositions,0);
			uint edge_converted = array_len(path1);
			NodeID *filters2 = get_filter_mql(path2,transpositions + edge_converted,1);
			uint filters1_len=array_len(filters1);
			uint filters2_len=array_len(filters2);
			heap_sort_mql(filters1,filters1_len);

			fill_customized_filter_mql(&e->dest->customized_filter,filters2);
			fill_customized_filter_mql(&e->dest->customized_filter,filters1);
			FILE *fp;
			fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
			fprintf(fp,"%s %llu\n",e->dest->alias,array_len(e->dest->customized_filter));
			fclose(fp);

			array_free(path1);
			array_free(path2);
			array_free(filters1);
			array_free(filters2);
		}
	}
	NodeID *src_filter=get_filter_mql(path,transpositions,1);
	NodeID *dest_filter=get_filter_mql(path,transpositions,0);
	src_filter=get_filter_mql(path,transpositions,1);
	dest_filter=get_filter_mql(path,transpositions,0);
	uint dest_filter_len=array_len(dest_filter);
	heap_sort_mql(dest_filter,dest_filter_len);

	fill_customized_filter_mql(&path[0]->src->customized_filter,src_filter);
	fill_customized_filter_mql(&path[pathLen-1]->dest->customized_filter,dest_filter);


	FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	fprintf(fp,"%s %llu src\n",path[0]->src->alias,array_len(path[0]->src->customized_filter));
	fprintf(fp,"%s %llu dest\n",path[pathLen-1]->dest->alias,array_len(path[pathLen-1]->dest->customized_filter));

	
	double time_used=simpletimer_end_mql();
	fprintf(fp,"%lfms used in customized_filter_mql\n",time_used);


	fclose(fp);
	array_free(src_filter);
	array_free(dest_filter);

}
void customized_filter_on_cycle_mql(QueryGraph *qg)
{
	DFS_mql(qg);
}


void build_customized_filter_on_cycle_mql(QGNode *n, int path_len, QGEdge ***path, bool *transpositions, QueryGraph *qg)
{

	simpletimer_start_mql();
	FILE *fp;
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	assert(array_len(*path)==path_len);
	uint firstPathIndex;
	for(firstPathIndex=0;firstPathIndex<path_len;++firstPathIndex)
	{
		fprintf(fp,"%s->%s,",(*path)[firstPathIndex]->src->alias,(*path)[firstPathIndex]->dest->alias);
	}fprintf(fp,"\n");
	fclose(fp);
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
	if(add_to_pathrecorder_mql(part_path))return;
	// Apply transpose.
	for(uint i = 0; i < part_path_len; i++) {
		QGEdge *e = part_path[i];
		if(part_transpositions[i]) QGEdge_Reverse(e);
	}
	int last_filter_index=-1;
	for(uint i=0;i<part_path_len;++i)
	{
		QGEdge **rotated_path = array_new(QGEdge *, part_path_len);
		bool rotated_transpositions[part_path_len];
		for(uint j=0;j<part_path_len;++j)
		{
			rotated_path=array_append(rotated_path,part_path[(j+i)%part_path_len]);
			rotated_transpositions[j]=part_transpositions[(j+i)%part_path_len];
		}
		if(part_path_len>1&&!_referred_entity(rotated_path[0]->src->alias))continue;
		//_should_divide_expression func not applicable
		if(rotated_path[0]->src->customized_filter==NULL&&last_filter_index>0)
		{
			QGEdge **starters_path = array_new(QGEdge *, part_path_len);
			bool starters_transpositions[part_path_len];
			for(int j=0;last_filter_index+j<i;++j)
			{
				starters_path=array_append(starters_path,part_path[last_filter_index+j]);
				starters_transpositions[j]=part_transpositions[last_filter_index+j];
			}
			NodeID *filters = get_filter_from_starters_mql(starters_path,starters_transpositions,part_path[last_filter_index]->src->customized_filter);
			fill_customized_filter_mql(&part_path[i]->src->customized_filter,filters);
			array_free(starters_path);
		}
		NodeID *filters = get_filter_on_cycle_mql(rotated_path,rotated_transpositions);
		fill_customized_filter_mql(&part_path[i]->src->customized_filter,filters);
		last_filter_index=i;
		array_free(rotated_path);
	}

	//undo  transpose.
	for(uint i = 0; i < part_path_len; i++) {
		QGEdge *e = part_path[i];
		if(part_transpositions[i]) QGEdge_Reverse(e);
	}
	fp=fopen("/home/qlma/customized-filter/outcount-redisgraph-mql","a+");
	double time_used=simpletimer_end_mql();
	fprintf(fp,"%lfms used in build_customized_filter_on_cycle_mql\n",time_used);
	fclose(fp);
}
// Construct algebraic expression form query graph.
AlgebraicExpression **AlgebraicExpression_FromQueryGraph
(
	QueryGraph *qg    // Query-graph to process
) {
	assert(qg);

	/* Construct algebraic expression(s) from query-graph.
	 * Trying to take advantage of long multiplications with as few
	 * transpose as possible we'll transform paths crossing the graph
	 * "diameter", these are guarantee to be the longest, although
	 * there might be situations in which these are not the most optimal paths
	 * to explore.
	 *
	 * Once a path been transformed it's removed from the query-graph and the process
	 * repeat itself. */

	/* A graph with no edges implies an empty algebraic expression
	 * the reasoning behind this decission is that algebraic expression
	 * represent graph traversals, no edges means no traversals. */
	AlgebraicExpression **exps = array_new(AlgebraicExpression *, 1);
	uint edge_count = QueryGraph_EdgeCount(qg);
	if(edge_count == 0) return exps;

	bool acyclic = IsAcyclicGraph(qg);
	QueryGraph *g = QueryGraph_Clone(qg);
	if(acyclic)
	{
		pathrecorder_init_mql();
		customized_filter_on_cycle_mql(g);
		int qgnode_len=array_len(qg->nodes);
		for(int i=0;i<qgnode_len;++i)
		{
			if(g->nodes[i]->customized_filter!=NULL)
				qg->nodes[i]->customized_filter=g->nodes[i]->customized_filter;
		}
		// customized_filter_on_cycle_mql(qg);
	}
	// As long as the query-graph isn't empty.
	while(QueryGraph_EdgeCount(g) > 0) {
		// Get leaf nodes at the deepest level.
		int depth;
		QGNode *n;
		if(acyclic) n = LongestPathTree(g, &depth); // Graph is a tree.
		else n = LongestPathGraph(g, &depth);       // Graph contains cycles.

		// Get a path of length level, allow closing a cycle if the graph is not acyclic.
		QGEdge **path = DFS(n, depth, !acyclic);
		uint path_len = array_len(path);
		assert(path_len == depth);

		/* TODO:
		 * In case path is a cycle, e.g. (b)-[]->(a)-[]->(b)
		 * make sure the first node on the path is referenced, _should_divide_expression(path, 0) is true.
		 * if this is not the case we will unnecessarily break the generated expression into 2 sub expressions
		 * while what we can do is simply rotate the cycle, (a)-[]->(b)-[]->(a)
		 * this is exactly the same only now we won't sub divide.
		 * Checking if path is a cycle done by testing the start and end node. */

		// Split path into sub paths.
		bool transpositions[path_len];
		_normalizePath(path, path_len, transpositions);

		customized_filter_mql(path,transpositions,qg);

		QGEdge ***paths = _Intermediate_Paths(path, qg);
		AlgebraicExpression **sub_exps = array_new(AlgebraicExpression *, 1);

		uint path_count = array_len(paths);
		uint edge_converted = 0;
		for(uint i = 0; i < path_count; i++) {
			// Construct expression.
			AlgebraicExpression *exp = _AlgebraicExpression_FromPath(paths[i], transpositions + edge_converted);
			edge_converted += array_len(paths[i]);
			sub_exps = array_append(sub_exps, exp);

			// if(i ==0 )
			// {
			// 	QGEdge **path=paths[0];
			// 	QGEdge *e =path[0];
			// 	if(e->src->customized_filter!=GrB_NULL){
			// 		exp = _AlgebraicExpression_MultiplyToTheLeft(AlgebraicExpression_NewOperand(e->src->customized_filter, true, e->src->alias, e->src->alias, NULL, e->src->label),exp);
			// 	}
			// }
			/* Remove exp[i] src label matrix (left most operand) as it's
			 * being used by exp[i-1] dest label matrix.
			 * (:A)-[:X]->(:B)-[:Y]->(:C)
			 * exp0: A * X * B
			 * exp1: B * Y * C
			 * should become
			 * exp0: A * X * B
			 * exp1: Y * C */
			if(i > 0) {
				// Make sure expression i follows previous expression.
				QGNode *src = QueryGraph_GetNodeByAlias(qg, AlgebraicExpression_Source(exp));
				if(src->label) {
					/* exp[i] shares a label matrix with exp[i-1]
					 * remove redundancy. */
					AlgebraicExpression *redundent = AlgebraicExpression_RemoveSource(&exp);
					AlgebraicExpression_Free(redundent);
				}
			}
			// Expression can not be empty.
			assert(AlgebraicExpression_OperandCount(exp) > 0);
		}

		sub_exps = _AlgebraicExpression_IsolateVariableLenExps(qg, sub_exps);

		uint sub_count = array_len(sub_exps);
		for(uint i = 0; i < sub_count; i++) {
			AlgebraicExpression *exp = sub_exps[i];
			// Add constructed expression to return value.
			exps = array_append(exps, exp);
		}

		// Remove path from graph.
		_RemovePathFromGraph(g, path);

		// Clean up
		for(uint i = 0; i < path_count; i++) array_free(paths[i]);
		array_free(path);
		array_free(paths);
		array_free(sub_exps);

		/* If original graph contained a cycle
		 * see now after we've removed a path if this is still the case. */
		if(!acyclic) acyclic = IsAcyclicGraph(g);
	}

	QueryGraph_Free(g);
	return exps;
}
