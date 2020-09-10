/*
* Copyright 2018-2020 Redis Labs Ltd. and Contributors
*
* This file is available under the Redis Labs Source Available License Agreement
*/

#include "gtest.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "../../src/util/rmalloc.h"
#include "../../src/query_ctx.h"
#include "../../src/graph/query_graph.h"
#include "../../src/filter_tree/filter_tree.h"
#include "../../src/ast/ast_build_filter_tree.h"
#include "../../src/arithmetic/algebraic_expression.h"
#include "../../src/execution_plan/optimizations/traverse_order.h"

#ifdef __cplusplus
}
#endif

class TraversalOrderingTest: public ::testing::Test {
  protected:
	static void SetUpTestCase() {
		// Use the malloc family for allocations
		Alloc_Reset();

		// Prepare thread-local variables
		ASSERT_TRUE(QueryCtx_Init());
	}

	static void TearDownTestCase() {
	}

	AST *_build_ast(const char *query) {
		cypher_parse_result_t *parse_result = cypher_parse(query, NULL, NULL, CYPHER_PARSE_ONLY_STATEMENTS);
		AST *ast = AST_Build(parse_result);
		return ast;
	}

	FT_FilterNode *build_filter_tree_from_query(const char *query) {
		AST *ast = _build_ast(query);
		return AST_BuildFilterTree(ast);
	}
};

TEST_F(TraversalOrderingTest, FilterFirst) {
	/* Given the ordered (left to right) set of algebraic expression:
	 * { [AB], [BC], [CD] }
	 * Which represents the traversal:
	 * (A)->(B)->(C)->(D)
	 * And a set of filters:
	 * C.V = X.
	 *
	 * We can reorder this set in such away
	 * that filters are applied as early as possible.
	 *
	 * Here are all of the possible permutations of the set
	 * in which the filter is applied at the earliest step:
	 * { [CB], [BA], [CD] } (D)<-(C)->(B)->(A) (2 transposes)
	 * { [CB], [CD], [BA] } (D)<-(C)->(B)->(A) (2 transposes)
	 * { [CD], [CB], [BA] } (A)<-(B)<-(C)->(D) (2 transposes) */

	FT_FilterNode *filters;
	QGNode *A = QGNode_New("A");
	QGNode *B = QGNode_New("B");
	QGNode *C = QGNode_New("C");
	QGNode *D = QGNode_New("D");

	QGEdge *AB = QGEdge_New(A, B, "E", "AB");
	QGEdge *BC = QGEdge_New(B, C, "E", "BC");
	QGEdge *CD = QGEdge_New(C, D, "E", "CD");

	QueryGraph *qg = QueryGraph_New(4, 3);

	QueryGraph_AddNode(qg, A);
	QueryGraph_AddNode(qg, B);
	QueryGraph_AddNode(qg, C);
	QueryGraph_AddNode(qg, D);
	QueryGraph_ConnectNodes(qg, A, B, AB);
	QueryGraph_ConnectNodes(qg, B, C, BC);
	QueryGraph_ConnectNodes(qg, C, D, CD);

	AlgebraicExpression *set[3];
	AlgebraicExpression *ExpAB = AlgebraicExpression_NewOperand(GrB_NULL, false, "A", "B", NULL, NULL);
	AlgebraicExpression *ExpBC = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "C", NULL, NULL);
	AlgebraicExpression *ExpCD = AlgebraicExpression_NewOperand(GrB_NULL, false, "C", "D", NULL, NULL);

	// { [AB], [BC], [CD] }
	set[0] = ExpAB;
	set[1] = ExpBC;
	set[2] = ExpCD;

	filters = build_filter_tree_from_query(
				  "MATCH (A)-[]->(B)-[]->(C)-[]->(D) WHERE A.val = 1 RETURN *");

	orderExpressions(qg, set, 3, filters, NULL);
	ASSERT_EQ(set[0], ExpAB);
	ASSERT_EQ(set[1], ExpBC);
	ASSERT_EQ(set[2], ExpCD);

	FilterTree_Free(filters);

	// { [AB], [BC], [CD] }
	set[0] = ExpAB;
	set[1] = ExpBC;
	set[2] = ExpCD;

	filters = build_filter_tree_from_query("MATCH (A)-[]->(B)-[]->(C)-[]->(D) WHERE B.val = 1 RETURN *");

	orderExpressions(qg, set, 3, filters, NULL);

	ASSERT_STREQ(AlgebraicExpression_Source(set[0]), "B");

	FilterTree_Free(filters);

	// { [AB], [BC], [CD] }
	set[0] = ExpAB;
	set[1] = ExpBC;
	set[2] = ExpCD;

	filters = build_filter_tree_from_query("MATCH (A)-[]->(B)-[]->(C)-[]->(D) WHERE C.val = 1 RETURN *");

	orderExpressions(qg, set, 3, filters, NULL);
	ASSERT_STREQ(AlgebraicExpression_Source(set[0]), "C");

	FilterTree_Free(filters);

	// { [AB], [BC], [CD] }
	set[0] = ExpAB;
	set[1] = ExpBC;
	set[2] = ExpCD;

	filters = build_filter_tree_from_query("MATCH (A)-[]->(B)-[]->(C)-[]->(D) WHERE D.val = 1 RETURN *");

	orderExpressions(qg, set, 3, filters, NULL);

	ASSERT_STREQ(AlgebraicExpression_Source(set[0]), "D");

	// Clean up.
	FilterTree_Free(filters);
	AlgebraicExpression_Free(ExpAB);
	AlgebraicExpression_Free(ExpBC);
	AlgebraicExpression_Free(ExpCD);
	QueryGraph_Free(qg);
}

TEST_F(TraversalOrderingTest, SingleOptimalArrangement) {
	/* Given the set of algebraic expressions:
	 * { [AB], [BC], [BD] }
	 * We represent the traversal:
	 * (A:L {v: 1})->(B)->(C), (B)->(D:L {v: 1})
	 * The optimal order of traversals should always be:
	 * { [AB], [BD], [BC] }
	 * Validate this for all input permutations.
	 */

	FT_FilterNode *filters;
	QGNode *A = QGNode_New("A");
	QGNode *B = QGNode_New("B");
	QGNode *C = QGNode_New("C");
	QGNode *D = QGNode_New("D");
	A->label = "L";
	D->label = "L";

	QGEdge *AB = QGEdge_New(A, B, "E", "AB");
	QGEdge *BC = QGEdge_New(B, C, "E", "BC");
	QGEdge *BD = QGEdge_New(B, D, "E", "BD");

	QueryGraph *qg = QueryGraph_New(4, 3);

	QueryGraph_AddNode(qg, A);
	QueryGraph_AddNode(qg, B);
	QueryGraph_AddNode(qg, C);
	QueryGraph_AddNode(qg, D);
	QueryGraph_ConnectNodes(qg, A, B, AB);
	QueryGraph_ConnectNodes(qg, B, C, BC);
	QueryGraph_ConnectNodes(qg, B, D, BD);

	AlgebraicExpression *ExpAB = AlgebraicExpression_NewOperand(GrB_NULL, false, "A", "B", NULL, "L");
	AlgebraicExpression *ExpBC = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "C", NULL, NULL);
	AlgebraicExpression *ExpBD = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "D", NULL, "L");

	filters = build_filter_tree_from_query(
				  "MATCH (A:L {v: 1})-[]->(B)-[]->(C), (B)-[]->(D:L {v: 1})) RETURN *");

	AlgebraicExpression *set[3] = {ExpAB, ExpBC, ExpBD};
	std::sort(set, set + 3);
	// Test every permutation of the set.
	do {
		AlgebraicExpression *tmp[3] = {set[0], set[1], set[2]};
		orderExpressions(qg, tmp, 3, filters, NULL);
		ASSERT_EQ(tmp[0], ExpAB);
		ASSERT_EQ(tmp[1], ExpBD);
		ASSERT_EQ(tmp[2], ExpBC);
	} while(std::next_permutation(set, set + 3));

	// Clean up.
	FilterTree_Free(filters);
	AlgebraicExpression_Free(ExpAB);
	AlgebraicExpression_Free(ExpBC);
	AlgebraicExpression_Free(ExpBD);
	QueryGraph_Free(qg);
}

static void _print_nodes(QGNode **nodes, int node_count) { // TODO tmp
	for(int i = 0; i < node_count; i ++) {
		printf("%s:%s, ", nodes[i]->alias, nodes[i]->label ? : "");
	}
	// for(int i = 0; i < node_count; i ++) {
	// for(int j = 0; j < node_count; j ++) {
	// if(nodes[j]->alias[0] == i + '0') printf("%s:%s, ", nodes[j]->alias, nodes[j]->label ? : "");
	// }
	// }
	printf("\n");
}

static void truth_table(bool *table, int n) {
	int idx = 0;
	for(int i = 0; i < 1 << n; i ++) {
		for(int j = 0; j < sizeof(int); j ++) {
			bool x = (i & (1 << j));
			table[idx++] = x;
		}
	}
}

TEST_F(TraversalOrderingTest, ValidateLabelScoring) {
	/* Given the graph with the structure:
	 * (A)->(B)->(C)->(D), (D)->(A), (B)->(D)
	 * Test all permutations of labeled nodes and
	 * validate that all produce an optimal scoring. */
	uint exp_count = 5;
	int node_count = 4;
	int edge_count = 5;
	QueryGraph *qg = QueryGraph_New(node_count, edge_count);

	// Build QGNodes
	QGNode *nodes[node_count];
	for(int i = 0; i < node_count + 1; i ++) {
		char *alias;
		asprintf(&alias, "%c", i + 'A');
		nodes[i] = QGNode_New(alias);
		QueryGraph_AddNode(qg, nodes[i]);
	}

	QGEdge *AB = QGEdge_New(NULL, NULL, "E", "AB");
	QueryGraph_ConnectNodes(qg, nodes[0], nodes[1], AB);
	QGEdge *BC = QGEdge_New(NULL, NULL, "E", "BC");
	QueryGraph_ConnectNodes(qg, nodes[1], nodes[2], BC);
	QGEdge *CD = QGEdge_New(NULL, NULL, "E", "CD");
	QueryGraph_ConnectNodes(qg, nodes[2], nodes[3], CD);
	QGEdge *CA = QGEdge_New(NULL, NULL, "E", "CA");
	QueryGraph_ConnectNodes(qg, nodes[3], nodes[0], CA);
	QGEdge *BD = QGEdge_New(NULL, NULL, "E", "BD");
	QueryGraph_ConnectNodes(qg, nodes[1], nodes[3], BD);

	AlgebraicExpression *ExpAB = AlgebraicExpression_NewOperand(GrB_NULL, false, "A", "B", NULL, NULL);
	AlgebraicExpression *ExpBC = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "C", NULL, NULL);
	AlgebraicExpression *ExpCD = AlgebraicExpression_NewOperand(GrB_NULL, false, "C", "D", NULL, NULL);
	AlgebraicExpression *ExpCA = AlgebraicExpression_NewOperand(GrB_NULL, false, "C", "A", NULL, NULL);
	AlgebraicExpression *ExpBD = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "D", NULL, NULL);

	// Store unmodified versions of operands, as ordering can modify them by introducing transposes.
	AlgebraicExpression *ExpAB_orig = AlgebraicExpression_NewOperand(GrB_NULL, false, "A", "B", NULL,
																	 NULL);
	AlgebraicExpression *ExpBC_orig = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "C", NULL,
																	 NULL);
	AlgebraicExpression *ExpCD_orig = AlgebraicExpression_NewOperand(GrB_NULL, false, "C", "D", NULL,
																	 NULL);
	AlgebraicExpression *ExpCA_orig = AlgebraicExpression_NewOperand(GrB_NULL, false, "C", "A", NULL,
																	 NULL);
	AlgebraicExpression *ExpBD_orig = AlgebraicExpression_NewOperand(GrB_NULL, false, "B", "D", NULL,
																	 NULL);
	AlgebraicExpression *orig_set[5] = {ExpAB_orig, ExpBC_orig, ExpCD_orig, ExpCA_orig, ExpBD_orig};

	int n = (2 << node_count) * 4; // 16 (rows) * 4 (cols)
	bool to_label[n];
	truth_table(to_label, node_count);
	AlgebraicExpression *set[5];
	memcpy(set, orig_set, exp_count * sizeof(AlgebraicExpression *));
	for(int i = 0; i < 1 << node_count; i ++) {
		// Label the appropriate nodes in the sequence.
		for(int j = 0; j < node_count; j ++) {
			if(to_label[i * 4 + j]) nodes[j]->label = "L";
			else nodes[j]->label = NULL;
		}

		// Order the expressions and obtain the score.
		orderExpressions(qg, set, exp_count, NULL, NULL);
		int first_score = score_arrangement(set, exp_count, qg, NULL, NULL);
		memcpy(set, orig_set, exp_count * sizeof(AlgebraicExpression *));

		// Test every permutation of the set.
		std::sort(set, set + exp_count);
		do {
			AlgebraicExpression *tmp[exp_count];
			memcpy(tmp, orig_set, exp_count * sizeof(AlgebraicExpression *));
			orderExpressions(qg, tmp, exp_count, NULL, NULL);
			int score = score_arrangement(tmp, exp_count, qg, NULL, NULL);
			if(score != first_score) { // TODO tmp
				printf("Scored %d, first score was %d\n", score, first_score);
			}
			ASSERT_EQ(score, first_score);
		} while(std::next_permutation(set, set + exp_count));
	}
	// Clean up.
	AlgebraicExpression_Free(ExpAB);
	AlgebraicExpression_Free(ExpBC);
	AlgebraicExpression_Free(ExpCD);
	AlgebraicExpression_Free(ExpCA);
	AlgebraicExpression_Free(ExpBD);
	QueryGraph_Free(qg);
}

