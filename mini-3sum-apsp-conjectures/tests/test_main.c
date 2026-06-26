#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "threesum.h"
#include "apsp.h"
#include "minplus.h"
#include "reduction.h"
#include "complexity.h"
#include "equivalence.h"
static int tr=0,tp=0;
#define T(n) do{tr++;printf("  TEST %s... ",n);}while(0)
#define P() do{printf("PASS\n");tp++;}while(0)
#define F(m) do{printf("FAIL: %s\n",m);}while(0)
#define C(c,m) do{if(c)P();else F(m);}while(0)
static void t1(void){T("3SUM naive");ts_elem_t v[]={-1,0,1,2,-2,3};ts_instance_t*i=ts_instance_create(v,6);ts_result_t r=ts_naive_O_n3(i);C(r.found&&r.count>=2,"triples");ts_result_free(&r);ts_instance_destroy(i);}
static void t2(void){T("3SUM hash");ts_elem_t v[]={-3,0,2,4,1,-4};ts_instance_t*i=ts_instance_create(v,6);ts_result_t a=ts_naive_O_n3(i),b=ts_quadratic_hash(i);C(a.count==b.count,"match");ts_result_free(&a);ts_result_free(&b);ts_instance_destroy(i);}
static void t3(void){T("3SUM two-ptr");ts_elem_t v[]={-4,-2,0,1,2,3,5};ts_instance_t*i=ts_instance_create(v,7);ts_result_t a=ts_naive_O_n3(i),b=ts_quadratic_two_pointer(i);C(a.count==b.count,"match");ts_result_free(&a);ts_result_free(&b);ts_instance_destroy(i);}
static void t4(void){T("3SUM sort+bs");ts_elem_t v[]={-2,-1,0,1,2,3};ts_instance_t*i=ts_instance_create(v,6);ts_result_t a=ts_naive_O_n3(i);ts_result_t b=ts_quadratic_sort_bsearch(i);C(a.found==b.found,"ok");ts_result_free(&a);ts_result_free(&b);ts_instance_destroy(i);}
static void t5(void){T("3SUM count");ts_elem_t v[]={-2,-1,0,1,2,3};ts_instance_t*i=ts_instance_create(v,6);int64_t c=ts_count_zero_sum_triples(i);ts_result_t r=ts_naive_O_n3(ts_instance_create(v,6));C(c==r.count,"match");ts_result_free(&r);ts_instance_destroy(i);}
static void t6(void){T("collinearity");double x[]={0,1,2,3,5},y[]={0,1,2,3,10};C(ts_collinearity_test(x,y,5),"collinear");}
static void t7(void){T("3SUM conjecture");ts_conjecture_t c=ts_conjecture_status();C(c.statement!=NULL,"defined");}
static void t8(void){T("APSP Floyd");apsp_adjacency_t*a=apsp_adjacency_create(4);apsp_adjacency_set_edge(a,0,1,3);apsp_adjacency_set_edge(a,0,2,8);apsp_adjacency_set_edge(a,1,2,4);apsp_adjacency_set_edge(a,1,3,2);apsp_adjacency_set_edge(a,2,3,1);apsp_graph_t*g=apsp_floyd_warshall(a);C(g&&!g->has_negative_cycle&&fabs(g->dist[3]-5)<1e-6,"dist(0,3)=5");apsp_graph_destroy(g);apsp_adjacency_destroy(a);}
static void t9(void){T("APSP Dijkstra");apsp_adjacency_t*a=apsp_adjacency_create(4);apsp_adjacency_set_edge(a,0,1,2);apsp_adjacency_set_edge(a,1,2,3);apsp_adjacency_set_edge(a,2,3,1);apsp_adjacency_set_edge(a,0,3,10);apsp_graph_t*g=apsp_dijkstra_all(a);C(fabs(g->dist[3]-6)<1e-6,"dist(0,3)=6");apsp_graph_destroy(g);apsp_adjacency_destroy(a);}
static void t10(void){T("Min-plus mult");mp_matrix_t*A=mp_matrix_create(3,3);mp_set(A,0,0,0);mp_set(A,0,1,2);mp_set(A,1,1,0);mp_set(A,1,2,1);mp_set(A,2,2,0);mp_matrix_t*C=mp_naive_multiply(A,A);C(C!=NULL,"ok");mp_matrix_destroy(A);mp_matrix_destroy(C);}
static void t11(void){T("Min-plus closure");mp_matrix_t*A=mp_matrix_create(4,4);mp_set(A,0,0,0);mp_set(A,0,1,3);mp_set(A,1,1,0);mp_set(A,1,2,1);mp_set(A,2,2,0);mp_set(A,2,3,2);mp_set(A,3,3,0);mp_matrix_t*D=mp_closure(A);C(mp_get(D,0,3)<=6,"ok");mp_matrix_destroy(A);mp_matrix_destroy(D);}
static void t12(void){T("Reduction info");reduction_info_t i=reduction_get_info(RED_3SUM_TO_APSP);C(i.name!=NULL,"exists");}
static void t13(void){T("Equivalence verify");C(reduction_verify_subcubic_equivalence(),"ok");}
static void t14(void){T("Complexity scaling");cx_measurement_t p[3]={{100,0.001,1000000,0},{200,0.004,8000000,0},{400,0.016,64000000,0}};cx_scaling_result_t*r=cx_estimate_exponent(p,3);C(r&&fabs(r->estimated_exponent-2)<0.5,"~2");cx_scaling_result_destroy(r);}
static void t15(void){T("Subcubic check");C(!cx_is_subcubic(3.0,100)&&cx_is_subcubic(2.9,1000),"checks");}
static void t16(void){T("Problem names");C(strcmp(eq_problem_name(EQ_APSP),"APSP")==0&&strcmp(eq_problem_name(EQ_3SUM),"3SUM")==0,"names");}
static void t17(void){T("Equiv classes");C(eq_are_equivalent(EQ_APSP,EQ_MIN_PLUS_PRODUCT)&&eq_are_equivalent(EQ_3SUM,EQ_COLLINEAR),"classes");}
static void t18(void){T("Hexagon");C(eq_verify_hexagon(),"ok");}
int main(void){printf("=== mini-3sum-apsp-conjectures Test Suite ===\n\n");t1();t2();t3();t4();t5();t6();t7();t8();t9();t10();t11();t12();t13();t14();t15();t16();t17();t18();printf("\n=== Results: %d/%d tests passed ===\n",tp,tr);return tp>=tr?0:1;}
