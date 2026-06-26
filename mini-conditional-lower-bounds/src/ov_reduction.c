/* ov_reduction.c -- Orthogonal Vectors Implementation */
#include "ov_reduction.h"
#include "condlb.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

VectorSet* vs_create(int capacity, int dim) {
    assert(capacity > 0 && dim > 0);
    VectorSet* vs = (VectorSet*)malloc(sizeof(VectorSet));
    if (!vs) return NULL;
    vs->vectors = (BinaryVector*)calloc((size_t)capacity, sizeof(BinaryVector));
    if (!vs->vectors) { free(vs); return NULL; }
    vs->n_vectors = 0; vs->dim = dim; vs->capacity = capacity;
    return vs;
}

void vs_free(VectorSet* set) {
    if (!set) return;
    if (set->vectors) {
        for (int i = 0; i < set->n_vectors; i++) free(set->vectors[i].bits);
        free(set->vectors);
    }
    free(set);
}

int vs_add_vector(VectorSet* set, const int* bits) {
    assert(set && bits);
    if (set->n_vectors >= set->capacity) {
        int nc = set->capacity * 2;
        BinaryVector* nv = (BinaryVector*)realloc(set->vectors, (size_t)nc * sizeof(BinaryVector));
        if (!nv) return -1;
        memset(nv + set->capacity, 0, (size_t)(nc - set->capacity) * sizeof(BinaryVector));
        set->vectors = nv; set->capacity = nc;
    }
    int idx = set->n_vectors++;
    set->vectors[idx].dim = set->dim; set->vectors[idx].id = idx;
    set->vectors[idx].bits = (int*)malloc((size_t)set->dim * sizeof(int));
    if (!set->vectors[idx].bits) return -1;
    memcpy(set->vectors[idx].bits, bits, (size_t)set->dim * sizeof(int));
    return idx;
}

uint64_t vs_dot_product(const VectorSet* set, int i, int j) {
    assert(set); assert(i >= 0 && i < set->n_vectors); assert(j >= 0 && j < set->n_vectors);
    uint64_t dot = 0;
    int* a = set->vectors[i].bits; int* b = set->vectors[j].bits;
    for (int k = 0; k < set->dim; k++) dot += (uint64_t)(a[k] & b[k]);
    return dot;
}

int ov_brute_force(const VectorSet* set, OrthogonalityResult* result) {
    assert(set); int n = set->n_vectors;
    if (result) result->found = 0;
    for (int i = 0; i < n; i++)
        for (int j = i + 1; j < n; j++)
            if (vs_dot_product(set, i, j) == 0) {
                if (result) { result->i=i; result->j=j; result->dot_product=0; result->orthogonal=1; result->found=1; }
                return 1;
            }
    return 0;
}

int ov_split_and_list(const VectorSet* set, double c_log_factor, OrthogonalityResult* result) {
    assert(set && c_log_factor > 0);
    int n = set->n_vectors, d = set->dim;
    int k = (int)(c_log_factor * log((double)n) / log(2.0));
    if (k < 2) k = 2; if (k > d) k = d;
    int gs = d / k, rem = d % k;
    int** po = (int**)malloc((size_t)n * sizeof(int*));
    if (!po) return 0;
    for (int i = 0; i < n; i++) {
        po[i] = (int*)calloc((size_t)k, sizeof(int));
        if (!po[i]) { for (int p=0;p<i;p++) free(po[p]); free(po); return 0; }
        int* bits = set->vectors[i].bits; int off = 0;
        for (int g = 0; g < k; g++) {
            int sz = gs + (g < rem ? 1 : 0); int sum = 0;
            for (int t = 0; t < sz; t++) sum += bits[off + t];
            po[i][g] = sum; off += sz;
        }
    }
    int found = 0;
    for (int i = 0; i < n && !found; i++)
        for (int j = i + 1; j < n && !found; j++) {
            int poss = 1;
            for (int g = 0; g < k; g++) {
                int sz = gs + (g < rem ? 1 : 0);
                if (po[i][g] == sz && po[j][g] == sz) { poss = 0; break; }
            }
            if (poss && vs_dot_product(set,i,j) == 0) {
                if (result) { result->i=i; result->j=j; result->dot_product=0; result->orthogonal=1; result->found=1; }
                found = 1;
            }
        }
    for (int i = 0; i < n; i++) free(po[i]);
    free(po); return found;
}

int ov_light_bulb(const VectorSet* set, OrthogonalityResult* result) {
    assert(set); int n = set->n_vectors, d = set->dim;
    int* sp = (int*)malloc((size_t)n * sizeof(int));
    if (!sp) return 0;
    for (int i = 0; i < n; i++) { int s=0; for (int k=0;k<d;k++) s+=set->vectors[i].bits[k]; sp[i]=s; }
    int found = 0;
    for (int i = 0; i < n && !found; i++) {
        if (sp[i] > d/2) continue;
        for (int j = i+1; j < n && !found; j++) {
            if (sp[j] > d/2) continue;
            if (vs_dot_product(set,i,j) == 0) {
                if (result) { result->i=i;result->j=j;result->dot_product=0;result->orthogonal=1;result->found=1; }
                found=1;
            }
        }
    }
    free(sp); return found;
}

int ov_fft_method(const VectorSet* set, OrthogonalityResult* result) {
    assert(set); int n=set->n_vectors, d=set->dim;
    if (d>20) return ov_brute_force(set,result);
    int h=d/2, h2=d-h, s1=1<<h, s2=1<<h2;
    int* c1=(int*)calloc((size_t)s1,sizeof(int));
    int* f1=(int*)malloc((size_t)s1*sizeof(int));
    if(!c1||!f1){free(c1);free(f1);return 0;}
    for(int i=0;i<s1;i++)f1[i]=-1;
    for(int i=0;i<n;i++){int m=0;for(int k=0;k<h;k++)if(set->vectors[i].bits[k])m|=(1<<k);c1[m]++;if(f1[m]==-1)f1[m]=i;}
    int* c2=(int*)calloc((size_t)s2,sizeof(int));
    int* f2=(int*)malloc((size_t)s2*sizeof(int));
    if(!c2||!f2){free(c1);free(f1);free(c2);free(f2);return 0;}
    for(int i=0;i<s2;i++)f2[i]=-1;
    for(int i=0;i<n;i++){int m=0;for(int k=0;k<h2;k++)if(set->vectors[i].bits[h+k])m|=(1<<k);c2[m]++;if(f2[m]==-1)f2[m]=i;}
    int found=0;
    for(int m1=0;m1<s1&&!found;m1++){if(!c1[m1])continue;
        for(int m2=0;m2<s2&&!found;m2++){if(!c2[m2])continue;
            if(vs_dot_product(set,f1[m1],f2[m2])==0){if(result){result->i=f1[m1];result->j=f2[m2];result->dot_product=0;result->orthogonal=1;result->found=1;}found=1;}
        }}
    free(c1);free(f1);free(c2);free(f2);return found;
}

VectorSet* vs_random(int n, int dim, double density) {
    assert(n>0&&dim>0&&density>=0.0&&density<=1.0);
    VectorSet* vs=vs_create(n,dim); if(!vs)return NULL;
    srand(42);
    for(int i=0;i<n;i++){int* bits=(int*)malloc((size_t)dim*sizeof(int));
        if(!bits){vs_free(vs);return NULL;}
        for(int k=0;k<dim;k++)bits[k]=((double)rand()/RAND_MAX<density)?1:0;
        vs_add_vector(vs,bits);free(bits);}
    return vs;
}

int vs_has_orthogonal_pair(const VectorSet* set){return ov_brute_force(set,NULL);}

int vs_find_all_pairs(const VectorSet* set, OrthogonalityResult* results, int max_results) {
    assert(set&&results); int n=set->n_vectors,count=0;
    for(int i=0;i<n&&count<max_results;i++)
        for(int j=i+1;j<n&&count<max_results;j++)
            if(vs_dot_product(set,i,j)==0){results[count].i=i;results[count].j=j;results[count].dot_product=0;results[count].orthogonal=1;results[count].found=1;count++;}
    return count;
}

void ov_to_edit_distance(const VectorSet* set, char** out_s1, int* len1, char** out_s2, int* len2) {
    assert(set); int n=set->n_vectors, d=set->dim;
    int gl=d+2; *len1=n*gl; *len2=n*gl;
    *out_s1=(char*)malloc((size_t)(*len1+1)); *out_s2=(char*)malloc((size_t)(*len2+1));
    if(!*out_s1||!*out_s2){free(*out_s1);free(*out_s2);return;}
    for(int i=0;i<n;i++){int base=i*gl;
        (*out_s1)[base]='V';(*out_s2)[base]='V';
        for(int k=0;k<d;k++){(*out_s1)[base+1+k]=set->vectors[i].bits[k]?'1':'0';(*out_s2)[base+1+k]=set->vectors[i].bits[k]?'0':'1';}
        (*out_s1)[base+1+d]='#';(*out_s2)[base+1+d]='#';}
    (*out_s1)[*len1]=0;(*out_s2)[*len2]=0;
}

void ov_to_lcs(const VectorSet* set, char** out_s1, int* len1, char** out_s2, int* len2) {
    assert(set); int n=set->n_vectors, d=set->dim;
    int bl=3*d; *len1=n*bl; *len2=n*bl;
    *out_s1=(char*)malloc((size_t)(*len1+1)); *out_s2=(char*)malloc((size_t)(*len2+1));
    if(!*out_s1||!*out_s2){free(*out_s1);free(*out_s2);return;}
    for(int i=0;i<n;i++){int base=i*bl;
        for(int k=0;k<d;k++){char c=set->vectors[i].bits[k]?'A':'B';
            (*out_s1)[base+3*k]=c;(*out_s1)[base+3*k+1]='X';(*out_s1)[base+3*k+2]='X';
            (*out_s2)[base+3*k]=set->vectors[i].bits[k]?'B':'A';(*out_s2)[base+3*k+1]='X';(*out_s2)[base+3*k+2]='X';}}
    (*out_s1)[*len1]=0;(*out_s2)[*len2]=0;
}

void ov_to_diameter(const VectorSet* set, int** out_adj_matrix, int* n_vertices) {
    assert(set); int n=set->n_vectors; *n_vertices=2*n+2; int N=*n_vertices;
    int* adj=(int*)calloc((size_t)(N*N),sizeof(int)); if(!adj)return;
    for(int i=0;i<n;i++){adj[0*N+(1+i)]=1;adj[(1+n+i)*N+(N-1)]=1;}
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)if(vs_dot_product(set,i,j)>0)adj[(1+i)*N+(1+n+j)]=1;
    *out_adj_matrix=adj;
}

void ov_to_subgraph_isomorphism(const VectorSet* set, int** out_host_adj, int* host_n, int** out_pat_adj, int* pat_n) {
    assert(set); int n=set->n_vectors, d=set->dim;
    *host_n=n+d+2; int HN=*host_n;
    int* host=(int*)calloc((size_t)(HN*HN),sizeof(int)); if(!host)return;
    *pat_n=4; int PN=*pat_n;
    int* pat=(int*)calloc((size_t)(PN*PN),sizeof(int)); if(!pat){free(host);return;}
    for(int i=0;i<n;i++)for(int k=0;k<d;k++)if(set->vectors[i].bits[k]){host[i*HN+(n+k)]=1;host[(n+k)*HN+i]=1;}
    pat[0*PN+1]=1;pat[1*PN+0]=1;pat[1*PN+2]=1;pat[2*PN+1]=1;pat[2*PN+3]=1;pat[3*PN+2]=1;pat[3*PN+0]=1;pat[0*PN+3]=1;
    *out_host_adj=host;*out_pat_adj=pat;
}

int ov_dynamic_string_matching_bound(const VectorSet* set, int* op_count, int* query_count) {
    assert(set); int n=set->n_vectors;
    if(op_count)*op_count=n; if(query_count)*query_count=n;
    return (set->dim>2.0*log((double)n))?1:0;
}

int ov_reachability_oracle_bound(const VectorSet* set, double* space_bound, double* query_time_bound) {
    assert(set); int n=set->n_vectors; double ln=log((double)n);
    if(space_bound)*space_bound=(double)n*ln;
    if(query_time_bound)*query_time_bound=ln;
    return (set->dim>2*(int)ln)?1:0;
}

int ov_regex_matching_bound(const VectorSet* set, double* lower_bound_exponent) {
    assert(set);
    if(lower_bound_exponent)*lower_bound_exponent=2.0;
    return 1;
}
