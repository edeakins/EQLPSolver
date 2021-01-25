/*
 * saucy.c
 * Searching for Automorphisms in Underlying CNF, yes?
 *
 * by Paul T. Darga <pdarga@umich.edu>
 * and Mark Liffiton <liffiton@umich.edu>
 * and Hadi Katebi <hadik@eecs.umich.edu>
 * 
 * Copyright (C) 2004, The Regents of the University of Michigan
 * See the LICENSE file for details.
 */


#include <stdlib.h> /* malloc, calloc, and free */
//#include <string.h> /* memcpy */
//#include <stdint.h>
#include <cstdlib> /* malloc, calloc, and free */
#include <cstring> /* memcpy */
#include <cstdint>
#include <iostream>

#include "saucy-equitable.h"

//using namespace std;

/* #define DCCOUNT(i, j) (s->dccount[(i)*s->e + (j)]) */
#define DCCOUNT(i, j) (s->dccount[(i)*s->wcount + (j)])
/* TODO: consider using unsigned ints/longs for everything */

struct coloring {
    int *lab;        /* Labelling of objects */
    int *unlab;      /* Inverse of lab */
    int *cfront;     /* Pointer to front of cells */
    int *clen;       /* Length of cells (defined for cfront's) */
    int *parent;     /* Contains parent cells from splitting */
};

struct saucy {
    /* Graph data */
    uint64_t n;           /* Size of domain */
    /*int e;*/           /* number of edges */
    uint64_t wcount;      /* number of edge colors */
    const int *adj;  /* Neighbors of k: edg[adj[k]]..edg[adj[k+1]] */
    const int *edg;  /* Actual neighbor data */
    const int *wght;  /* Actual edge colors */
    const int *dadj; /* Fanin neighbor indices, for digraphs */
    const int *dedg; /* Fanin neighbor data, for digraphs */
    const int *dwght; /* Fanin weight data, for digraphs */
    void *arg;       /* Opaque client data */

    /* Coloring data */
    struct coloring left, right;
    int *nextnon;    /* Forward next-nonsingleton pointers */ 
    int *prevnon;    /* Backward next-nonsingleton pointers */

    /* Refinement: inducers */
    char *indmark;   /* Induce marks */
    int *ninduce;    /* Nonsingletons that might induce refinement */
    int *sinduce;    /* Singletons that might induce refinement */
    int nninduce;    /* Size of ninduce stack */
    int nsinduce;    /* Size of sinduce stack */

    /* Refinement: marked cells */
    int *clist;      /* List of cells marked for refining */
    int csize;       /* Number of cells in clist */

    /* Refinement: workspace */
    char *stuff;     /* Bit vector, but one char per bit */
    int *wstuff;     /* weight data */
    int *ccount;     /* Number of connections to refining cell */
    int *dccount;     /* Number of connections to refining cell */
    int *bucket;     /* Workspace */
    int *count;      /* Num vertices with same adj count to ref cell */
    int *junk;       /* More workspace */
    int *diffL;   /* Weight workspace */
    int *gamma;      /* Working permutation */
    int *conncnts;   /* Connection counts for cell fronts */

    /* Search data */
    int lev;         /* Current search tree level */
    int anc;         /* Level of greatest common ancestor with zeta */
    int *anctar;     /* Copy of target cell at anc */
    int kanctar;     /* Location within anctar to iterate from */
    int *start;      /* Location of target at each level */
    int indmin;      /* Used for group size computation */
    int match;       /* Have we not diverged from previous left? */

    /* Search: orbit partition */
    int *theta;      /* Running approximation of orbit partition */
    int *thsize;     /* Size of cells in theta, defined for mcrs */
    int *thnext;     /* Next rep in list (circular list) */
    int *thprev;     /* Previous rep in list */
    int *threp;      /* First rep for a given cell front */
    int *thfront;    /* The cell front associated with this rep */

    /* Search: split record */
    int *splitwho;   /* List of where splits occurred */
    int *splitfrom;  /* List of cells which were split */
    int *splitlev;   /* Where splitwho/from begins for each level */
    int nsplits;     /* Number of splits at this point */

    /* Search: differences from leftmost */
    char *diffmark;  /* Marked for diff labels */
    int *diffs;      /* List of diff labels */
    int *difflev;    /* How many labels diffed at each level */
    int ndiffs;      /* Current number of diffs */
    int *undifflev;  /* How many diff labels fixed at each level */
    int nundiffs;    /* Current number of diffs in singletons (fixed) */
    int *unsupp;     /* Inverted diff array */
    int *specmin;    /* Speculated mappings */
    int *pairs;      /* Not-undiffed diffs that can make two-cycles */
    int *unpairs;    /* Indices into pairs */
    int npairs;      /* Number of such pairs */
    int *diffnons;   /* Diffs that haven't been undiffed */
    int *undiffnons; /* Inverse of that */
    int ndiffnons;   /* Number of such diffs */

    /* Polymorphic functions */
    saucy_consumer *consumer;
    int (*split)(struct saucy *, struct coloring *, int, int);
    int (*is_automorphism)(struct saucy *);
    int (*ref_singleton)(struct saucy *, struct coloring *, int);
    int (*ref_nonsingle)(struct saucy *, struct coloring *, int);

     /* Statistics structure */
    struct saucy_stats *stats;
};

static int
array_find_min(const int *a, int n)
{
    const int *start = a, *end = a + n, *min = a;
    while (++a != end) {
        if (*a < *min) min = a;
    }
    return min - start;
}

/*
static int
list_search( int *list, int key, int len, int *ndx )
{
    for( *ndx = 0; *ndx < len; ++(*ndx)){
        if( list[*ndx] == key ) 
            return 1;
    }

    return 0;
}
*/

static void
swap(int *a, int x, int y)
{
    int tmp = a[x];
    a[x] = a[y];
    a[y] = tmp;
}

static void
sift_up(int *a, int k)
{
    int p;
    do {
        p = k / 2;
        if (a[k] <= a[p]) {
            return;
        }
        else {
            swap(a, k, p);
            k = p;
        }
    } while (k > 1);
}

static void
sift_down(int *a, int n)
{
    int p = 1, k = 2;
    while (k <= n) {
        if (k < n && a[k] < a[k+1]) ++k;
        if (a[p] < a[k]) {
            swap(a, p, k);
            p = k;
            k = 2 * p;
        }
        else {
            return;
        }
    }
}

static void
heap_sort(int *a, int n)
{
    int i;
    for (i = 1; i < n; ++i) {
        sift_up(a-1, i+1);
    }
    --i;
    while (i > 0) {
        swap(a, 0, i);
        sift_down(a-1, i--);
    }
}

static void
insertion_sort(int *a, int n)
{
    int i, j, k;
    for (i = 1; i < n; ++i) {
        k = a[i];
        for (j = i; j > 0 && a[j-1] > k; --j) {
            a[j] = a[j-1];
        }
        a[j] = k;
    }
}

static int
partition(int *a, int n, int m)
{
    int f = 0, b = n;
    for (;;) {
        while (a[f] <= m) ++f;
        do  --b; while (m <= a[b]);
        if (f < b) {
            swap(a, f, b);
            ++f;
        }
        else break;
    }
    return f;
}

static int
log_base2(int n)
{
    int k = 0;
    while (n > 1) {
        ++k;
        n >>= 1;
    }
    return k;
}

static int
median(int a, int b, int c)
{
    if (a <= b) {
        if (b <= c) return b;
        if (a <= c) return c;
        return a;
    }
    else {
        if (a <= c) return a;
        if (b <= c) return c;
        return b;
    }
}

static void
introsort_loop(int *a, int n, int lim)
{
    int p;
    while (n > 16) {
        if (lim == 0) {
            heap_sort(a, n);
            return;
        }
        --lim;
        p = partition(a, n, median(a[0], a[n/2], a[n-1]));
        introsort_loop(a + p, n - p, lim);
        n = p;
    }
}

static void
introsort(int *a, int n)
{
    introsort_loop(a, n, 2 * log_base2(n));
    insertion_sort(a, n);
}

static int
do_find_min(struct coloring *c, int t)
{
    return array_find_min(c->lab + t, c->clen[t] + 1) + t;
}

static int
find_min(struct saucy *s, int t)
{
    return do_find_min(&s->right, t);
}

static void
set_label(struct coloring *c, int index, int value)
{
    c->lab[index] = value;
    c->unlab[value] = index;
}

static void
swap_labels(struct coloring *c, int a, int b)
{
    int tmp = c->lab[a];
    set_label(c, a, c->lab[b]);
    set_label(c, b, tmp);
}

static void
move_to_back(struct saucy *s, struct coloring *c, int k)
{
    int cf = c->cfront[k];
    int cb = cf + c->clen[cf];
    int offset = s->conncnts[cf]++;

    /* Move this connected label to the back of its cell */
    swap_labels(c, cb - offset, c->unlab[k]);

    /* Add it to the cell list if it's the first one swapped */
    if (!offset) s->clist[s->csize++] = cf;
}

static void
data_mark(struct saucy *s, struct coloring *c, int k)
{
    int cf = c->cfront[k];

    /* Move connects to the back of nonsingletons */
    if (c->clen[cf]) move_to_back(s, c, k);
}

static void
data_count(struct saucy *s, struct coloring *c, int k)
{
    int cf = c->cfront[k];

    /* Move to back and count the number of connections */
    if (c->clen[cf] && !s->ccount[k]++) move_to_back(s, c, k);
}

static int
check_mapping(struct saucy *s, const int *adj, const int *edg, const int *wght, int k)
{
    int i, gk, ret = 1;

    /* Mark gamma of neighbors */
    for (i = adj[k]; i != adj[k+1]; ++i) {
        s->stuff[s->gamma[edg[i]]] = 1;
        /* wstuff should only need n spots since a vertex */
        /* can connect to at most n vertices */
        s->wstuff[s->gamma[edg[i]]] = wght[i];
    }

    /* Check neighbors of gamma */
    gk = s->gamma[k];
    for (i = adj[gk]; ret && i != adj[gk+1]; ++i) {
        ret = s->stuff[edg[i]];
        /* TODO: verify that this is a valid test for weight data */
        ret = ret && (wght[i] == s->wstuff[edg[i]]);
    }

    /* Clear out bit vector before we leave */
    for (i = adj[k]; i != adj[k+1]; ++i) {
        s->stuff[s->gamma[edg[i]]] = 0;
        s->wstuff[s->gamma[edg[i]]] = 0;
    }

    return ret;
}

static int
is_undirected_automorphism(struct saucy *s)
{
    int i, j;

    for (i = 0; i < s->ndiffs; ++i) {
        j = s->unsupp[i];
        if (!check_mapping(s, s->adj, s->edg, s->wght, j)) return 0;
    }
    return 1;
}

static int
is_directed_automorphism(struct saucy *s)
{
    int i, j;

    for (i = 0; i < s->ndiffs; ++i) {
        j = s->unsupp[i];
        if (!check_mapping(s, s->adj, s->edg, s->wght, j)) return 0;
        if (!check_mapping(s, s->dadj, s->dedg, s->dwght, j)) return 0;
    }
    return 1;
}

static void
add_induce(struct saucy *s, struct coloring *c, int who)
{
    if (!c->clen[who]) {
        s->sinduce[s->nsinduce++] = who;
    }
    else {
        s->ninduce[s->nninduce++] = who;
    }
    s->indmark[who] = 1;
}

static void
fix_fronts(struct coloring *c, int cf, int ff)
{
    int i, end = cf + c->clen[cf];
    for (i = ff; i <= end; ++i) {
        c->cfront[c->lab[i]] = cf;
    }
}

static void
array_indirect_sort(int *a, const int *b, int n)
{
    int h, i, j, k;

    /* Shell sort, as implemented in nauty, (C) Brendan McKay */
    j = n / 3;
    h = 1;
    do { h = 3 * h + 1; } while (h < j);

    do {
        for (i = h; i < n; ++i) {
            k = a[i];
            for (j = i; b[a[j-h]] > b[k]; ) {
                a[j] = a[j-h];
                if ((j -= h) < h) break;
            }
            a[j] = k;
        }
        h /= 3;
    } while (h > 0);
}

static int
at_terminal(struct saucy *s)
{
    return s->nsplits == s->n;
}

static void
add_diffnon(struct saucy *s, int k)
{
    /* Only add if we're in a consistent state */
    if (s->ndiffnons == -1) return;

    s->undiffnons[k] = s->ndiffnons;
    s->diffnons[s->ndiffnons++] = k;
}

static void
remove_diffnon(struct saucy *s, int k)
{
    int j;

    if (s->undiffnons[k] == -1) return;

    j = s->diffnons[--s->ndiffnons];
    s->diffnons[s->undiffnons[k]] = j;
    s->undiffnons[j] = s->undiffnons[k];

    s->undiffnons[k] = -1;
}

static void
add_diff(struct saucy *s, int k)
{
    if (!s->diffmark[k]) {
        s->diffmark[k] = 1;
        s->diffs[s->ndiffs++] = k;
        add_diffnon(s, k);
    }
}

static int
is_a_pair(struct saucy *s, int k)
{
    return s->unpairs[k] != -1;
}

static int
in_cell_range(struct coloring *c, int ff, int cf)
{
    int cb = cf + c->clen[cf];
    return cf <= ff && ff <= cb;
}

static void
add_pair(struct saucy *s, int k)
{
    if (s->npairs != -1) {
        s->unpairs[k] = s->npairs;
        s->pairs[s->npairs++] = k;
    }
}

static void
eat_pair(struct saucy *s, int k)
{
    int j;
    j = s->pairs[--s->npairs];
    s->pairs[s->unpairs[k]] = j;
    s->unpairs[j] = s->unpairs[k];
    s->unpairs[k] = -1;
}

static void
pick_all_the_pairs(struct saucy *s)
{
    int i;
    for (i = 0; i < s->npairs; ++i) {
        s->unpairs[s->pairs[i]] = -1;
    }
    s->npairs = 0;
}

static void
clear_undiffnons(struct saucy *s)
{
    int i;
    for (i = 0 ; i < s->ndiffnons ; ++i) {
        s->undiffnons[s->diffnons[i]] = -1;
    }
}

static void
fix_diff_singleton(struct saucy *s, int cf)
{
    int r = s->right.lab[cf];
    int l = s->left.lab[cf];
    int rcfl;

    if (!s->right.clen[cf] && r != l) {

        /* Make sure diff is marked */
        add_diff(s, r);

        /* It is now undiffed since it is singleton */
        ++s->nundiffs;
        remove_diffnon(s, r);

        /* Mark the other if not singleton already */
        rcfl = s->right.cfront[l];
        if (s->right.clen[rcfl]) {
            add_diff(s, l);

            /* Check for pairs */
            if (in_cell_range(&s->right, s->left.unlab[r], rcfl)) {
                add_pair(s, l);
            }
        }
        /* Otherwise we might be eating a pair */
        else if (is_a_pair(s, r)) {
            eat_pair(s, r);
        }
    }
}

static void
fix_diff_subtract(struct saucy *s, int cf, const int *a, const int *b)
{
    int i, k;
    int cb = cf + s->right.clen[cf];

    /* Mark the contents of the first set */
    for (i = cf; i <= cb; ++i) {
        s->stuff[a[i]] = 1;
    }

    /* Add elements from second set not present in the first */
    for (i = cf; i <= cb; ++i) {
        k = b[i];
        if (!s->stuff[k]) add_diff(s, k);
    }

    /* Clear the marks of the first set */
    for (i = cf; i <= cb; ++i) {
        s->stuff[a[i]] = 0;
    }
}

static void
fix_diffs(struct saucy *s, int cf, int ff)
{
    int min;

    /* Check for singleton cases in both cells */
    fix_diff_singleton(s, cf);
    fix_diff_singleton(s, ff);

    /* If they're both nonsingleton, do subtraction on smaller */
    if (s->right.clen[cf] && s->right.clen[ff]) {
        min = s->right.clen[cf] < s->right.clen[ff] ? cf : ff;
        fix_diff_subtract(s, min, s->left.lab, s->right.lab);
        fix_diff_subtract(s, min, s->right.lab, s->left.lab);
    }
}

static void
clear_parent(struct saucy *s, struct coloring *c)
{   
    int i;
    for (int i = 0; i < s->n; ++i)
        c->parent[i] = -1;
}

static void
split_color(struct coloring *c, int cf, int ff)
{
    int cb, fb;

    /* Fix lengths */
    fb = ff - 1;
    cb = cf + c->clen[cf];
    c->clen[cf] = fb - cf;
    c->clen[ff] = cb - ff;
    (c->parent[cf]) > -1 ? c->parent[ff] = c->parent[cf] : c->parent[ff] = cf; 

    /* Fix cell front pointers */
    fix_fronts(c, ff, ff);
}

static void
split_common(struct saucy *s, struct coloring *c, int cf, int ff)
{
    split_color(c, cf, ff);

    /* Add to refinement */
    if (s->indmark[cf] || c->clen[ff] < c->clen[cf]) {
        add_induce(s, c, ff);
    }
    else {
        add_induce(s, c, cf);
    }
}

static int
split_left(struct saucy *s, struct coloring *c, int cf, int ff)
{
    /* Record the split */
    s->splitwho[s->nsplits] = ff;
    s->splitfrom[s->nsplits] = cf;
    ++s->nsplits;

    /* Do common splitting tasks */
    split_common(s, c, cf, ff);

    /* Always succeeds */
    return 1;
}

static int
split_init(struct saucy *s, struct coloring *c, int cf, int ff)
{
    split_left(s, c, cf, ff);

    /* Maintain nonsingleton list for finding new targets */
    if (c->clen[ff]) {
        s->prevnon[s->nextnon[cf]] = ff;
        s->nextnon[ff] = s->nextnon[cf];
        s->prevnon[ff] = cf;
        s->nextnon[cf] = ff;
    }
    if (!c->clen[cf]) {
        s->nextnon[s->prevnon[cf]] = s->nextnon[cf];
        s->prevnon[s->nextnon[cf]] = s->prevnon[cf];
    }

    /* Always succeeds */
    return 1;
}

static int
split_other(struct saucy *s, struct coloring *c, int cf, int ff)
{
    int k = s->nsplits;

    /* Verify the split with init */
    if (s->splitwho[k] != ff || s->splitfrom[k] != cf
            || k >= s->splitlev[s->lev]) {
        return 0;
    }
    ++s->nsplits;

    /* Do common splitting tasks */
    split_common(s, c, cf, ff);

    /* Fix differences with init */
    fix_diffs(s, cf, ff);

    /* If we got this far we succeeded */
    return 1;
}

static int
refine_cell(struct saucy *s, struct coloring *c,
    int (*refine)(struct saucy *, struct coloring *, int))
{
    int i, cf, ret = 1;

    /*
     * The connected list must be consistent.  This is for
     * detecting mappings across nodes at a given level.  However,
     * at the root of the tree, we never have to map with another
     * node, so we lack this consistency constraint in that case.
     */
    if (s->lev > 1) introsort(s->clist, s->csize);

    /* Now iterate over the marked cells */
    for (i = 0; ret && i < s->csize; ++i) {
        cf = s->clist[i];
        ret = refine(s, c, cf);
    }

    /* Clear the connected marks */
    for (i = 0; i < s->csize; ++i) {
        cf = s->clist[i];
        s->conncnts[cf] = 0;
    }
    s->csize = 0;
    return ret;
}

static int
maybe_split(struct saucy *s, struct coloring *c, int cf, int ff)
{
    return cf == ff ? 1 : s->split(s, c, cf, ff);
}

static int
ref_single_cell(struct saucy *s, struct coloring *c, int cf)
{
    int zcnt = c->clen[cf] + 1 - s->conncnts[cf];
    return maybe_split(s, c, cf, cf + zcnt);
}

static int
ref_singleton(struct saucy *s, struct coloring *c,
    const int *adj, const int *edg, const int *wght, int cf)
{
    int i, j, k = c->lab[cf], temp;
    int wcount = 0;
    int ret = 1;

    for( i = adj[k]; i != adj[k+1]; i++ ){
        /*
        if( !list_search( s->diffL, wght[i], wcount, &ndx ) ){ 
            s->diffL[ndx] = wght[i];
            DCCOUNT(0, ndx) = 1;
            DCCOUNT(1, ndx) = i;
            ++wcount;
        }else{
            ++DCCOUNT(0, ndx);
            temp = DCCOUNT(0, ndx);
            DCCOUNT(temp, ndx) = i;
        }
        */
        //TODO: fix DCCOUNT, supposed to be 2d array, but not 2d with only one edge weight 
        if( !DCCOUNT(0, wght[i]) ){
            s->diffL[wcount] = wght[i];
            wcount++;
        }
        DCCOUNT(0, wght[i])++;
        temp = DCCOUNT(0, wght[i]);
        DCCOUNT(temp, wght[i]) = i;
    }

    for( j = 0; j < wcount && ret; j++ ){
        for( i = 1; i <= DCCOUNT(0, s->diffL[j]); i++ ){
            data_mark( s, c, edg[ DCCOUNT(i, s->diffL[j]) ] );
            DCCOUNT(i, s->diffL[j]) = 0;
        }

        ret = ret && refine_cell(s, c, ref_single_cell);
        
        DCCOUNT(0, s->diffL[j]) = 0;
        s->diffL[j] = 0;
    }

    return ret;
}

static int
ref_singleton_directed(struct saucy *s, struct coloring *c, int cf)
{
    return ref_singleton(s, c, s->adj, s->edg, s->wght, cf)
        && ref_singleton(s, c, s->dadj, s->dedg, s->dwght, cf);
}

static int
ref_singleton_undirected(struct saucy *s, struct coloring *c, int cf)
{
    return ref_singleton(s, c, s->adj, s->edg, s->wght, cf);
}

static int
ref_nonsingle_cell( struct saucy *s, struct coloring *c, int cf )
{
    int cnt, i, cb, nzf, ff, fb, bmin, bmax;

    /* Find the front and back */
    cb = cf + c->clen[cf];
    nzf = cb - s->conncnts[cf] + 1;

    /* Prepare the buckets */
    ff = nzf;
    cnt = s->ccount[c->lab[ff]];
    s->count[ff] = bmin = bmax = cnt;
    s->bucket[cnt] = 1;

    /* Iterate through the rest of the vertices */
    while (++ff <= cb) {
        cnt = s->ccount[c->lab[ff]];

        /* Initialize intermediate buckets */
        while (bmin > cnt) s->bucket[--bmin] = 0;
        while (bmax < cnt) s->bucket[++bmax] = 0;

        /* Mark this count */
        ++s->bucket[cnt];
        s->count[ff] = cnt;
    }

    /* If they all had the same count, bail */
    if (bmin == bmax && cf == nzf) return 1;
    ff = fb = nzf;

    /* Calculate bucket locations, sizes */
    for (i = bmin; i <= bmax; ++i, ff = fb) {
        if (!s->bucket[i]) continue;
        fb = ff + s->bucket[i];
        s->bucket[i] = fb;
    }

    /* Repair the partition nest */
    for (i = nzf; i <= cb; ++i) {
        s->junk[--s->bucket[s->count[i]]] = c->lab[i];
    }
    for (i = nzf; i <= cb; ++i) {
        set_label(c, i, s->junk[i]);
    }

    /* Split; induce */
    for (i = bmax; i > bmin; --i) {
        ff = s->bucket[i];
        if (ff && !s->split(s, c, cf, ff)) return 0;
    }

    /* If there was a zero area, then there's one more cell */
    return maybe_split(s, c, cf, s->bucket[bmin]);
}

static int
ref_nonsingle(struct saucy *s, struct coloring *c,
    const int *adj, const int *edg, const int *wght, int cf)
{
    int i, j, k, temp;
    int wcount = 0;
    int ret = 1;
    const int cb = cf + c->clen[cf];
    const int size = cb - cf + 1;

    /* Double check for nonsingles which became singles later */
    if (cf == cb) {
        return ref_singleton(s, c, adj, edg, wght, cf);
    }

    /* Establish connected list */
    /* junk holds the cell used for refining the others */
    /* this cell doesn't change in the course of refinement */
    memcpy(s->junk, c->lab + cf, size * sizeof(int));
    for( i = 0; i < size; ++i ){
        k = s->junk[i];
        for( j = adj[k]; j != adj[k+1]; ++j ){
            if( !DCCOUNT(0, wght[j]) ){
                s->diffL[wcount] = wght[j];
                ++wcount;
            }
            ++DCCOUNT(0, wght[j]);
            temp = DCCOUNT(0, wght[j]);
            DCCOUNT(temp, wght[j]) = j;
            
            /*
            if( !list_search( s->diffL, wght[j], wcount, &ndx ) ){
                s->diffL[ndx] = wght[j];
                DCCOUNT(0, ndx) = 1;
                DCCOUNT(1, ndx) = j;
                ++wcount; 
            }else{
                ++DCCOUNT(0, ndx);
                temp = DCCOUNT(0, ndx);
                DCCOUNT(temp, ndx) = j;
            }
            */
        }
    }

    for( j = 0; j < wcount && ret; ++j ){
        for( i = 1; i <= DCCOUNT(0, s->diffL[j]); ++i ){
            data_count( s, c, edg[ DCCOUNT(i, s->diffL[j]) ] );
        }

        /* Refine the cells we're connected to with weight j */
        ret = ret && refine_cell( s, c, ref_nonsingle_cell );

        /* Clear the counts, including weight counts */
        for( i = 1; i <= DCCOUNT(0, s->diffL[j]); ++i ){
            s->ccount[edg[ DCCOUNT(i, s->diffL[j]) ]] = 0;
            DCCOUNT(i, s->diffL[j]) = 0;
        }
        DCCOUNT(0, s->diffL[j]) = 0;
        s->diffL[j] = 0;
    }

    return ret;
}

static int
ref_nonsingle_directed(struct saucy *s, struct coloring *c, int cf)
{
    /* added weight data as a parameter - Schrock */
    return ref_nonsingle(s, c, s->adj, s->edg, s->wght, cf)
        && ref_nonsingle(s, c, s->dadj, s->dedg, s->dwght, cf);
}

static int
ref_nonsingle_undirected(struct saucy *s, struct coloring *c, int cf)
{
    /* added weight data as a parameter - Schrock */
    return ref_nonsingle(s, c, s->adj, s->edg, s->wght, cf);
}

static void
clear_refine(struct saucy *s)
{
    int i;
    for (i = 0; i < s->nninduce; ++i) {
        s->indmark[s->ninduce[i]] = 0;
    }
    for (i = 0; i < s->nsinduce; ++i) {
        s->indmark[s->sinduce[i]] = 0;
    }
    s->nninduce = s->nsinduce = 0;
}

static int
refine(struct saucy *s, struct coloring *c)
{
    int front;

    /* Keep going until refinement stops */
    while (1) {

        /* If discrete, bail */
        if (at_terminal(s)) {
            clear_refine(s);
            return 1;
        };

        /* Look for something else to refine on */
        if (s->nsinduce) {
            front = s->sinduce[--s->nsinduce];
            s->indmark[front] = 0;
            if (!s->ref_singleton(s, c, front)) break;
        }
        else if (s->nninduce) {
            front = s->ninduce[--s->nninduce];
            s->indmark[front] = 0;
            if (!s->ref_nonsingle(s, c, front)) break;
        }
        else {
            return 1;
        };
    }

    clear_refine(s);
    return 0;
}

static int
descend(struct saucy *s, struct coloring *c, int target, int min, struct eq_part *eq_steps, int idx)
{
    int back = target + c->clen[target];
    int ret;

    /* Count this node */
    ++s->stats->nodes;

    /* Move the minimum label to the back */
    eq_steps[idx].target = target;
    swap_labels(c, min, back);

    /* Split the cell */
    s->difflev[s->lev] = s->ndiffs;
    s->undifflev[s->lev] = s->nundiffs;
    ++s->lev;
    s->split(s, c, target, back);

    /* Now go and do some work */
    ret = refine(s, c);

    /* This is the new enhancement in saucy 3.0 */
    if (c == &s->right && ret) {
            int i, j, v, sum1, sum2, xor1, xor2;
        for (i = s->nsplits - 1; i > s->splitlev[s->lev-1]; --i) {
            v = c->lab[s->splitwho[i]];
            sum1 = xor1 = 0;
            for (j = s->adj[v]; j < s->adj[v+1]; j++) {
                sum1 += c->cfront[s->edg[j]];
                xor1 ^= c->cfront[s->edg[j]];
            }
            v = s->left.lab[s->splitwho[i]];
            sum2 = xor2 = 0;
            for (j = s->adj[v]; j < s->adj[v+1]; j++) {
                sum2 += s->left.cfront[s->edg[j]];
                xor2 ^= s->left.cfront[s->edg[j]];
            }
            if ((sum1 != sum2) || (xor1 != xor2)) {
                ret = 0;
                break;
            }
            v = c->lab[s->splitfrom[i]];
            sum1 = xor1 = 0;
            for (j = s->adj[v]; j < s->adj[v+1]; j++) {
                sum1 += c->cfront[s->edg[j]];
                xor1 ^= c->cfront[s->edg[j]];
            }
            v = s->left.lab[s->splitfrom[i]];
            sum2 = xor2 = 0;
            for (j = s->adj[v]; j < s->adj[v+1]; j++) {
                sum2 += s->left.cfront[s->edg[j]];
                xor2 ^= s->left.cfront[s->edg[j]];
            }
            if ((sum1 != sum2) || (xor1 != xor2)) {
                ret = 0;
                break;
            }
        }
    }

    //eq_steps[idx].target = c->lab[min];
    memcpy( eq_steps[idx].labels, c->lab, s->n*sizeof(int) );
    memcpy( eq_steps[idx].fronts, c->cfront, s->n*sizeof(int) );
    memcpy( eq_steps[idx].parents, c->parent, s->n*sizeof(int) );
    clear_parent(s, c);
    //memcpy( eq_steps[idx].lens, c->clen, s->n*sizeof(int) );

    return ret;
}

static int
descend_leftmost( struct saucy *s, struct eq_part *eq_steps )
{
    int target, min, idx=0;
    
    eq_steps[idx].target = -1;
    memcpy( eq_steps[idx].labels, s->left.lab, s->n*sizeof(int) );
    memcpy( eq_steps[idx].fronts, s->left.cfront, s->n*sizeof(int) );
    //memcpy( eq_steps[idx].lens, s->left.clen, s->n*sizeof(int) );
    
    /* Keep going until we're discrete */
    while (!at_terminal(s)) {
        idx++;
        target = s->nextnon[-1];
        min = target;
        s->start[s->lev] = target;
        s->splitlev[s->lev] = s->nsplits;
        if (!descend(s, &s->left, target, min, eq_steps, idx)) return 0;
    }
    s->splitlev[s->lev] = s->n;
    return 1;
}

void
saucy_search(
    struct saucy *s,
    const struct saucy_graph *g,
    int directed,
    const int *colors,
    saucy_consumer *consumer,
    void *arg,
    struct saucy_stats *stats,
    struct eq_part *eq_steps)
{
    int i, j, max = 0;

    /* Save client information */
    s->stats = stats;
    s->arg = arg;
    s->consumer = consumer;

    /* Save graph information */
    s->n = g->n;
    /*s->e = g->e;*/
    s->wcount = g->w;
    s->adj = g->adj;
    s->edg = g->edg;
    s->wght = g->wght;
    s->dadj = g->adj + g->n + 1;
    s->dedg = g->edg + g->e;
    s->dwght = g->wght + g->e;

    /* Polymorphism */
    if (directed) {
        s->is_automorphism = is_directed_automorphism;
        s->ref_singleton = ref_singleton_directed;
        s->ref_nonsingle = ref_nonsingle_directed;
    }
    else {
        s->is_automorphism = is_undirected_automorphism;
        s->ref_singleton = ref_singleton_undirected;
        s->ref_nonsingle = ref_nonsingle_undirected;
    }

    /* Initialize scalars */
    s->indmin = 0;
    s->lev = s->anc = 1;
    s->ndiffs = s->nundiffs = s->ndiffnons = 0;

    /* The initial orbit partition is discrete */
    for (i = 0; i < s->n; ++i) {
        s->theta[i] = i;
    }

    /* The initial permutation is the identity */
    for (i = 0; i < s->n; ++i) {
        s->gamma[i] = i;
    }

    /* Initially every cell of theta has one element */
    for (i = 0; i < s->n; ++i) {
        s->thsize[i] = 1;
    }

    /* Every theta rep list is singleton */
    for (i = 0; i < s->n; ++i) {
        s->thprev[i] = s->thnext[i] = i;
    }

    /* We have no pairs yet */
    s->npairs = 0;
    for (i = 0; i < s->n; ++i) {
        s->unpairs[i] = -1;
    }

    /* Ensure no stray pointers in undiffnons, which is checked by removed_diffnon() */
    for (i = 0; i < s->n; ++i) {
        s->undiffnons[i] = -1;
    }

    /* Initialize stats */
    s->stats->grpsize_base = 1.0;
    s->stats->grpsize_exp = 0;
    s->stats->nodes = 1;
    s->stats->bads = s->stats->gens = s->stats->support = 0;

    /* Prepare for refinement */
    s->nninduce = s->nsinduce = 0;
    s->csize = 0;

    /* Count cell sizes */
    for (i = 0; i < s->n; ++i) {
        s->ccount[colors[i]]++;
        if (max < colors[i]) max = colors[i];
    }
    s->nsplits = max + 1;

    /* Build cell lengths */
    s->left.clen[0] = s->ccount[0] - 1;
    for (i = 0; i < max; ++i) {
        s->left.clen[s->ccount[i]] = s->ccount[i+1] - 1;
        s->ccount[i+1] += s->ccount[i];
    }

    /* Build the label array */
    for (i = 0; i < s->n; ++i) {
        set_label(&s->left, --s->ccount[colors[i]], i);
    }

    /* Clear out ccount */
    for (i = 0; i <= max; ++i) {
        s->ccount[i] = 0;
    }

    /* Update refinement stuff based on initial partition */
    for (i = 0; i < s->n; i += s->left.clen[i]+1) {
        add_induce(s, &s->left, i);
        fix_fronts(&s->left, i, i);
    }

    /* Prepare lists based on cell lengths */
    for (i = 0, j = -1; i < s->n; i += s->left.clen[i] + 1) {
        if (!s->left.clen[i]) continue;
        s->prevnon[i] = j;
        s->nextnon[j] = i;
        j = i;
    }

    /* Fix the end */
    s->prevnon[s->n] = j;
    s->nextnon[j] = s->n;
    /* Preprocessing after initial coloring */
    
    s->split = split_init; 
    clear_parent(s, &s->left);
    refine(s, &s->left);

    /* Descend along the leftmost branch and compute zeta */
    descend_leftmost( s, eq_steps );
}

static int *ints(int n) { return (int *)malloc(n * sizeof(int)); }
static int *zeros(uint64_t n) { return (int *)calloc(n, sizeof(int)); }
static char *bits(int n) { return (char *)calloc(n, sizeof(char)); }

struct saucy *
saucy_alloc(int n, int w)
{
    struct saucy *s = (struct saucy *)malloc(sizeof(struct saucy));
    if (s == NULL) return NULL;

    s->ninduce = ints(n);
    s->sinduce = ints(n);
    s->indmark = bits(n);
    s->left.cfront = zeros(n);
    s->left.clen = ints(n);
    s->right.cfront = zeros(n);
    s->right.clen = ints(n);
    s->stuff = bits(n+1);
    s->wstuff = ints(n+1);
    s->bucket = ints(n+2);
    s->count = ints(n+1);
    s->ccount = zeros(n);
    /* s->dccount = zeros(e*e); */
    /* TODO: some improvement likely available if number of 
             different weights recorded, s->dccount = zeros(w*n) */
    //s->dccount = zeros(w*n);
    s->dccount = zeros((n*n+1)*w);
    /* s->diffL = (double *)malloc(n * n * sizeof(double)); */
    /* s->diffL = (int *)calloc(e, sizeof(int)); */ /* TODO: figure out if this is right */
    /* TODO: this becomes s->diffL = zeros(w) if number of different
             weights is recorded */
    //s->diffL = zeros(n);
    s->diffL = zeros(w);
    s->clist = ints(n);
    s->nextnon = ints(n+1) + 1;
    s->prevnon = ints(n+1);
    s->anctar = ints(n);
    s->start = ints(n);
    s->gamma = ints(n);
    s->junk = ints(n);
    s->theta = ints(n);
    s->thsize = ints(n);
    s->left.lab = ints(n);
    s->left.unlab = ints(n);
    s->left.parent = ints(n);
    s->right.lab = ints(n);
    s->right.unlab = ints(n);
    s->splitwho = ints(n);
    s->splitfrom = ints(n);
    s->splitlev = ints(n+1);
    s->unsupp = ints(n);
    s->conncnts = zeros(n);
    s->diffmark = bits(n);
    s->diffs = ints(n);
    s->difflev = ints(n);
    s->undifflev = ints(n);
    s->specmin = ints(n);
    s->thnext = ints(n);
    s->thprev = ints(n);
    s->threp = ints(n);
    s->thfront = ints(n);
    s->pairs = ints(n);
    s->unpairs = ints(n);
    s->diffnons = ints(n);
    s->undiffnons = ints(n);

    if (s->ninduce && s->sinduce && s->left.cfront && s->left.clen
        && s->right.cfront && s->right.clen
        && s->stuff && s->bucket && s->count && s->ccount
        && s->wstuff && s->diffL && s->dccount
        && s->clist && s->nextnon-1 && s->prevnon
        && s->start && s->gamma && s->theta && s->left.unlab
        && s->right.lab && s->right.unlab
        && s->left.lab && s->splitwho && s->junk
        && s->splitfrom && s->splitlev && s->thsize
        && s->unsupp && s->conncnts && s->anctar
        && s->diffmark && s->diffs && s->indmark
        && s->thnext && s->thprev && s->threp && s->thfront
        && s->pairs && s->unpairs && s->diffnons && s->undiffnons
        && s->difflev && s->undifflev && s->specmin)
    {
        return s;
    }
    else {
        saucy_free(s);
        return NULL;
    }
}

void
saucy_free(struct saucy *s)
{
    free(s->undiffnons);
    free(s->diffnons);
    free(s->unpairs);
    free(s->pairs);
    free(s->thfront);
    free(s->threp);
    free(s->thnext);
    free(s->thprev);
    free(s->specmin);
    free(s->anctar);
    free(s->thsize);
    free(s->undifflev);
    free(s->difflev);
    free(s->diffs);
    free(s->diffmark);
    free(s->conncnts);
    free(s->unsupp);
    free(s->splitlev);
    free(s->splitfrom);
    free(s->splitwho);
    free(s->right.unlab);
    free(s->right.lab);
    free(s->left.unlab);
    free(s->left.lab);
    free(s->theta);
    free(s->junk);
    free(s->gamma);
    free(s->start);
    free(s->prevnon);
    free(s->nextnon-1);
    free(s->clist);
    free(s->ccount);
    free(s->dccount);
    free(s->count);
    free(s->diffL);
    free(s->bucket);
    free(s->stuff);
    free(s->wstuff);
    free(s->right.clen);
    free(s->right.cfront);
    free(s->left.clen);
    free(s->left.cfront);
    free(s->indmark);
    free(s->sinduce);
    free(s->ninduce);
    free(s);
}
