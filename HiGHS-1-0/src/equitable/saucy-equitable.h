#ifndef SAUCY_H
#define SAUCY_H

#define SAUCY_VERSION "2.0"
#include <limits>
#include "HighsSimpleDec.h"

typedef int saucy_consumer(int, const int *, int, int *, void *);

struct saucy;

struct saucy_stats {
    double grpsize_base;
    int grpsize_exp;
    int levels;
    int nodes;
    int bads;
    int gens;
    int support;
    int iter;
};

struct saucy_graph {
    int n;
    int nCols;
    int e;
    int w;
    int *adj;
    int *edg;
    int *wght;
};

struct eq_part {
    int nsplits;
    int target;
    int *labels;
    int *fronts;
    int *parents;
};

struct saucy *saucy_alloc(long long int n, long long int w, std::string model_name);

void saucy_search(
    struct saucy *s,
    const struct saucy_graph *graph,
    int directed,
    const int *colors,
    saucy_consumer *consumer,
    void *arg,
    struct saucy_stats *stats,
    struct eq_part *eq_steps);

void saucy_free(struct saucy *s);

#endif
