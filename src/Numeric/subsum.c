
/* ======================================================================
	         SUBSUM.c David Pisinger 1994,1995, revised 1999
   ====================================================================== */

/* This C-code solves the subset sum problem. It was presented as
 * subroutine in 
 *
 *   D. Pisinger
 *   An exact algorithm for large multiple knapsack problems
 *   to appear European Journal of Operational Research
 *
 * where it was used to tighten the capacity constraints and to 
 * derive lower bounds by splitting a solution to the surrogate 
 * relaxed problem.
 * Further details on the project can also be found in
 *
 *   D. Pisinger
 *   Algorithms for Knapsack Problems
 *   Report 95/1, DIKU, University of Copenhagen
 *   Universitetsparken 1
 *   DK-2100 Copenhagen
 *
 * The algorithm may be used for academic, non-commercial purposes
 * only.
 * -------------------------------------------------------------------
 * The present code is a callable routine which solves a Subset-sum
 * Problem:
 *
 *           maximize   \sum_{j=1}^{n} w_{j} x_{j}
 *           subject to \sum_{j=1}^{n} w_{j} x_{j} \leq c
 *                      x_{j} \in \{0,1\}, j = 1,\ldots,n
 *
 * The decomp algorithm is called as
 *
 *          z = decomp(n, w, x, c)
 *
 * where w[], x[] are arrays of integers. The optimal objective
 * value is returned in z, and x[] gives the solution vector.
 * If you need a different interface for your algorithm, decomp
 * may easily be adapted to your own datastructures since all tables
 * are copied to the internal representation.
 *
 * Different types should be defined as follows:
 *
 *    itype     should be sufficiently large to hold the weights
 *    stype     should be sufficient to hold sum of weights
 *    ptype     should hold the product of an stype and itype
 *
 * The code has been tested on a hp9000/735, and conforms with the
 * ANSI-C standard.
 *
 * Errors and questions are refered to:
 *
 *   David Pisinger, associate professor
 *   DIKU, University of Copenhagen,
 *   Universitetsparken 1,
 *   DK-2100 Copenhagen.
 *   e-mail: pisinger@diku.dk
 *   fax: +45 35 32 14 01
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <limits.h>
#include <string.h>
#include <math.h>


/* ======================================================================
				   macros
   ====================================================================== */

#define TRUE  1
#define FALSE 0

#define SWAP(a,b)    {register item t;t=*(a);*(a)=*(b);*(b)=t;}
#define SIZE(a)      ((int) (((a)->lset+1)-(a)->fset))


/* ======================================================================
				 type declarations
   ====================================================================== */

typedef int   boolean; /* boolean                 */
typedef long  ntype;   /* number of items         */
typedef long  itype;   /* item weights            */
typedef long  stype;    /* sum of profit or weight */

/* item record */
typedef struct irec {
  itype     w;         /* weight of item */
  int       *x;        /* solution variable */
} item;

/* state record */
typedef struct srec {
  stype     wsum;      /* weight sum of state */
  item      *mod;      /* by modifying mod */
} state;

/* set of states */
typedef struct {
  ntype    size;       /* set size */
  state    *fset;      /* first element in set */
  state    *lset;      /* last element in set */
} wgtset;


/* all problem specific information */
typedef struct {
  ntype    n;          /* number of items */
  itype    type;       /* type of problem */

  item     *fitem;     /* first item */
  item     *litem;     /* last item */
  stype    c;          /* capacity of knapsack */
  stype    z;          /* found solution */

  /* information corresponding to Subset-sum Problem */
  item     *s;           /* s+1 is first item enumme. */
  item     *t;           /* t-1 is last item enummer. */
  item     *b;           /* break item for smallest c */
  stype    wsumall;      /* weight sum of all items   */
  stype    wsumb;        /* weight sum up to b        */
  wgtset   da;           /* set of partial vectors    */
  wgtset   db;           /* set of partial vectors    */
  state    *opt1, *opt2; /* states giving current z   */
} allinfo;


/* ======================================================================
                                  error
   ====================================================================== */

void error(char *str, ...)
{
  va_list args;

  va_start(args, str);
  vprintf(str, args); printf("\n");
  va_end(args);
  printf("Program is terminated !!!\n\n");
  exit(-1);
}


/* ======================================================================
				  palloc
   ====================================================================== */

void pfree(void *p)
{
  if (p == NULL) error("freeing null");
  free(p);
}


void * palloc(long no, long sz)
{
  long *p, size;

  size = sz * no;
  if (size == 0) size = 1;
  if (size != (size_t) size) error("alloc too big %ld", size);
  p = malloc(size);
  if (p == NULL) error("no memory size %ld", size);
  return p;
}


/* ======================================================================
                                  wmultiply
   ====================================================================== */

void wadd(allinfo *a, item *t)
{
  register state *i, *j, *k;
  register stype c, w, iw, jw, kw;
  register item *mod;
  state *r;
 
  /* initialize limits */
  r = palloc(2*a->db.size + 1, sizeof(state));
  i = a->db.fset;
  j = a->db.fset;
  k = r;
  w = t->w;
  c = a->c;
  mod = t;

  /* add large state at end, and copy first state */
  (a->db.lset+1)->wsum = c + 1;
  *k = *i; i++;

  /* now merge sets: i, j indices to each set, k index to merged set */
  for (iw = i->wsum, jw = j->wsum + w, kw = k->wsum;;) {
    if (iw <= jw) {
      if (iw > c) break;
      k++; kw = iw; *k = *i;
      i++; iw = i->wsum;
    } else {
      if (jw != kw) { k++; kw = k->wsum = jw; k->mod = mod; }
      j++; jw = j->wsum + w;
    }
  }

  /* save limits */
  pfree(a->db.fset);
  a->db.fset = r;
  a->db.lset = k;
  a->db.size = SIZE(&(a->db));
}


void wsub(allinfo *a, item *s)
{
  register state *i, *j, *k;
  register stype c, w, iw, jw, kw;
  register item *mod;
  state *r;

  /* initialize limits */
  r = palloc(2*a->da.size + 1, sizeof(state));
  i = a->da.fset;
  j = a->da.fset;
  k = r;
  w = s->w;
  /* c = a->c - (a->wsumall - a->wsumb); */
  c = a->z - (a->wsumall - a->wsumb); 
  if (c < 0) c = 0; /* never subtract more than can improve solution */
  mod = s;

  /* add small state at end, and copy first state */
  (a->da.lset+1)->wsum = c-1;
  *k = *i; i++;

  /* now merge sets: i, j indices to each set, k index to merged set */
  for (iw = i->wsum, jw = j->wsum - w, kw = k->wsum;;) {
    if (iw >= jw) {
      if (iw < c) break;
      k++; kw = iw; *k = *i;
      i++; iw = i->wsum;
    } else {
      if (jw != kw) { k++; kw = k->wsum = jw; k->mod = mod; }
      j++; jw = j->wsum - w;
    }
  }

  /* save limits */
  pfree(a->da.fset);
  a->da.fset = r;
  a->da.lset = k;
  a->da.size = SIZE(&(a->da));
}


/* ======================================================================
			      definesolution
   ====================================================================== */

void definesolution(allinfo *a)
{
  register item *i, *m, *l;
  register state *j, *k;
  register stype ps, ws, w0, c;

  /* first take items up to b2 */
  for (i = a->fitem, m = a->b; i != m; i++) *(i->x) = 1;

  /* backtrack set da */
  for (j = a->opt1, ws = a->opt1->wsum, w0 = a->wsumb; ws != w0; j--) {
    if (j->wsum == ws) { i = j->mod; *(i->x) = 0; ws += i->w; }
  }

  /* backtrack set db */
  for (j = a->opt2, ws = a->opt2->wsum; ws != 0; j--) {
    if (j->wsum == ws) { i = j->mod; *(i->x) = 1; ws -= i->w; }
  }
}


/* ======================================================================
				  reducewgtset
   ====================================================================== */

void reducewgtset(allinfo *a)
{
  register state *i, *j, *k, *l;
  register stype ws, c, maxw;
  state *i1, *k1;

  /* check if knapsack is filled */
  c = a->c;

  i = a->da.fset; j = a->da.lset+1;
  k = a->db.fset; l = a->db.lset+1;
  for (maxw = -1, l->wsum = c+1; i != j; i++) {
    ws = i->wsum + k->wsum;
    while (ws <= c) {
      if (ws > maxw) { maxw = ws; i1 = i; k1 = k; }
      k++; if (k > l) error("bad value of k");
      ws = i->wsum + k->wsum;
    }
  }
  a->z = maxw; a->opt1 = i1; a->opt2 = k1;
}


/* ======================================================================
				 findbreak
   ====================================================================== */

void findbreak(allinfo *a)
{
  register item *i, *m;
  register stype r, wtot;
  stype wsum;
  state *k;

  /* find break item for the current knapsack */
  r = a->c;
  for (i = a->fitem, m = a->litem+1; i != m; i++) {
    if (i->w > r) break;
    r -= i->w;
  }
  wsum = a->c - r;

  /* initialize limits */
  a->b     = i;
  a->s     = i-1;
  a->t     = i;
  a->wsumb = wsum;

  /* find total weight sum */
  for (i = a->b, m = a->litem+1, wtot = wsum; i != m; i++) wtot += i->w;
  a->wsumall = wtot; 

  /* initialize table of partial vectors */
  a->da.size = 1;
  a->da.fset = palloc(a->da.size + 1, sizeof(state));
  a->da.lset = a->da.fset;

  a->db.size = 1;
  a->db.fset = palloc(a->db.size + 1, sizeof(state));
  a->db.lset = a->db.fset;

  /* initialize first partial vector */
  k = a->da.fset; k->wsum = wsum; k->mod = NULL;
  k = a->db.fset; k->wsum = 0;    k->mod = NULL;
}


/* ======================================================================
                                  copyproblem
   ====================================================================== */

void copyproblem(item *f, item *l, int *w, int *x)
{
  register item *i, *m;
  register int *ww, *xx;

  for (i = f, m = l+1, ww = w, xx = x; i != m; i++, ww++, xx++) {
    *xx = 0; i->w = *ww; i->x = xx;
  }
}


/* ======================================================================
				  partition
   ====================================================================== */

void partition(allinfo *a)
{
  stype c;
  item *i;

  findbreak(a);
  for (;;) {
    reducewgtset(a);
    if (a->z == a->c) break;
    if ((a->s < a->fitem) && (a->t > a->litem)) break;

    i = a->t;
    if (i <= a->litem) { wadd(a,i); a->t++; }

    i = a->s;
    if (i >= a->fitem) { wsub(a,i); a->s--; }
  }
  definesolution(a);
  pfree(a->da.fset);
  pfree(a->db.fset);
}


/* ======================================================================
				decomp
   ====================================================================== */

int decomp(int n, int *w, int *x, int c)
{
  allinfo a;
  item *tab;
  item *i;

  /* allocate space for internal representation */
  tab = (item *) palloc(n, sizeof(item));
  a.fitem = &tab[0]; a.litem = &tab[n-1];
  copyproblem(a.fitem, a.litem, w, x);

  /* copy information */
  a.c = c;
  a.z = -1;
  a.n = n;

  /* run horowitz-sahni decomposition */
  partition(&a);

  /* return objective value */
  pfree(tab);
  return a.z;
}


/* ======================================================================
                                end
   ====================================================================== */

