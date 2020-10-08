#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rstat.h>

/// Simplified version of gsl_rstat_workspace: for mean and variance only
typedef struct
{
  double mean;     /* current mean */
  double M2;       /* M_k = sum_{i=1..n} [ x_i - mean_n ]^k */
  double M3;
  double M4;
  size_t n;        /* number of data points added */
} mv_workspace;

char global_abort = 0;

double active_x;
gsl_rstat_workspace *active_rstat_p = NULL;

/// Custom error handler that attempts to dump state but does not tell the process to exit or otherwise interfere.
void handler(const char * reason, const char * file, int line, int gsl_errno) {
  global_abort = 1;
  printf("GSL error handler: %s, %s:%d, %s (%u)\n", reason, file, line, gsl_strerror(gsl_errno), gsl_errno);
  printf("active_x = %.8e, active_rstat_p = %p\n", active_x, active_rstat_p);
  if (active_rstat_p) {
    gsl_rstat_workspace *p1 = active_rstat_p;
    printf("min=%.6e, max=%.6e, mean=%.6e, M2=%.6e, M3=%.6e, M4=%.6e, n=%lu\n",
            p1->min, p1->max, p1->mean, p1->M2, p1->M3, p1->M4, p1->n
    );

    gsl_rstat_quantile_workspace *p2 = p1->median_workspace_p;
    printf("p=%.6e, n=%lu\n", p2->p, p2->n);
    // q[0], p->q[1], p->q[2], p->q[3], p->q[4]
    for (size_t i=0; i < 5; i++) {
      printf("[%lu] q=%.6e, npos=%d, np=%.6e, dnp=%.6e\n",
              i, p2->q[i], p2->npos[i], p2->np[i], p2->dnp[i]
      );
    }
    fflush(stdout);
  }
}

/// Simplified gsl_rstat_add (stripped code not needed for mean/variance)
int mv_add(const double x, mv_workspace *w) {
  double delta = x - w->mean;
  double delta_n, delta_nsq, term1, n;

  n = (double) ++(w->n);
  delta_n = delta / n;
  delta_nsq = delta_n * delta_n;
  term1 = delta * delta_n * (n - 1.0);
  w->mean += delta_n;
  w->M4 += term1 * delta_nsq * (n * n - 3.0 * n + 3.0) +
          6.0 * delta_nsq * w->M2 - 4.0 * delta_n * w->M3;
  w->M3 += term1 * delta_n * (n - 2.0) - 3.0 * delta_n * w->M2;
  w->M2 += term1;

  return GSL_SUCCESS;
}

/// Taken straight from rstat.c -- just changed struct
double mv_variance(mv_workspace *w) {
  if (w->n > 1) {
    double n = (double) w->n;
    return (w->M2 / (n - 1.0));
  }
  return 0.0;
}


/// Wrapper functions that apply the input value to each implementation.
/// Add also updates the global pointer to the active gsl_rstat_workspace
/// so that it can be accessed by the GSL error handler for debugging.
int add(double x, gsl_rstat_workspace *rstat_p, mv_workspace *mv_p) {
  active_rstat_p = rstat_p;
  active_x = x;
  mv_add(x, mv_p);
  int retval = gsl_rstat_add(x, rstat_p);
  if (retval != GSL_SUCCESS) {
    printf("gsl_rstat_add returned failure status: %d\n", retval);
    fflush(stdout);
  }
  return retval;
}

void print_stats(size_t i, gsl_rstat_workspace *rstat_p, mv_workspace *mv_p) {
  printf("i=%lu, gsl_mean=%.8e, gsl_var=%.8e, mv_mean=%.8e, mv_var=%.8e\n",
          i, gsl_rstat_mean(rstat_p), gsl_rstat_variance(rstat_p), mv_p->mean, mv_variance(mv_p));
}


/// Alternating huge/infinite values as inputs apparently makes gsl_rstat_quantile_add
/// blow up on rquartile.c:136 with GSL_ERROR ("invalid input argument x", GSL_EINVAL);
void death1() {
  gsl_rstat_workspace *rstat_p = gsl_rstat_alloc();
  mv_workspace mv = { 0 };
  double next = INFINITY;
  for (size_t i=0; i < 10; i++) {
    if (i == 9) {
      printf("Goodbye, cruel world!\n");
      fflush(stdout);
    }
    add(next, rstat_p, &mv);
    print_stats(i, rstat_p, &mv);
    next = -next;
  }
  gsl_rstat_free(rstat_p);
}


void death2() {
  gsl_rstat_workspace *rstat_p = gsl_rstat_alloc();
  mv_workspace mv = { 0 };

  global_abort = 0;
  const size_t N = 16 * 1024 * 1024;
  const double xs[] = {
    1e308, 1.2e308, 1.5e308, 1.6e308, 1.7e308, 1.7e+308,
    1e308, 1.2e308, 1.5e308, 1.6e308, 1.7e308, 1.7e+308,
    -1e308, -1.2e308, -1.5e308, -1.6e308, -1.7e308, -1.7e+308,
    //8.99E+307, -8.99E+307,
    //(double)(18428729675200000000UL - 18328729675200000000UL)
    // double range is 2.3E-308 to 1.7E+308
    //HUGE_VALF, -HUGE_VALF,
    //-INFINITY, INFINITY,
    //0, -1e9, +1e9, -5e100, +5e100,
  };
  const size_t xs_len = sizeof(xs)/sizeof(xs[0]);

  size_t j = 0;
  for (size_t i = 0; i < N; ++i) {
    if (i < 25 || i % 8192 == 0) {
      print_stats(i, rstat_p, &mv);
      //printf("%lu %.8e %.8e\n", i+1, xs[j], gsl_rstat_mean(rstat_p));
    }
    add(xs[j], rstat_p, &mv);
    if (global_abort) {
      printf("Stopping at i=%lu since error handler set abort flag.\n", i);
      break;
    }
    if (++j >= xs_len) {
      j = 0;
    }
  }
  gsl_rstat_free(rstat_p);
}

int main()
{
  gsl_set_error_handler(&handler);
  death1();
  death2();
  return 0;
}
