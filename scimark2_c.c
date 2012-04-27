// SciMark 2.0 benchmark merged in a single C file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define VERSION 2.0

#ifndef NULL
#define NULL 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


typedef struct {
    int m[17];
    int seed;
    int i;                      /* originally = 4 */
    int j;                      /* originally =  16 */
    int /*boolean */ haveRange; /* = false; */
    double left;                /* = 0.0; */
    double right;               /* = 1.0; */
    double width;               /* = 1.0; */
} Random_struct, *Random;

Random new_Random_seed(int seed);
double Random_nextDouble(Random R);
void Random_delete(Random R);
double *RandomVector(int N, Random R);
double **RandomMatrix(int M, int N, Random R);

double LU_num_flops(int N) {
    //  rougly 2/3*N^3 */
    double Nd = (double)N;
    return (2.0 * Nd * Nd * Nd / 3.0);
}

void LU_copy_matrix(int M, int N, double **lu, double **A) {
    int i, j;

    for (i = 0; i < M; i++)
        for (j = 0; j < N; j++)
            lu[i][j] = A[i][j];
}

int LU_factor(int M, int N, double **A, int *pivot) {
    int minMN = M < N ? M : N;
    int j = 0;

    for (j = 0; j < minMN; j++) {
        // find pivot in column j and  test for singularity.
        int jp = j;
        int i;
        double t = fabs(A[j][j]);
        for (i = j + 1; i < M; i++) {
            double ab = fabs(A[i][j]);
            if (ab > t) {
                jp = i;
                t = ab;
            }
        }

        pivot[j] = jp;

        /* jp now has the index of maximum element  */
        /* of column j, below the diagonal          */
        if (A[jp][j] == 0)
            return 1;           /* factorization failed because of zero pivot */


        if (jp != j) {
            /* swap rows j and jp */
            double *tA = A[j];
            A[j] = A[jp];
            A[jp] = tA;
        }

        if (j < M - 1) {        /* compute elements j+1:M of jth column  */
            /* note A(j,j), was A(jp,p) previously which was */
            /* guarranteed not to be zero (Label #1)         */

            double recp = 1.0 / A[j][j];
            int k;
            for (k = j + 1; k < M; k++)
                A[k][j] *= recp;
        }

        if (j < minMN - 1) {
            /* rank-1 update to trailing submatrix:   E = E - x*y; */
            /* E is the region A(j+1:M, j+1:N) */
            /* x is the column vector A(j+1:M,j) */
            /* y is row vector A(j,j+1:N)        */

            int ii;
            for (ii = j + 1; ii < M; ii++) {
                double *Aii = A[ii];
                double *Aj = A[j];
                double AiiJ = Aii[j];
                int jj;
                for (jj = j + 1; jj < N; jj++)
                    Aii[jj] -= AiiJ * Aj[jj];
            }
        }
    }
    return 0;
}

#define PI  3.1415926535897932

/*-----------------------------------------------------------------------*/

static int int_log2(int n);

double FFT_num_flops(int N) {
    double Nd = (double)N;
    double logN = (double)int_log2(N);
    return (5.0 * Nd - 2) * logN + 2 * (Nd + 1);
}

static int int_log2(int n) {
    int k = 1;
    int log = 0;
    for ( /*k=1 */ ; k < n; k *= 2, log++);
    if (n != (1 << log)) {
        printf("FFT: Data length is not a power of 2!: %d ", n);
        exit(1);
    }
    return log;
}

void FFT_bitreverse(int N, double *data) {
    /* This is the Goldrader bit-reversal algorithm */
    int n = N / 2;
    int nm1 = n - 1;
    int i = 0;
    int j = 0;
    for (; i < nm1; i++) {

        /*int ii = 2*i; */
        int ii = i << 1;

        /*int jj = 2*j; */
        int jj = j << 1;

        /* int k = n / 2 ; */
        int k = n >> 1;

        if (i < j) {
            double tmp_real = data[ii];
            double tmp_imag = data[ii + 1];
            data[ii] = data[jj];
            data[ii + 1] = data[jj + 1];
            data[jj] = tmp_real;
            data[jj + 1] = tmp_imag;
        }

        while (k <= j) {
            /*j = j - k ; */
            j -= k;

            /*k = k / 2 ;  */
            k >>= 1;
        }
        j += k;
    }
}


static void FFT_transform_internal(int N, double *data, int direction) {
    int n = N / 2;
    int bit = 0;
    int logn;
    int dual = 1;

    if (n == 1)
        return;                 /* Identity operation! */
    logn = int_log2(n);

    if (N == 0)
        return;

    /* bit reverse the input data for decimation in time algorithm */
    FFT_bitreverse(N, data);

    /* apply fft recursion */
    /* this loop executed int_log2(N) times */
    for (bit = 0; bit < logn; bit++, dual *= 2) {
        double w_real = 1.0;
        double w_imag = 0.0;
        int a;
        int b;

        double theta = 2.0 * direction * PI / (2.0 * (double)dual);
        double s = sin(theta);
        double t = sin(theta / 2.0);
        double s2 = 2.0 * t * t;

        for (a = 0, b = 0; b < n; b += 2 * dual) {
            int i = 2 * b;
            int j = 2 * (b + dual);

            double wd_real = data[j];
            double wd_imag = data[j + 1];

            data[j] = data[i] - wd_real;
            data[j + 1] = data[i + 1] - wd_imag;
            data[i] += wd_real;
            data[i + 1] += wd_imag;
        }

        /* a = 1 .. (dual-1) */
        for (a = 1; a < dual; a++) {
            /* trignometric recurrence for w-> exp(i theta) w */
            {
                double tmp_real = w_real - s * w_imag - s2 * w_real;
                double tmp_imag = w_imag + s * w_real - s2 * w_imag;
                w_real = tmp_real;
                w_imag = tmp_imag;
            }
            for (b = 0; b < n; b += 2 * dual) {
                int i = 2 * (b + a);
                int j = 2 * (b + a + dual);

                double z1_real = data[j];
                double z1_imag = data[j + 1];

                double wd_real = w_real * z1_real - w_imag * z1_imag;
                double wd_imag = w_real * z1_imag + w_imag * z1_real;

                data[j] = data[i] - wd_real;
                data[j + 1] = data[i + 1] - wd_imag;
                data[i] += wd_real;
                data[i + 1] += wd_imag;
            }
        }
    }
}

void FFT_transform(int N, double *data) {
    FFT_transform_internal(N, data, -1);
}

void FFT_inverse(int N, double *data) {
    int n = N / 2;
    double norm = 0.0;
    int i = 0;
    FFT_transform_internal(N, data, +1);

    /* Normalize */
    norm = 1 / ((double)n);
    for (i = 0; i < N; i++)
        data[i] *= norm;
}

double SOR_num_flops(int M, int N, int num_iterations) {
    double Md = (double)M;
    double Nd = (double)N;
    double num_iterD = (double)num_iterations;

    return (Md - 1) * (Nd - 1) * num_iterD * 6.0;
}

void SOR_execute(int M, int N, double omega, double **G, int num_iterations) {
    double omega_over_four = omega * 0.25;
    double one_minus_omega = 1.0 - omega;

    /* update interior points */
    int Mm1 = M - 1;
    int Nm1 = N - 1;
    int p, i, j;
    double *Gi;
    double *Gim1;
    double *Gip1;

    for (p = 0; p < num_iterations; p++) {
        for (i = 1; i < Mm1; i++) {
            Gi = G[i];
            Gim1 = G[i - 1];
            Gip1 = G[i + 1];
            for (j = 1; j < Nm1; j++)
                Gi[j] = omega_over_four * (Gim1[j] + Gip1[j] + Gi[j - 1] + Gi[j + 1]) +
                        one_minus_omega * Gi[j];
        }
    }
}


/**
 Estimate Pi by approximating the area of a circle.

 How: generate N random numbers in the unit square, (0,0) to (1,1)
 and see how are within a radius of 1 or less, i.e.
 <pre>

 sqrt(x^2 + y^2) < r

 </pre>
  since the radius is 1.0, we can square both sides
  and avoid a sqrt() computation:
  <pre>

    x^2 + y^2 <= 1.0

  </pre>
  this area under the curve is (Pi * r^2)/ 4.0,
  and the area of the unit of square is 1.0,
  so Pi can be approximated by
  <pre>
                # points with x^2+y^2 < 1
     Pi =~      --------------------------  * 4.0
                     total # points
  </pre>
*/
static const int SEED = 113;

double MonteCarlo_num_flops(int Num_samples) {
    /* 3 flops in x^2+y^2 and 1 flop in random routine */
    return ((double)Num_samples) * 4.0;
}

double MonteCarlo_integrate(int Num_samples) {
    Random R = new_Random_seed(SEED);
    int under_curve = 0;
    int count;

    for (count = 0; count < Num_samples; count++) {
        double x = Random_nextDouble(R);
        double y = Random_nextDouble(R);
        if (x * x + y * y <= 1.0)
            under_curve++;
    }
    Random_delete(R);
    return ((double)under_curve / Num_samples) * 4.0;
}


/* static const int mdig = 32; */
#define MDIG 32

  /* static const int one = 1; */
#define ONE 1

static const int m1 = (ONE << (MDIG - 2)) + ((ONE << (MDIG - 2)) - ONE);
static const int m2 = ONE << MDIG / 2;

  /* For mdig = 32 : m1 =          2147483647, m2 =      65536
     For mdig = 64 : m1 = 9223372036854775807, m2 = 4294967296
  */

                                /* move to initialize() because  */
                                /* compiler could not resolve as */
                                /*   a constant.                 */

static /*const */ double dm1;   /*  = 1.0 / (double) m1; */

/* private methods */

static void initialize(Random R, int seed);

Random new_Random_seed(int seed) {
    Random R = (Random) malloc(sizeof(Random_struct));
    initialize(R, seed);
    R->left = 0.0;
    R->right = 1.0;
    R->width = 1.0;
    R->haveRange = 0; /*false */
    return R;
}

Random new_Random(int seed, double left, double right) {
    Random R = (Random) malloc(sizeof(Random_struct));
    initialize(R, seed);
    R->left = left;
    R->right = right;
    R->width = right - left;
    R->haveRange = 1;           /* true */
    return R;
}

void Random_delete(Random R) {
    free(R);
}

/* Returns the next random number in the sequence.  */

double Random_nextDouble(Random R) {
    int k;
    int I = R->i;
    int J = R->j;
    int *m = R->m;

    k = m[I] - m[J];
    if (k < 0)
        k += m1;
    R->m[J] = k;

    if (I == 0)
        I = 16;
    else
        I--;
    R->i = I;

    if (J == 0)
        J = 16;
    else
        J--;
    R->j = J;

    if (R->haveRange)
        return R->left + dm1 * (double)k *R->width;
    else
        return dm1 * (double)k;

}


// PRIVATE METHODS

static void initialize(Random R, int seed) {
    int jseed, k0, k1, j0, j1, iloop;
    dm1 = 1.0 / (double)m1;
    R->seed = seed;

    if (seed < 0)
        seed = -seed;           /* seed = abs(seed) */
    jseed = (seed < m1 ? seed : m1);    /* jseed = min(seed, m1) */
    if (jseed % 2 == 0)
        --jseed;
    k0 = 9069 % m2;
    k1 = 9069 / m2;
    j0 = jseed % m2;
    j1 = jseed / m2;
    for (iloop = 0; iloop < 17; ++iloop) {
        jseed = j0 * k0;
        j1 = (jseed / m2 + j0 * k1 + j1 * k0) % (m2 / 2);
        j0 = jseed % m2;
        R->m[iloop] = j0 + m2 * j1;
    }
    R->i = 4;
    R->j = 16;

}

double *RandomVector(int N, Random R) {
    int i;
    double *x = (double *)malloc(sizeof(double) * N);
    for (i = 0; i < N; i++)
        x[i] = Random_nextDouble(R);
    return x;
}

double **RandomMatrix(int M, int N, Random R) {
    int i, j;

    /* allocate matrix */
    double **A = (double **)malloc(sizeof(double *) * M);

    if (A == NULL)
        return NULL;

    for (i = 0; i < M; i++) {
        A[i] = (double *)malloc(sizeof(double) * N);
        if (A[i] == NULL) {
            free(A);
            return NULL;
        }
        for (j = 0; j < N; j++)
            A[i][j] = Random_nextDouble(R);
    }
    return A;
}


typedef struct {
    int running;                /* boolean */
    double last_time;
    double total;

} *Stopwatch, Stopwatch_struct;

double seconds() {
    return ((double)clock()) / (double)CLOCKS_PER_SEC;
}

void Stopwtach_reset(Stopwatch Q) {
    Q->running = 0;             /* false */
    Q->last_time = 0.0;
    Q->total = 0.0;
}

Stopwatch new_Stopwatch(void) {
    Stopwatch S = (Stopwatch) malloc(sizeof(Stopwatch_struct));
    if (S == NULL)
        return NULL;
    Stopwtach_reset(S);
    return S;
}

void Stopwatch_delete(Stopwatch S) {
    if (S != NULL)
        free(S);
}

/* Start resets the timer to 0.0; use resume for continued total */

void Stopwatch_start(Stopwatch Q) {
    if (!(Q->running)) {
        Q->running = 1;         /* true */
        Q->total = 0.0;
        Q->last_time = seconds();
    }
}

/**
    Resume timing, after stopping.  (Does not wipe out accumulated times.)
*/

void Stopwatch_resume(Stopwatch Q) {
    if (!(Q->running)) {
        Q->last_time = seconds();
        Q->running = 1;         /*true */
    }
}

void Stopwatch_stop(Stopwatch Q) {
    if (Q->running) {
        Q->total += seconds() - Q->last_time;
        Q->running = 0;         /* false */
    }
}

double Stopwatch_read(Stopwatch Q) {
    if (Q->running) {
        double t = seconds();
        Q->total += t - Q->last_time;
        Q->last_time = t;
    }
    return Q->total;
}


    /* multiple iterations used to make kernel have roughly
       same granulairty as other Scimark kernels. */

double SparseCompRow_num_flops(int N, int nz, int num_iterations) {
    /* Note that if nz does not divide N evenly, then the
       actual number of nonzeros used is adjusted slightly.
     */
    int actual_nz = (nz / N) * N;
    return ((double)actual_nz) * 2.0 * ((double)num_iterations);
}


/* computes  a matrix-vector multiply with a sparse matrix
held in compress-row format.  If the size of the matrix
in MxN with nz nonzeros, then the val[] is the nz nonzeros,
with its ith entry in column col[i].  The integer vector row[]
is of size M+1 and row[i] points to the begining of the
ith row in col[].
*/
void SparseCompRow_matmult(int M, double *y, double *val, int *row,
                           int *col, double *x, int NUM_ITERATIONS) {
    int reps, r, i;
    for (reps = 0; reps < NUM_ITERATIONS; reps++) {
        for (r = 0; r < M; r++) {
            double sum = 0.0;
            int rowR = row[r];
            int rowRp1 = row[r + 1];
            for (i = rowR; i < rowRp1; i++)
                sum += x[col[i]] * val[i];
            y[r] = sum;
        }
    }
}





double **new_Array2D_double(int M, int N) {
    int i = 0;
    int failed = 0;

    double **A = (double **)malloc(sizeof(double *) * M);
    if (A == NULL)
        return NULL;

    for (i = 0; i < M; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
        if (A[i] == NULL) {
            failed = 1;
            break;
        }
    }

    /* if we didn't successfully allocate all rows of A      */
    /* clean up any allocated memory (i.e. go back and free  */
    /* previous rows) and return NULL                        */
    if (failed) {
        i--;
        for (; i <= 0; i--)
            free(A[i]);
        free(A);
        return NULL;
    } else
        return A;
}

void Array2D_double_delete(int M, int N, double **A) {
    int i;
    if (A == NULL)
        return;

    for (i = 0; i < M; i++)
        free(A[i]);
    free(A);
}

void Array2D_double_copy(int M, int N, double **B, double **A) {
    int remainder = N & 3;      /* N mod 4; */
    int i = 0;
    int j = 0;

    for (i = 0; i < M; i++) {
        double *Bi = B[i];
        double *Ai = A[i];
        for (j = 0; j < remainder; j++)
            Bi[j] = Ai[j];
        for (j = remainder; j < N; j += 4) {
            Bi[j] = Ai[j];
            Bi[j + 1] = Ai[j + 1];
            Bi[j + 2] = Ai[j + 2];
            Bi[j + 3] = Ai[j + 3];
        }
    }
}


double kernel_measureFFT(int N, double mintime, Random R) {
    /* initialize FFT data as complex (N real/img pairs) */

    int twoN = 2 * N;
    double *x = RandomVector(twoN, R);
    long cycles = 1;
    Stopwatch Q = new_Stopwatch();
    int i = 0;
    double result = 0.0;

    while (1) {
        Stopwatch_start(Q);
        for (i = 0; i < cycles; i++) {
            FFT_transform(twoN, x);     /* forward transform */
            FFT_inverse(twoN, x);       /* backward transform */
        }
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= mintime)
            break;

        cycles *= 2;

    }
    /* approx Mflops */

    result = FFT_num_flops(N) * cycles / Stopwatch_read(Q) * 1.0e-6;
    Stopwatch_delete(Q);
    free(x);
    return result;
}

double kernel_measureSOR(int N, double min_time, Random R) {
    double **G = RandomMatrix(N, N, R);
    double result = 0.0;

    Stopwatch Q = new_Stopwatch();
    int cycles = 1;
    while (1) {
        Stopwatch_start(Q);
        SOR_execute(N, N, 1.25, G, cycles);
        Stopwatch_stop(Q);

        if (Stopwatch_read(Q) >= min_time)
            break;

        cycles *= 2;
    }
    /* approx Mflops */

    result = SOR_num_flops(N, N, cycles) / Stopwatch_read(Q) * 1.0e-6;
    Stopwatch_delete(Q);
    Array2D_double_delete(N, N, G);
    return result;
}

double kernel_measureMonteCarlo(double min_time, Random R) {
    double result = 0.0;
    Stopwatch Q = new_Stopwatch();

    int cycles = 1;
    while (1) {
        Stopwatch_start(Q);
        MonteCarlo_integrate(cycles);
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= min_time)
            break;

        cycles *= 2;
    }
    /* approx Mflops */
    result = MonteCarlo_num_flops(cycles) / Stopwatch_read(Q) * 1.0e-6;
    Stopwatch_delete(Q);
    return result;
}

double kernel_measureSparseMatMult(int N, int nz, double min_time, Random R) {
    /* initialize vector multipliers and storage for result */
    /* y = A*y;  */
    double *x = RandomVector(N, R);
    double *y = (double *)malloc(sizeof(double) * N);
    double result = 0.0;

    // initialize square sparse matrix
    //
    // for this test, we create a sparse matrix with M/nz nonzeros
    // per row, with spaced-out evenly between the begining of the
    // row to the main diagonal.  Thus, the resulting pattern looks
    // like
    //             +-----------------+
    //             +*                +
    //             +***              +
    //             +* * *            +
    //             +** *  *          +
    //             +**  *   *        +
    //             +* *   *   *      +
    //             +*  *   *    *    +
    //             +*   *    *    *  +
    //             +-----------------+
    //
    // (as best reproducible with integer artihmetic)
    // Note that the first nr rows will have elements past
    // the diagonal.

    int nr = nz / N;            /* average number of nonzeros per row  */
    int anz = nr * N;           /* _actual_ number of nonzeros         */

    double *val = RandomVector(anz, R);
    int *col = (int *)malloc(sizeof(int) * nz);
    int *row = (int *)malloc(sizeof(int) * (N + 1));
    int r = 0;
    int cycles = 1;

    Stopwatch Q = new_Stopwatch();

    row[0] = 0;
    for (r = 0; r < N; r++) {
        /* initialize elements for row r */

        int rowr = row[r];
        int step = r / nr;
        int i = 0;

        row[r + 1] = rowr + nr;
        if (step < 1)
            step = 1;           /* take at least unit steps */

        for (i = 0; i < nr; i++)
            col[rowr + i] = i * step;
    }

    while (1) {
        Stopwatch_start(Q);
        SparseCompRow_matmult(N, y, val, row, col, x, cycles);
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= min_time)
            break;

        cycles *= 2;
    }

    /* approx Mflops */
    result = SparseCompRow_num_flops(N, nz, cycles) / Stopwatch_read(Q) * 1.0e-6;

    Stopwatch_delete(Q);
    free(row);
    free(col);
    free(val);
    free(y);
    free(x);
    return result;
}

double kernel_measureLU(int N, double min_time, Random R) {
    double **A = NULL;
    double **lu = NULL;
    int *pivot = NULL;

    Stopwatch Q = new_Stopwatch();
    double result = 0.0;
    int i = 0;
    int cycles = 1;

    if ((A = RandomMatrix(N, N, R)) == NULL)
        exit(1);
    if ((lu = new_Array2D_double(N, N)) == NULL)
        exit(1);
    if ((pivot = (int *)malloc(N * sizeof(int))) == NULL)
        exit(1);

    while (1) {
        Stopwatch_start(Q);
        for (i = 0; i < cycles; i++) {
            Array2D_double_copy(N, N, lu, A);
            LU_factor(N, N, lu, pivot);
        }
        Stopwatch_stop(Q);
        if (Stopwatch_read(Q) >= min_time)
            break;

        cycles *= 2;
    }
    /* approx Mflops */
    result = LU_num_flops(N) * cycles / Stopwatch_read(Q) * 1.0e-6;

    Stopwatch_delete(Q);
    free(pivot);
    Array2D_double_delete(N, N, lu);
    Array2D_double_delete(N, N, A);
    return result;
}


const double RESOLUTION_DEFAULT = 2.0;  /* secs (normally 2.0) */
const int RANDOM_SEED = 101010;

/* default: small (cache-contained) problem sizes */

const int FFT_SIZE = 1024;      /* must be a power of two */
const int SOR_SIZE = 100;       /* NxN grid */
const int SPARSE_SIZE_M = 1000;
const int SPARSE_SIZE_nz = 5000;
const int LU_SIZE = 100;

/* large (out-of-cache) problem sizes */

const int LG_FFT_SIZE = 1048576;        /* must be a power of two */
const int LG_SOR_SIZE = 1000;   /*  NxN grid  */
const int LG_SPARSE_SIZE_M = 100000;
const int LG_SPARSE_SIZE_nz = 1000000;
const int LG_LU_SIZE = 1000;

void print_banner() {
    printf ("**                                                              **\n");
    printf ("** SciMark2 Numeric Benchmark, see http://math.nist.gov/scimark **\n");
    printf ("** for details. (Results can be submitted to pozo@nist.gov)     **\n");
    printf ("**                                                              **\n");
}

int main(int argc, char *argv[]) {
    /* default to the (small) cache-contained version */
    double min_time = RESOLUTION_DEFAULT;
    int FFT_size = FFT_SIZE;
    int SOR_size = SOR_SIZE;
    int Sparse_size_M = SPARSE_SIZE_M;
    int Sparse_size_nz = SPARSE_SIZE_nz;
    int LU_size = LU_SIZE;

    /* run the benchmark */

    double res[6] = { 0.0 };
    Random R = new_Random_seed(RANDOM_SEED);

    if (argc > 1) {
        int current_arg = 1;

        if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0) {
            fprintf(stderr, "Usage: [-large] [minimum_time]\n");
            exit(0);
        }

        if (strcmp(argv[1], "-large") == 0) {
            FFT_size = LG_FFT_SIZE;
            SOR_size = LG_SOR_SIZE;
            Sparse_size_M = LG_SPARSE_SIZE_M;
            Sparse_size_nz = LG_SPARSE_SIZE_nz;
            LU_size = LG_LU_SIZE;
            current_arg++;
        }

        if (current_arg < argc) {
            min_time = atof(argv[current_arg]);
        }
    }

    print_banner();
    printf("Using %10.2f seconds min time per kenel.\n", min_time);

    res[1] = kernel_measureFFT(FFT_size, min_time, R);
    res[2] = kernel_measureSOR(SOR_size, min_time, R);
    res[3] = kernel_measureMonteCarlo(min_time, R);
    res[4] = kernel_measureSparseMatMult(Sparse_size_M, Sparse_size_nz, min_time, R);
    res[5] = kernel_measureLU(LU_size, min_time, R);
    res[0] = (res[1] + res[2] + res[3] + res[4] + res[5]) / 5;

    /* print out results  */
    printf("Composite Score:        %8.2f\n", res[0]);
    printf("FFT             Mflops: %8.2f    (N=%d)\n", res[1], FFT_size);
    printf("SOR             Mflops: %8.2f    (%d x %d)\n",
           res[2], SOR_size, SOR_size);
    printf("MonteCarlo:     Mflops: %8.2f\n", res[3]);
    printf("Sparse matmult  Mflops: %8.2f    (N=%d, nz=%d)\n", res[4],
           Sparse_size_M, Sparse_size_nz);
    printf("LU              Mflops: %8.2f    (M=%d, N=%d)\n", res[5],
           LU_size, LU_size);

    Random_delete(R);
    return 0;
}