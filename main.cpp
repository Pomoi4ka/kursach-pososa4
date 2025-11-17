#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cfloat>

#ifdef _OPENMP
#include <omp.h>
#endif

const float TOLERANCE = 1e-12;
const size_t PAIR_COUNT = 2;
const size_t RECT_VERT_COUNT = 4;

const bool TRACE = false;

typedef size_t pair_index[PAIR_COUNT];
typedef int (*compar_f)(void *, const void *, const void *);

static void swap(void *x, void *y, size_t n);
static void quicksort(void *xs, size_t x_size, size_t x_count,
               compar_f cmp, void *cmp_ctx);
static float normalized_atan2(float y, float x);
static bool mat_solve_homogeneous_sys(struct comp_ctx &ctx, float es[],
                                      float ans[], size_t rows, size_t cols,
                                      float tol);
static bool mat_solve_sys(struct comp_ctx &ctx, float es[], float ans[],
                          size_t rows, size_t cols, float tol);

struct comp_ctx {
    std::ofstream *prot;
    float eps;

    short level;
    short pad;

    void func(const char *);

    void ln();
    void write(const char *);
    void write(float);
    void write(int);
    void write(const float[], size_t);
    void write(const struct pointf[], size_t);
    void write(const float mat[], size_t, size_t);
    void write(struct pointf);
    void write(struct rectf);
};

struct pointf {
    float x, y;

    pointf sub(pointf a) const;
    float dot(pointf a) const;
    float len() const;
    float sqlen() const;
    float cosine(pointf a) const;
};

struct rectf {
    pointf vs[RECT_VERT_COUNT];
    comp_ctx *ctx;

    void clockwisify();
    int intersection_points(rectf const &a, pointf ps[]) const;
    bool has_point(pointf) const;
};

struct line_coeff {
    float k, b, c;
    comp_ctx *ctx;

    static line_coeff from_points(comp_ctx *, pointf, pointf);
    bool intersection(line_coeff, pointf&);
};

struct rectf_list {
    typedef rectf item;

    item *items;
    size_t count;
    size_t cap;

#ifdef _OPENMP
    omp_lock_t lock;
#endif

    void alloc(size_t);
    void add(rectf);
};

struct pointf_list {
    typedef pointf item;

    item *items;
    size_t count;
    size_t cap;

    void alloc(size_t);
    void add(pointf);
};

float pointf::cosine(pointf a) const
{
    return fabsf(dot(a))/sqrtf(sqlen()*a.sqlen());
}

float pointf::len() const
{
    return sqrtf(sqlen());
}

float pointf::sqlen() const
{
    return x*x + y*y;
}

float pointf::dot(pointf a) const
{
    return x * a.x + y * a.y;
}

pointf pointf::sub(pointf a) const
{
    pointf result = {x - a.x, y - a.y};
    return result;
}

static float normalized_atan2(float y, float x)
{
    const unsigned sign_mask = 0x80000000;
    const float b = 0.596227f;

    unsigned ux_s  = sign_mask & (unsigned &)x;
    unsigned uy_s  = sign_mask & (unsigned &)y;

    float q = (float)((~ux_s & uy_s) >> 29 | ux_s >> 30);

    float bxy_a = fabsf(b * x * y);
    float num = bxy_a + y * y;
    float atan_1q =  num / (x * x + bxy_a + num);

    unsigned uatan_2q = (ux_s ^ uy_s) | (unsigned &)atan_1q;
    return q + (float &)uatan_2q;
}

void comp_ctx::ln()
{
    if (!TRACE) return;
    write("\n");
    pad = level;
}

void comp_ctx::func(const char *name)
{
    if (!TRACE) return;
    pad = level = 0;
    write(name);
    write(":");
    level = 1;
    ln();
}

void comp_ctx::write(const float fs[], size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) write(", ");
        write(fs[i]);
    }
}

void comp_ctx::write(const pointf ps[], size_t n)
{
    for (size_t i = 0; i < n; ++i) {
        if (i > 0) write(", ");
        write(ps[i]);
    }
}

void comp_ctx::write(pointf p)
{
    write("(");
    write(p.x);
    write(", ");
    write(p.y);
    write(")");
}

void comp_ctx::write(rectf r)
{
    if (!TRACE) return;
    write("rectf {");
    level++; ln();
    write(r.vs, RECT_VERT_COUNT/2); write(","); ln();
    write(r.vs + RECT_VERT_COUNT/2, RECT_VERT_COUNT/2);
    level--; ln();
    write("}");
}

void comp_ctx::write(const float mat[], size_t rows, size_t cols)
{
    if (!TRACE) return;
    write("matrix [");
    level++; ln();
    for (size_t i = 0; i < rows; ++i) {
        if (i) ln();
        write("[");
        for (size_t j = 0; j < cols; ++j) {
            *prot << std::setw(13);
            write(mat[cols*i + j]);
        }
        write("]");
    }
    level--; ln();
    write("]");
}

void comp_ctx::write(float a)
{
    if (!TRACE) return;
    for (short i = 0; i < pad; ++i) *prot << "  ";
    *prot << a;
    pad = 0;
}

void comp_ctx::write(const char *text)
{
    if (!TRACE) return;
    for (short i = 0; i < pad; ++i) *prot << "  ";
    *prot << text;
    pad = 0;
}

void rectf_list::alloc(size_t n) {
    if (cap > count + n) return;

    while (cap <= count + n) {
        if (cap) cap *= 2;
        else     cap  = 1;
    }

    item *new_items = new item[cap];
    item *end       = new_items + count;
    for (item *a = new_items, *b = items; a != end; ++a, ++b)
        *a = *b;

    delete[] items;
    items = new_items;
}

void rectf_list::add(rectf p) {
#ifdef _OPENMP
    omp_set_lock(&lock);
#endif
    alloc(1);
    items[count++] = p;
#ifdef _OPENMP
    omp_unset_lock(&lock);
#endif
}

void pointf_list::alloc(size_t n) {
    if (cap > count + n) return;

    while (cap <= count + n) {
        if (cap) cap *= 2;
        else     cap  = 1;
    }

    item *new_items = new item[cap];
    item *end       = new_items + count;
    for (item *a = new_items, *b = items; a != end; ++a, ++b)
        *a = *b;

    delete[] items;
    items = new_items;
}

void pointf_list::add(pointf p) {
    alloc(1);
    items[count++] = p;
}

pointf mean_point(pointf ps[], size_t n)
{
    float x = 0;
    float y = 0;
    for (size_t i = 0; i < n; ++i) {
        x += ps[i].x;
        y += ps[i].y;
    }

    pointf mean = {x/n, y/n};
    return mean;
}

void rectf::clockwisify()
{
    pointf mean = mean_point(vs, RECT_VERT_COUNT);

    pointf vp[RECT_VERT_COUNT] = {
        vs[0].sub(mean),
        vs[1].sub(mean),
        vs[2].sub(mean),
        vs[3].sub(mean),
    };

    float as[RECT_VERT_COUNT] = {
        normalized_atan2(vp[0].y, vp[0].x),
        normalized_atan2(vp[1].y, vp[1].x),
        normalized_atan2(vp[2].y, vp[2].x),
        normalized_atan2(vp[3].y, vp[3].x),
    };
    ctx->func("rectf::clockwisify");
    ctx->write("средняя точка = ");
    ctx->write(mean);
    ctx->ln();
    ctx->write(*this);

    for (bool finished = false; !finished; ) {
        finished = true;
        for (size_t i = 1; i < RECT_VERT_COUNT; ++i) {
            if (as[i-1] >= as[i]) continue;
            finished = false;
            swap(&vs[i-1], &vs[i], sizeof *vs);
            swap(&as[i-1], &as[i], sizeof *as);
        }
    }
    ctx->write(" -> ");
    ctx->write(*this);
    ctx->ln();
}

bool rectf::has_point(pointf a) const
{
    const size_t rows = 2;
    const size_t cols = 3;
    const size_t mats = 2;

    float es[mats][rows*cols] = {
        {
            vs[1].x - vs[0].x, vs[2].x - vs[0].x, a.x - vs[0].x,
            vs[1].y - vs[0].y, vs[2].y - vs[0].y, a.y - vs[0].y,
        },
        {
            vs[2].x - vs[0].x, vs[3].x - vs[0].x, a.x - vs[0].x,
            vs[2].y - vs[0].y, vs[3].y - vs[0].y, a.y - vs[0].y,
        }
    };

    ctx->func("rectf::has_point");

    for (size_t i = 0; i < mats; ++i) {
        float xs[rows];
        float &alpha = xs[0], &beta = xs[1];

        bool exist = mat_solve_sys(*ctx, es[i], xs, rows, cols, TOLERANCE);
        if (!exist)                                continue;
        if (alpha < -ctx->eps || beta < -ctx->eps) continue;
        if (alpha + beta > 1 + ctx->eps)           continue;
        return true;
    }
    return false;
}

static size_t mat_pivot_row(const float es[], size_t rows, size_t cols,
                            size_t row, size_t col)
{
    size_t mi, i;
    float max = FLT_MIN;

    for (mi = i = row; i < rows; ++i) {
        float e = es[cols * i + col];
        if (fabsf(e) <= max) continue;
        max = e;
        mi = i;
    }
    return mi;
}

static size_t mat_reduce(float es[], size_t rows, size_t cols, float tol)
{
    size_t pivots = 0;

    // Я бы мог запихнуть es, rows, cols в структуру, но компилятор слишком
    // тупой чтобы увидеть, что у меня размер массивов весзде известен и
    // swap можно скопилить через пару инструкций, а не делать целый цикл
    // копирующий по одному float.
    for (size_t row = 0, col = 0; col < cols && row < rows; ++col) {
        size_t sel = mat_pivot_row(es, rows, cols, row, col);
        float pivot = es[cols * sel + col];
        if (fabsf(pivot) <= tol) continue;
        swap(&es[row * cols], &es[sel * cols], sizeof *es * cols);

        float invpivot = 1.f/pivot;
        for (size_t i = 0; i < cols; ++i) es[row * cols + i] *= invpivot;
        for (size_t i = 0; i < rows; ++i) {
            if (i == row) continue;
            float factor = es[i * cols + col];
            for (size_t j = 0; j < cols; ++j) {
                es[i * cols + j] -= factor * es[row * cols + j];
            }
        }
        ++row, ++pivots;
    }

    return pivots;
}

static bool mat_solve_sys(comp_ctx &ctx, float es[], float ans[], size_t rows,
                          size_t cols, float tol)
{
    ctx.func("matrixf::solve_sys");
    ctx.write(es, rows, cols); ctx.ln();
    mat_reduce(es, rows, cols, tol);
    ctx.write("Результат: "); ctx.write(es, rows, cols); ctx.ln();
    for (size_t i = 0; i < rows; ++i) {
        float b = es[i * cols + cols - 1];
        float sum = 0;
        for (size_t j = 0; j < cols - 1; sum += es[cols * i + j++])
            ;;
        if (fabsf(sum) <= tol && fabsf(b) > tol) {
            ctx.write("Есть противоречия => решений нет");
            ctx.ln();
            return false;
        }
        ans[i] = b;
    }
    ctx.write("Решение: "); ctx.write(ans, rows); ctx.ln();
    return true;
}

static bool mat_solve_homogeneous_sys(comp_ctx &ctx, float es[], float ans[],
                                      size_t rows, size_t cols, float tol)
{
    int pivots = mat_reduce(es, rows, cols, tol);
    int free_cols = cols - pivots;

    for (size_t i = 0; i < cols; ++i) ans[i] = 0.0;

    for (size_t i = cols - free_cols; i < cols; ++i) {
        ans[i] = 1.0;
        for (int j = 0; j < pivots; ++j)
            ans[j] += -es[cols * j + i];
    }

    return !free_cols;
}

static void swap(void *x, void *y, size_t n)
{
    char *a = reinterpret_cast<char *>(x);
    char *b = reinterpret_cast<char *>(y);
    char  t, *end = a + n;
    while (a != end) t = *a, *a = *b, *b = t, ++b, ++a;
}

static void quicksort_impl(void *items, size_t size, long low,
                    long high, compar_f cmp, void *cmp_ctx)
{
    char *xs    = reinterpret_cast<char *>(items);
    char *pivot = xs + high - size;
    long  lo    = low - size;
    long  hi    = high - size;

    if (low >= high) return;

    for (;;) {
        do lo += size; while (cmp(cmp_ctx, xs + lo, pivot) < 0);
        if (hi) do hi -= size; while (hi && cmp(cmp_ctx, xs + hi, pivot) > 0);
        if (lo >= hi) break;
        swap(xs + lo, xs + hi, size);
    }

    swap(xs + lo, pivot, size);
    quicksort_impl(xs, size, low, lo, cmp, cmp_ctx);
    quicksort_impl(xs, size, lo + size, high, cmp, cmp_ctx);
}

static void quicksort(void *xs, size_t x_size, size_t x_count,
                      compar_f cmp, void *cmp_ctx)
{
    quicksort_impl(xs, x_size, 0, x_count * x_size, cmp, cmp_ctx);
}

static int point_compare_clockwise(void *cx, const void *pp1,
                                   const void *pp2)
{
    pointf const &p1 = *reinterpret_cast<const pointf *>(pp1);
    pointf const &p2 = *reinterpret_cast<const pointf *>(pp2);

    pointf const &mean = *reinterpret_cast<pointf *>(cx);

    pointf a = p1.sub(mean);
    pointf b = p2.sub(mean);

    if (normalized_atan2(a.y, a.x) < normalized_atan2(b.y, b.x))
        return -1;
    return 1;
}

static void sort_points_clockwise(pointf ps[], size_t n)
{
    pointf mean = mean_point(ps, n);
    quicksort(ps, sizeof *ps, n, point_compare_clockwise, &mean);
}

bool line_coeff::intersection(line_coeff l, pointf &p)
{
    const size_t rows = 2;
    const size_t cols = 3;

    float es[] = {
        k  , b  , -c  ,
        l.k, l.b, -l.c
    };

    return mat_solve_sys(*ctx, es, (float *)&p, rows, cols, TOLERANCE);
}

line_coeff line_coeff::from_points(comp_ctx *ctx, pointf p1, pointf p2)
{
    line_coeff line;

    const size_t rows = 2;
    const size_t cols = 3;

    float es[] = {
        p1.x, p1.y, 1,
        p2.x, p2.y, 1,
    };

    mat_solve_homogeneous_sys(*ctx, es, (float *)&line, rows, cols, TOLERANCE);
    line.ctx = ctx;
    return line;
}

int rectf::intersection_points(rectf const &a, pointf ps[]) const
{
    int n = 0;
    line_coeff ls[2][RECT_VERT_COUNT];

    for (size_t i = 0; i < RECT_VERT_COUNT; ++i) {
        pointf p1, p2;
        p1 = vs[i];
        p2 = vs[(i + 1) % RECT_VERT_COUNT];
        ls[0][i] = line_coeff::from_points(ctx, p1, p2);

        p1 = a.vs[i];
        p2 = a.vs[(i + 1) % RECT_VERT_COUNT];
        ls[1][i] = line_coeff::from_points(ctx, p1, p2);
    }

    for (size_t i = 0; i < RECT_VERT_COUNT; ++i) {
        if (has_point(a.vs[i])) ps[n++] = a.vs[i];
        if (a.has_point(vs[i])) ps[n++] = vs[i];

        for (size_t j = 0; j < RECT_VERT_COUNT; ++j) {
            pointf p;
            if (!ls[0][i].intersection(ls[1][j], p)) continue;
            if (!has_point(p) || !a.has_point(p)) continue;
            ps[n++] = p;
        }
    }
    return n;
}

static float polygon_area(pointf ps[], size_t n)
{
    float sum = 0;
    for (size_t j = n - 1, i = 0; i < n; j = i++) {
        sum += (ps[i].x + ps[j].x) * (ps[i].y - ps[j].y);
    }
    return fabsf(sum) * 0.5;
}

static void add_rect(comp_ctx &c, pointf_list const &ps,
                     rectf_list &rs, size_t ids[RECT_VERT_COUNT])
{
    rectf r = {{}, &c};
    size_t n;
    for (n = 0; n < RECT_VERT_COUNT; n++) {
        r.vs[n] = ps.items[ids[n]];
    }

    r.clockwisify();

    pointf p12 = r.vs[0].sub(r.vs[1]);
    pointf p14 = r.vs[0].sub(r.vs[3]);
    pointf p23 = r.vs[1].sub(r.vs[2]);
    pointf p43 = r.vs[3].sub(r.vs[2]);

    const float eps = c.eps;

    if (fabsf(p12.len() - p23.len()) <= eps) return;

    if (p12.cosine(p14) >= eps) return;
    if (p23.cosine(p43) >= eps) return;
    if (p14.cosine(p43) >= eps) return;
    if (p12.cosine(p23) >= eps) return;
    rs.add(r);
}

static void find_rects(comp_ctx &ccx, pointf_list const &ps, rectf_list &rs)
{
#ifdef _OPENMP
    if (TRACE) {
        std::cerr << "Трассировка при параллельных "
            "вычслениях невозможна" << std::endl;
        return;
    }

#pragma omp target teams distribute parallel for
#endif
    for (size_t i = 0; i < ps.count; ++i) {
        for (size_t j = i + 1; j < ps.count; ++j) {
            for (size_t k = j + 1; k < ps.count; ++k) {
                for (size_t l = k + 1; l < ps.count; ++l) {
                    size_t ids[] = {i, j, k, l};
                    add_rect(ccx, ps, rs, ids);
                }
            }
        }
    }
}

static bool file_read_number(std::ifstream &f, float &x)
{
    f >> std::noskipws;
    while (!f.eof()) {
        if (f >> x) return true;
        f.clear(), f.get();
    }
    return false;
}

static void file_problem_err(const char *path)
{
    std::cerr << "ОШИБКА: Проблема с файлом " << path << std::endl;
}

static bool read_input_file(const char *path, comp_ctx &ctx, pointf_list &ps)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        file_problem_err(path);
        return false;
    }
    if (!file_read_number(f, ctx.eps)) {
        std::cerr << "ОШИБКА: неверный формат файла" << std::endl;
        return false;
    }

    while (!f.eof()) {
        float cs[2];

        for (size_t i = 0; i < sizeof cs / sizeof *cs; ++i) {
            if (file_read_number(f, cs[i])) continue;
            if (i == 0) break;
            std::cerr << "ОШИБКА: неверный формат файла"
                      << std::endl;
            return false;
        }
        pointf p = {cs[0], cs[1]};

        ctx.func("read_input_file");
        ctx.write("Прочитана точка ");
        ctx.write(p);
        ctx.ln();
        ps.add(p);
    }

    return true;
}

static bool find_max_intersect_area(comp_ctx &ctx, rectf_list const &rs,
                                    pair_index ids)
{
    const size_t MIN_RECT_COUNT = 2;
    if (rs.count < MIN_RECT_COUNT) {
        std::cerr << "ОШИБКА: "
            "невозможно найти пару прямоугольников: "
                  << rs.count << " элементов ("
                  << MIN_RECT_COUNT << " минимум)" << std::endl;
        return false;
    }

    float max_area = FLT_MIN;

    ids[0] = 0;
    ids[1] = 1;

#ifdef _OPENMP
    omp_lock_t lock;
    omp_init_lock(&lock);

#pragma omp target teams distribute parallel for    \
    map(from: max_area, ids) schedule(guided)
#endif
    for (size_t i = 0; i < rs.count; ++i) {
        for (size_t j = i + 1; j < rs.count; ++j) {
            pointf points[64];
            int n = rs.items[i].intersection_points
                (rs.items[j], points);
            if (n <= 2) continue;
            sort_points_clockwise(points, n);
            float area = polygon_area(points, n);
#ifdef _OPENMP
            omp_set_lock(&lock);
#endif
            if (area > max_area) {
                max_area = area;
                ids[0] = i;
                ids[1] = j;
            }
#ifdef _OPENMP
            omp_unset_lock(&lock);
#endif
        }
    }
    return true;
}

static void print_verts(rectf_list const &rects, pair_index max)
{
    std::cout << "Вершины найденной пары прямоугольников с "
              << "максимальной площадью пересечения:" << std::endl;
    for (size_t i = 0; i < PAIR_COUNT; ++i) {
        for (size_t j = 0; j < RECT_VERT_COUNT; ++j) {
            if (j > 0) std::cout << " ";
            std::cout << "(" << rects.items[max[i]].vs[j].x
                      << ", " << rects.items[max[i]].vs[j].y
                      << ")";
        }
        std::cout << std::endl;
    }
}

int main()
{
    const char *input = "input.txt", *prot_path = "protocol.txt";

    std::ofstream prot;
    comp_ctx ctx = {};
    pointf_list points = {};
    rectf_list rects = {};
    pair_index max;

#ifdef _OPENMP
    omp_init_lock(&rects.lock);
#endif

    if (TRACE) {
        prot.open(prot_path);
        if (!prot) {
            file_problem_err(prot_path);
            return 1;
        }
        ctx.prot = &prot;
    }

    if (!read_input_file(input, ctx, points)) return 1;
    find_rects(ctx, points, rects);

    if (!find_max_intersect_area(ctx, rects, max)) return 1;
    print_verts(rects, max);

    delete[] rects.items;
    delete[] points.items;
    return 0;
}
