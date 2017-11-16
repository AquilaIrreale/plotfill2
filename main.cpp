#include <cstdio>
#include <cstdlib>
#include <cstddef>

#include <algorithm>
#include <utility>
#include <vector>
#include <queue>
#include <unordered_set>

extern "C" {
#include <png.h>
}

#include "clipper.hpp"

#define DP_TRESHOLD     0.3 // PU
#define INLINE_TRESHOLD 0.8 // PU
#define PEN_SIZE        20  // PU (0.5 mm)

#define BIGPRIME 32416190071ULL

#define ACTION(s) fputs(s "... ", stderr)
#define DONE()    fputs("DONE\n", stderr)

namespace std {
    template <> struct hash<ClipperLib::IntPoint> {
        size_t operator ()(const ClipperLib::IntPoint &p) const {
            std::hash<ClipperLib::cInt> hasher;
            return (BIGPRIME + hasher(p.X)) * BIGPRIME + hasher(p.Y);
        }
    };
}

using namespace std;
using namespace ClipperLib;

void read_png(FILE *fd, vector<vector<bool> > &mat, int &res);

void clean_image(vector<vector<bool> > &mat);
bool mat_test(vector<vector<bool> > &mat, size_t i, size_t j);
inline bool is_border(vector<vector<bool> > &mat, size_t i, size_t j);
void flood_fill(vector<vector<bool> > &mat, size_t i, size_t j, bool fill);
void follow_contour(vector<vector<bool> > &mat, size_t i, size_t j, Path &ret);
void identify_track(vector<vector<bool> > &mat, size_t i, size_t j, Paths &ret);
inline cInt px_to_pu(cInt px, int res);
void remove_duplicates(Path &ret); // TESTED
void douglas_peucker(Path &ret, double treshold); // TESTED
void line_simplify(Path &path);
int poly_orientation(Path &poly);
void hpgl_print(Path &path);

int main(int argc, char *argv[])
{
    const int    pen_size = PEN_SIZE;
    const double treshold = DP_TRESHOLD;
    bool fill = false;

    if (argc < 2) {
        fputs("plotfill2: usage: plotfill filename\n", stderr);
        return 0;
    }


    FILE *fd = fopen(argv[1], "rb");

    if (!fd) {
        fputs("plotfill2: ", stderr);
        perror(argv[1]);
        return 0;
    }

    ACTION("=> Reading file");

    int res;
    vector<vector<bool> > mat;
    read_png(fd, mat, res);
    if (!res) {
        fputs("DPI error\n", stderr);
        return 1;
    }

    fclose(fd);

    DONE();
    ACTION("=> Image pre-processing");

    clean_image(mat);

    DONE();

    // DEBUG
    /*
    for (auto v = mat.rbegin(); v != mat.rend(); ++v) {
        for (auto p = v->begin(); p != v->end(); ++p) {
            putchar(*p ? '#' : '.');
            putchar(' ');
        }

        putchar('\n');
    }

    unsigned x, y;
    scanf("%u%u", &x, &y);

    Path test;
    follow_contour(mat, y, x, test);

    for (auto &p : test) {
        printf("(%llu, %llu)\n", p.X, p.Y);
    }

    return 0;
    */
    // DEBUG END


    ACTION("=> Processing tracks");
    putc('\n', stderr);

    unsigned t = 0;

    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            if (!is_border(mat, i, j)) {
                continue;
            }

            fprintf(stderr, "=> Found track %u at (%lu, %lu)\n", t++, j, i);

            Paths track;
            identify_track(mat, i, j, track);

            for (auto &path : track) {
                for (auto &p : path) {
                    p.X = px_to_pu(p.X, res);
                    p.Y = px_to_pu(p.Y, res);
                }

                /*
                transform(
                    path.begin(),
                    path.end(),
                    path.begin(),

                    [&res](IntPoint &p) {
                        return IntPoint(px_to_pu(p.X, res),
                                        px_to_pu(p.Y, res));
                    }
                );
                */

                remove_duplicates(path);
                line_simplify(path);
                douglas_peucker(path, treshold);
            }

            // DEBUG
            /*
            Paths offp;
            ClipperOffset z;
            z.AddPaths(track, jtRound, etClosedPolygon);
            z.Execute(offp, -pen_size/2);
            if (offp.empty()) {
                puts(":(");
                return 0;
            }
            for (auto &p : offp) {
                hpgl_print(p);
            }

            return 0;
            */
            //

            ACTION("=> Offsetting");
            putc('\n', stderr);

            Paths offseted;
            ClipperOffset co;
            co.AddPaths(track, jtRound, etClosedPolygon);
            int off = -pen_size/2;
            int l = 0;
            while (co.Execute(offseted, off), !offseted.empty()) {
                fprintf(stderr, "   -> Layer %d\n", l++);
                for (auto &path : offseted) {
                    if (path.front() != path.back()) {
                        path.push_back(path.front());
                    }

                    hpgl_print(path);
                }

                off -= pen_size;

                if (!fill) {
                    break;
                }
            }

            if (l == 0) {
                fputs("   -> Cannot offset (too narrow?)\n", stderr);
            }
        }
    }

    puts("PU\n0,0;\n");

    return 0;
}

void read_png(FILE *fd, vector<vector<bool> > &mat, int &res)
{
    png_byte color_type;
    png_byte bit_depth;
    png_bytepp row_p;

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING,
                                             NULL, NULL, NULL);
    png_infop info = png_create_info_struct(png);

    png_init_io(png, fd);
    png_read_info(png, info);

    int w      = png_get_image_width(png, info);
    int h      = png_get_image_height(png, info);
    color_type = png_get_color_type(png, info);
    bit_depth  = png_get_bit_depth(png, info);

    if (bit_depth > 8) {
        png_set_strip_16(png);
    }

    if (color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(png);
    }

    if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(png);
    }

    if(color_type == PNG_COLOR_TYPE_GRAY ||
        color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {

        png_set_gray_to_rgb(png);
    }

    if (color_type & PNG_COLOR_MASK_ALPHA) {
        png_set_strip_alpha(png);
    }

    png_read_update_info(png, info);

    res = png_get_pixels_per_meter(png, info);

    row_p = (png_bytepp)malloc(sizeof(png_bytep) * h);
    for (int i = 0; i < h; i++) {
        row_p[i] = (png_bytep)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_p);

    mat.clear();
    mat.resize(h, vector<bool>(w, false));
    for (int i = 0; i < h; i++) {
        for (int j = 0; j < w; j++) {
            for (int c = 0; c < 3; c++) {
                if (row_p[i][j*3 + c] < 255) {
                    mat[h-1-i][j] = true;
                    break;
                }
            }
        }

        free(row_p[i]);

        mat[i].shrink_to_fit();
    }

    free(row_p);

    mat.shrink_to_fit();

    png_destroy_info_struct(png, &info);
}

void clean_image(vector<vector<bool> > &mat)
{
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            if (mat[i][j]) {
                if ((!mat_test(mat, i+1, j) && !mat_test(mat, i-1, j))
                ||  (!mat_test(mat, i, j+1) && !mat_test(mat, i, j-1))) {
                    
                    mat[i][j] = false;
                }
            } else {
                if ((mat_test(mat, i+1, j) && mat_test(mat, i-1, j))
                ||  (mat_test(mat, i, j+1) && mat_test(mat, i, j-1))) {

                    mat[i][j] = true;
                }
            }
        }
    }
}

bool mat_test(vector<vector<bool> > &mat, size_t i, size_t j)
{
    if (i < 0 || j < 0
    ||  i >= mat.size()
    ||  j >= mat[0].size()) {
        return false;
    }

    return mat[i][j];
}

bool is_border(vector<vector<bool> > &mat, size_t i, size_t j)
{
    if (!mat_test(mat, i, j)) {
        return false;
    }

    int neighb = 0;
    if (mat_test(mat, i-1, j)) neighb++;
    if (mat_test(mat, i+1, j)) neighb++;
    if (mat_test(mat, i, j-1)) neighb++;
    if (mat_test(mat, i, j+1)) neighb++;

    switch (neighb) {
    case 3:
        return true;

    case 4:
    case 1:
    case 0:
        return false;

    case 2:
    default:
        break;
    }

    if ((mat_test(mat, i+1, j) && mat_test(mat, i-1, j))
    ||  (mat_test(mat, i, j+1) && mat_test(mat, i, j-1))) {
        return false;
    }

    return true;
}

void flood_fill(vector<vector<bool> > &mat, size_t i, size_t j, bool fill)
{
    if (mat.empty()
    ||  mat[0].empty()
    ||  i >= mat.size()
    ||  j >= mat[0].size()) {

        return;
    }

    queue<pair<size_t, size_t> > q;
    q.emplace(i, j);
    while (!q.empty()) {
        size_t y = q.front().first;
        size_t w = q.front().second;
        size_t e = q.front().second;

        q.pop();

        if (mat[y][w] == fill) {
            continue;
        }

        while (w > 0              && mat[y][w-1] != fill) w--;
        while (e < mat[0].size() - 1 && mat[y][e+1] != fill) e++;

        for (size_t x = w; x <= e; x++) {
            mat[y][x] = fill;

            if (y   > 0          && mat[y-1][x] != fill) q.emplace(y-1, x);
            if (y+1 < mat.size() && mat[y+1][x] != fill) q.emplace(y+1, x);
        }
    }
}

// DELETE?
/*
void gen_neighbors(vector<vector<bool> > &mat, size_t i, size_t j,
                   vector<pair<size_t, size_t> > &ret)
{
    ret.clear();

    bool top    = i >= mat.size() - 1;
    bool right  = j >= mat[0].size() - 1;
    bool bottom = i == 0;
    bool left   = j == 0;

    if (!top) {
        ret.emplace(i+1, j);
    }

    if (!right) {
        if (!top) {
            ret.emplace(i+1, j+1);
        }

        ret.emplace(i, j+1);

        if (!bottom) {
            ret.emplace(i-1, j+1);
        }
    }

    if (!bottom) {
        ret.emplace(i-1, j);
    }

    if (!left) {
        if (!bottom) {
            ret.emplace(i-1, j-1);
        }

        ret.emplace(i, j-1);

        if (!top) {
            ret.emplace(i+1, j);
        }
    }
}
*/

bool is_contour(vector<vector<bool> > & mat, IntPoint &to, IntPoint &from)
{
    if (to.X < 0 || to.Y < 0 ||
        (size_t) to.X > mat[0].size() || (size_t) to.Y > mat.size() ||
        from.X < 0 || from.Y < 0 ||
        (size_t) from.X > mat[0].size() || (size_t) from.Y > mat.size()) {

        return false;
    }

    cInt dx = to.X - from.X;
    cInt dy = to.Y - from.Y;

    cInt ax, bx, ay, by;
    if (llabs(dx) == 1 && dy == 0) {
        ax = min(to.X, from.X);
        bx = ax;
        ay = to.Y;
        by = to.Y - 1;
    } else if (llabs(dy) == 1 && dx == 0) {
        ax = to.X;
        bx = to.X - 1;
        ay = min(to.Y, from.Y);
        by = ay;
    } else {
        return false;
    }

    /*
    bool a, b;
    if (ax >= 0 && (size_t)ax < mat.size() &&
        ay >= 0 && (size_t)ay < mat[0].size()) {
        
        a = mat[ay][ax];
    } else {
        a = false;
    }

    if (bx >= 0 && (size_t)bx < mat.size() &&
        by >= 0 && (size_t)by < mat[0].size()) {
        
        b = mat[by][bx];
    } else {
        b = false;
    }
    */

    return mat_test(mat, ay, ax) != mat_test(mat, by, bx);
}

bool on_edge(vector<vector<bool> > &mat, IntPoint &p)
{
    array<IntPoint, 4> px;
    px[0].X = p.X;
    px[0].Y = p.Y;
    px[1].X = p.X;
    px[1].Y = p.Y - 1;
    px[2].X = p.X - 1;
    px[2].Y = p.Y - 1;
    px[3].X = p.X - 1;
    px[3].Y = p.Y;

    size_t count = 0;
    for (auto &it : px) {
        if (it.Y >= 0 && (size_t) it.Y < mat.size() &&
            it.X >= 0 && (size_t) it.X < mat[0].size() &&
            mat[it.Y][it.X]) {

            count++;
        }
    }

    return count > 0 && count < 4;
}

void gen_neighbors(IntPoint p, array<IntPoint, 4> &ret)
{
    ret[0].X = p.X + 0;
    ret[0].Y = p.Y + 1;
    ret[1].X = p.X + 1;
    ret[1].Y = p.Y + 0;
    ret[2].X = p.X - 0;
    ret[2].Y = p.Y - 1;
    ret[3].X = p.X - 1;
    ret[3].Y = p.Y - 0;
}

void gen_corners(IntPoint p, array<IntPoint, 4> &ret)
{
    ret[0].X = p.X;
    ret[0].Y = p.Y;
    ret[1].X = p.X;
    ret[1].Y = p.Y + 1;
    ret[2].X = p.X + 1;
    ret[2].Y = p.Y + 1;
    ret[3].X = p.X + 1;
    ret[3].Y = p.Y;
}

void follow_contour(vector<vector<bool> > &mat,
                    size_t i, size_t j, Path &ret)
{
    ret.clear();

    IntPoint base(-1, -1);

    array<IntPoint, 4> corners;
    gen_corners(IntPoint(j, i), corners);

    for (IntPoint &p : corners) {
        if (on_edge(mat, p)) {
            base = p;
            break;
        }
    }

    if (base.X < 0 || base.Y < 0) {
        return;
    }

    ret.push_back(base);

    array<IntPoint, 4> neighbors;
    gen_neighbors(base, neighbors);

    array<IntPoint, 4>::iterator p;
    for (p = neighbors.begin(); p != neighbors.end(); ++p) {
        if (is_contour(mat, *p, base)) {
            break;
        }
    }

    if (p == neighbors.end()) {
        ret.clear();
        return;
    }

    ret.push_back(*p);
    
    while (*p != base) {
        IntPoint cur = *p;
        gen_neighbors(cur, neighbors);
        
        for (p = neighbors.begin(); p != neighbors.end(); ++p) {
            if (*p != ret.rbegin()[1] && is_contour(mat, *p, cur)) {
                break;
            }
        }

        if (p == neighbors.end()) {
            ret.clear();
            return;
        }

        ret.push_back(*p);
    }
}

/*
void follow_contour(vector<vector<bool> > &mat,
                    size_t i, size_t j, Path &ret)
{
    vector<pair<size_t, size_t> > neighbors;
    gen_neighbors(mat, i, j, neighbors);
    auto p = neighbors.begin();
    while (mat[p->first][p->second]) {
        ++p;
    }

    if (p == neighbors.begin()) {
        p = neighbors.end();
    } else {
        --p;
    }

    path.clear();
    path.push_back(IntPoint(j, i));

    do {
    }
}
*/

bool touches_contour(Path &contour, IntPoint p)
{
    array<IntPoint, 4> corners;
    gen_corners(p, corners);

    for (auto &c : corners) {
        if (find(contour.begin(),
                 contour.end(), c) != contour.end()) {

            return true;
        }
    }

    return false;
}

void contour_to_pixels(Path &contour,
                       unordered_set<IntPoint> &contour_pixels)
{
    contour_pixels.clear();

    if (contour.front() != contour.back()) {
        contour.push_back(contour.front());
    }

    for (Path::iterator a = contour.begin(),
                        b = contour.begin() + 1;
         b != contour.end();
         ++a, ++b) {

        cInt dx = b->X - a->X;
        cInt dy = b->Y - a->Y;

        if (dx > 0 && dy > 0) {
            exit(1);
        }

        if (dx == 0 && dy == 0) {
            continue;
        }

        if (dx == 1) {
            contour_pixels.insert(IntPoint(a->X, a->Y - 1));
        } else if (dx == -1) {
            contour_pixels.insert(IntPoint(b->X, b->Y));
        } else if (dy == 1) {
            contour_pixels.insert(IntPoint(a->X, a->Y));
        } else if (dy == -1) {
            contour_pixels.insert(IntPoint(b->X - 1, b->Y));
        } else {
            exit(2);
        }
    }
}

void identify_track(vector<vector<bool> > &mat,
                    size_t i, size_t j, Paths &ret)
{
    ret.clear();

    Path ext_contour;
    follow_contour(mat, i, j, ext_contour);

    if (poly_orientation(ext_contour) < 0) {
        reverse(ext_contour.begin(), ext_contour.end());
    }

    ret.push_back(ext_contour);

    unordered_set<IntPoint> contour_pixels;
    contour_to_pixels(ext_contour, contour_pixels);

    queue<IntPoint> q;
    for (const IntPoint &p : contour_pixels) {
        q.push(p);
    }

    while (!q.empty()) {
        IntPoint base = q.front();
        q.pop();

        if (!mat_test(mat, base.Y, base.X + 1)) {
            // base is on a left side
            continue;
        }

        for (IntPoint p(base.X + 1, base.Y);
             contour_pixels.count(p) == 0;
             p.X++) {
            
            if (!is_border(mat, p.Y, p.X)) {
                continue;
            }

            fprintf(stderr, "   -> Found hole at (%lld, %lld)\n",
                    p.X, p.Y);

            /*
            fprintf(stderr, "%lu, N:%lu, E:%lu, S:%lu, W:%lu\n",
                    contour_pixels.count(IntPoint(p.X, p.Y)),
                    contour_pixels.count(IntPoint(p.X, p.Y+1)),
                    contour_pixels.count(IntPoint(p.X+1, p.Y)),
                    contour_pixels.count(IntPoint(p.X, p.Y-1)),
                    contour_pixels.count(IntPoint(p.X-1, p.Y)));
            */

            Path int_contour;
            follow_contour(mat, p.Y, p.X, int_contour);

            if (poly_orientation(int_contour) > 0) {
                reverse(int_contour.begin(), int_contour.end());
            }

            unordered_set<IntPoint> int_pixels;
            contour_to_pixels(int_contour, int_pixels);

            for (const IntPoint &p : int_pixels) {
                q.push(p);
            }

            contour_pixels.insert(int_pixels.begin(),
                                  int_pixels.end());

            ret.push_back(int_contour);

            break;
        }
    }

    /*
    for (IntPoint &p : ext_contour_pixels) {
        for (IntPoint cur(p.X + 1, p.Y);
             !touches_contour(ext_contour, cur);
             cur.X++) {

            if (!mat[cur.Y][cur.X + 1]) {
                continue;
            }

            Path int_contour;
            follow_contour(mat, cur.Y, cur.X, int_contour);

            // SAFETY
            if (poly_orientation(int_contour) > 0) {
                reverse(int_contour.begin(), int_contour.end());
            }
            // SAFETY

            ret.push_back(int_contour);

            flood_fill(mat, cur.Y, cur.X + 1);
        }
    }
    */

    flood_fill(mat, i, j, false);
}

inline cInt px_to_pu(cInt px, int res)
{
    return px * 40000 / res;
}

void remove_duplicates(Path &ret)
{
    if (ret.size() < 2) {
        return;
    }

    IntPoint cur(-1, -1);
    Path::iterator j = ret.begin(); 
    for (auto i = ret.begin(); i != ret.end(); ++i) {
        if (*i != cur) {
            *j++ = *i;
            cur = *i;
        }
    }

    ret.erase(j, ret.end());
    ret.shrink_to_fit();
}

double dist(IntPoint &p1, IntPoint &p2)
{
    double x = p2.X - p1.X;
    double y = p2.Y - p1.Y;

    return fabs(x*x + y*y);
}

double dist(IntPoint &p1, IntPoint &p2, IntPoint &q)
{
    if (p1 == p2) {
        return dist(p1, q);
    }

    double xb = p2.X - p1.X;
    double yb = p2.Y - p1.Y;
    double xc =  q.X - p1.X;
    double yc =  q.Y - p1.Y;

    return fabs(xb*yc - xc*yb) / sqrt(xb*xb + yb*yb);
}

void douglas_peucker_r(Path::iterator l,
                       Path::iterator r,
                       double treshold)
{
    if (l == r || l == r-1) {
        return;
    }

    Path::iterator pivot;
    double max_dist = -1;

    for (Path::iterator i = l+1; i != r; i++) {
        double cur_dist = dist(*l, *r, *i);
        if (cur_dist  > max_dist) {
            max_dist = cur_dist;
            pivot = i;
        }
    }

    if (max_dist >= treshold) {
        douglas_peucker_r(l, pivot, treshold);
        douglas_peucker_r(pivot, r, treshold);
    } else {
        fill(l+1, r, IntPoint(-1, -1));
    }
}

void douglas_peucker(Path &ret, double treshold)
{
    if (ret.size() < 2) {
        return;
    }

    if (ret.front() != ret.back()) {
        ret.push_back(ret.front());
    }

    douglas_peucker_r(ret.begin(), ret.end()-1, treshold);

    Path::iterator new_end = remove_if(ret.begin(),
                                       ret.end(),
                                       [](IntPoint &p)
                                       {return p.X == -1;});
    
    ret.erase(new_end, ret.end());
    ret.shrink_to_fit();
}

bool is_inline(Path::iterator a, Path::iterator b)
{
    if (a == b || a + 1 == b) {
        return true;
    }

    for (Path::iterator i = a + 1; i != b; ++i) {
        if (dist(*a, *b, *i) > INLINE_TRESHOLD) {
            return false;
        }
    }

    return true;
}

void line_simplify(Path &path)
{
    if (path.size() < 3) {
        return;
    }

    Path::iterator a = path.begin();
    Path::iterator b = a;

    while (a + 1 != path.end()) {
        while (b != path.end() && is_inline(a, b)) {
            ++b;
        }

        --b;

        if (a != b && a + 1 != b) {
            fill(a + 1, b, IntPoint(-1, -1));
        }

        a = b;
    }

    Path::iterator new_end = remove_if(path.begin(),
                                       path.end(),
                                       [](IntPoint &p)
                                       {return p.X == -1;});
    
    path.erase(new_end, path.end());

    Path junction(3);
    junction[0] = path.rbegin()[1];
    junction[1] = path.front();
    junction[2] = path[1];

    if (is_inline(junction.begin(), junction.end())) {
        path.pop_back();
        path.front() = path.back();
    }

    path.shrink_to_fit();
}

// Positive = clockwise
int poly_orientation(Path &poly)
{
    if (poly.size() < 2) {
        return 0;
    }

    cInt ret(0);

    auto i = poly.begin();
    auto j = poly.begin() + 1;
    for (; j != poly.end(); ++i, ++j) {
        ret += (j->X - i->X) * (j->Y + i->Y);
    }

    ret += (i->X - poly[0].X) * (i->Y + poly[0].Y);

    if (ret > 0) {
        return 1;
    } else if (ret < 0) {
        return -1;
    } else {
        return 0;
    }
}

void hpgl_print(Path &path)
{
    printf("PU\n%lld,%lld;\nPD", path.front().X, path.front().Y);

    Path::iterator p;
    for (p = path.begin() + 1; p != path.end(); ++p) {
        printf("\n%lld,%lld", p->X, p->Y);

        if (p + 1 != path.end()) {
            putchar(',');
        } else {
            putchar(';');
        }
    }

    putchar('\n');
}


