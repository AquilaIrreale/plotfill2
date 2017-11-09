#include <cstdio>
#include <cstdlib>
#include <cstddef>

#include <algorithm>
#include <utility>
#include <vector>
#include <queue>

extern "C" {
#include <png.h>
}

#include "clipper.hpp"

using namespace std;
using namespace ClipperLib;

void read_png(FILE *fd, vector<vector<bool> > &mat, int &res);

void clean_image(vector<vector<bool> > &mat);
bool mat_test(vector<vector<bool> > &mat, size_t i, size_t j);
inline bool is_border(vector<vector<bool> > &mat, size_t i, size_t j);
void flood_reverse(vector<vector<bool> > &mat, size_t i, size_t j);
void follow_contour(vector<vector<bool> > &mat, size_t i, size_t j, Path &ret);
void identify_track(vector<vector<bool> > &mat, size_t i, size_t j, Paths &ret);
inline cInt px_to_pu(cInt px, int res);
void remove_duplicates(Path &ret); // TESTED
void douglas_peucker(Path &ret, double treshold); // TESTED
int poly_orientation(Path &poly);
void hpgl_print(Path &path);

int main(int argc, char *argv[])
{
    const int    pen_size = 20; // 0.5 mm
    const double treshold = 1.0;

    if (argc < 2) {
        fputs("plotfill: usage: plotfill filename\n", stderr);
        return 0;
    }

    FILE *fd = fopen(argv[1], "rb");

    if (!fd) {
        fputs("plotfill: ", stderr);
        perror(argv[1]);
        return 0;
    }

    int res;
    vector<vector<bool> > mat;
    read_png(fd, mat, res);
    if (!res) {
        fputs("DPI error\n", stderr);
        return 1;
    }

    fclose(fd);

    clean_image(mat);

    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            if (!is_border(mat, i, j)) {
                continue;
            }

            Paths track;
            identify_track(mat, i, j, track);

            for (auto &path : track) {
                for (auto &p : path) {
                    p.X = px_to_pu(p.X, res);
                    p.Y = px_to_pu(p.Y, res);
                }

                remove_duplicates(path);
                douglas_peucker(path, treshold);
            }

            Paths offseted;
            ClipperOffset co;
            co.AddPaths(track, jtSquare, etClosedPolygon);
            int off = -pen_size/2;
            while (co.Execute(offseted, off), !offseted.empty()) {
                for (auto &path : offseted) {
                    hpgl_print(path);
                }

                off -= pen_size;
            }
        }
    }

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
    }

    free(row_p);

    png_destroy_info_struct(png, &info);
}

void clean_image(vector<vector<bool> > &mat)
{
    for (size_t i = 0; i < mat.size(); i++) {
        for (size_t j = 0; j < mat[i].size(); j++) {
            if ((!mat_test(mat, i+1, j) && !mat_test(mat, i-1, j))
            ||  (!mat_test(mat, i, j+1) && !mat_test(mat, i, j-1))) {
                
                mat[i][j] = false;
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

void flood_reverse(vector<vector<bool> > &mat, size_t i, size_t j)
{
    if (mat.empty()
    ||  mat[0].empty()
    ||  i >= mat.size()
    ||  j >= mat[0].size()) {

        return;
    }

    bool fill = !mat[i][j];

    queue<pair<size_t, size_t> > q;
    q.emplace(i, j);
    while (!q.empty()) {
        pair<size_t, size_t> cur = q.front();
        q.pop();

        size_t y = cur.first;
        size_t w = cur.second;
        size_t e = cur.second;

        while (w > 0              && mat[y][w] != fill) w--;
        while (e < mat.size() - 1 && mat[y][e] != fill) e++;

        for (size_t x = w; x <= e; x++) {
            mat[y][x] = fill;

            if (y   > 0          && mat[y-1][x] != fill) q.emplace(y-1, x);
            if (y+1 < mat.size() && mat[y+1][x] != fill) q.emplace(y+1, x);
        }
    }
}

void follow_contour(vector<vector<bool> > &mat, size_t i, size_t j, Path &ret)
{
}

void identify_track(vector<vector<bool> > &mat, size_t i, size_t j, Paths &ret)
{
    ret.clear();

    Path ext_contour;
    follow_contour(mat, i, j, ext_contour);

    ret.push_back(ext_contour);

    for (IntPoint &p : ext_contour) {
        
        for (IntPoint cur(p.X + 1, p.Y);
             find(ext_contour.begin(),
                  ext_contour.end(),
                  cur) == ext_contour.end();
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

            flood_reverse(mat, cur.Y, cur.X + 1);
        }
    }

    flood_reverse(mat, i, j);
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
}

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
}


