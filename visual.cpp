#include <raylib.h>
#include <raymath.h>

#include <cfloat>
#include <unordered_map>
#include <vector>
#include <algorithm>

#define main main_
#include "main.cpp"
#undef main

const Color BACKGROUND_COLOR = BLACK;

static Vector2 screen;
static std::vector<size_t> dragging;

struct Cam {
    Vector2 pos = {};
    float _scale = 1;
    float target_scale = 1;

    void update();
    Vector2 to_screen(Vector2) const;
    Vector2 to_world(Vector2 p) const;
    void scale_add(float);
};

void Cam::update()
{
    float t = _scale/target_scale * GetFrameTime() * 30;
    _scale = Lerp(_scale, target_scale, t);
}

void Cam::scale_add(float ts)
{
    target_scale += ts;
    target_scale = Clamp(target_scale, 1, FLT_MAX);
}

Vector2 Cam::to_world(Vector2 p) const
{
    return (p - screen/2)/_scale - pos;
}

Vector2 Cam::to_screen(Vector2 p) const
{
    return (p + pos) * _scale + screen/2;
}

void draw_axis(Cam cam)
{
    auto c = cam.to_screen({0, 0});
    DrawLine(0, c.y, screen.x, c.y, WHITE);
    DrawLine(c.x, 0, c.x, screen.y, WHITE);
}

void draw_points(Cam cam, pointf_dynarr const &ps)
{
    for (size_t i = 0; i < ps.count; ++i) {
        auto &p = ps.items[i];
        Vector2 mpos = cam.to_world(GetMousePosition());

        Rectangle r = {p.x, p.y, 1, 1};
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) &&
            CheckCollisionPointRec(mpos, r)) {
            dragging.push_back(i);
        }

        (Vector2&)r = cam.to_screen((Vector2&)r);
        r.width = r.height = cam._scale;
        if (std::find_if(dragging.begin(), dragging.end(), [i](size_t &id) {
            return id == i;
        }) != dragging.end()) {
            if (dragging.size() > 1) {
                mpos = GetMouseDelta()/cam._scale;
                if (IsKeyDown(KEY_LEFT_SHIFT)) {
                    std::vector<pointf> dragged;

                    for (auto i : dragging) dragged.push_back(ps.items[i]);
                    pointf m = mean_point(dragged.data(), dragged.size());
                    pointf mp = m.sub((pointf&)mpos);
                    float cos = mp.cosine(m);
                    float ang = acosf(cos);
                    float sin = sinf(ang);
                    p = p.sub(m);
                    p = {p.x*cos + p.y*sin + m.x, p.x*-sin + p.y*cos + m.y};
                } else {
                    (Vector2&)p += mpos;
                }
            } else {
                if (IsKeyDown(KEY_LEFT_SHIFT))
                    mpos = {truncf(mpos.x) + 0.5f, truncf(mpos.y) + 0.5f};
                (Vector2&)p = mpos - Vector2{0.5, 0.5};
            }
        }
        DrawRectangleRec(r, WHITE);
    }
}

bool operator==(pointf a, pointf b)
{
    return a.x == b.x && a.y == b.y;
}

void draw_rects_and_inters(Cam cam, pointf_dynarr const &ps)
{
    rectf_dynarr rs = {};
    comp_ctx ctx;
    ctx.eps = 1e-2;
    ctx.tol = 1e-6;

    find_rects(ctx, ps, rs);
    char buf[256];
    sprintf(buf, "%zu", rs.count);
    DrawText(buf, 0, 0, 20, WHITE);
    for (size_t i = 0; i < rs.count; ++i) {
        auto r = rs.items[i];
        auto vs = r.vs;

        Vector2 mpos = cam.to_world(GetMousePosition()) - Vector2{0.5, 0.5};
        if (IsMouseButtonPressed(MOUSE_LEFT_BUTTON) &&
            r.has_point((pointf&)mpos)) {
            // I know this is suck
            for (size_t j = 0; j < ps.count; ++j) {
                for (size_t k = 0; k < RECT_VERT_COUNT; ++k) {
                    if (ps.items[j] == vs[k]) dragging.push_back(j);
                }
            }
        }

        pointf verts[] = {
            vs[3],
            vs[0],
            vs[1],
            vs[2],
            vs[3],
        };
        for (auto &v : verts)
            (Vector2&)v = cam.to_screen((Vector2&)v) + Vector2{0.5, 0.5}*cam._scale;
        Color c = RED;
        c.a = 0xaa;
        DrawTriangleStrip((Vector2*)verts, std::size(verts), c);
    }

    for (size_t i = 0; i < rs.count; ++i) {
        for (size_t j = i + 1; j < rs.count; ++j) {
            pointf ps[30];
            rectf &a = rs.items[i], &b = rs.items[j];
            int n = a.intersection_points(b, ps);
            for (int k = 0; k < n; ++k) {
                pointf p = ps[k].sub({-0.5, -0.5});
                Vector2 sc = cam.to_screen((Vector2&)p);
                Color c = BLUE;
                c.a = 0xaa;
                DrawCircleV(sc, 0.2*cam._scale, c);
            }
        }
    }

    delete[] rs.items;
}

int main()
{
    pointf_dynarr points{};
    Cam cam;
    cam._scale = 50;

    InitWindow(800, 600, "Visual");

    SetWindowState(FLAG_WINDOW_RESIZABLE);
    SetTargetFPS(60);

    while (!WindowShouldClose()) {
        BeginDrawing();
        screen = {(float)GetScreenWidth(), (float)GetScreenHeight()};

        cam.update();

        ClearBackground(BACKGROUND_COLOR);
        cam.scale_add(GetMouseWheelMove()*0.8);
        if (IsMouseButtonDown(MOUSE_BUTTON_RIGHT))
            cam.pos += GetMouseDelta()/cam._scale;

        if (IsMouseButtonPressed(MOUSE_BUTTON_MIDDLE)) {
            auto p = cam.to_world(GetMousePosition()) - Vector2{0.5, 0.5};
            points.add((pointf&)p);
        }

        if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT)) {
            dragging.clear();
        }

        draw_axis(cam);
        draw_points(cam, points);
        draw_rects_and_inters(cam, points);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
