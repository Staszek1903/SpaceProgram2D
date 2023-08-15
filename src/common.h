#ifndef COMMON_H
#define COMMON_H

#include "raylib.h"
#include "rlgl.h"
#include "gui.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include <functional>
#include <cmath> 
#include <memory>
#include <array>
#include <assert.h>
#include <limits>
#include <algorithm>

typedef std::vector<Vector2> ArrayV2;

//--------------------------------------
// CONSTS
//--------------------------------------
constexpr int screenWidth = 800;
constexpr int screenHeight = 450;
constexpr float EPSILON = 0.0001;

//--------------------------------------
// STRUCTS
//--------------------------------------
enum PartID{
    POD,
    FUELTANK,
    LARGEENGINE,
    SMALLENGINE,
    SEPARATOR,
    COUNT
};

enum PartType{
    NONE,
    ENGINE,
    TANK
};

struct PartItem
{
    PartID id;
    Vector2 position = {0.0,0.0};
    float scale = 1.0;
    Color color = GRAY;//{GRAY.r,GRAY.g,GRAY.b, 128};
    int links_joined = {0};

    PartType type = {PartType::NONE};
    int info_index = {-1};
};

bool operator==(const PartItem & a, const PartItem & b){
    return (&a == &b);
}

struct EngineInfo{
    float max_thrust;
    float isp;
    float exhaust_width;
    float exhaust_offset;
};

struct EnginePart{
    int part_item_index = {-1};
    float throttle = {0};
    bool active = {false};
};

struct FueltankInfo{
    float capacity = {0};
};

struct FueltankPart{
    int part_item_index = {-1};
    float fuel_amount = {1.0};
};

struct Vessel{
    std::vector<PartItem> parts{};
    std::vector<EnginePart> engines{};
    std::vector<FueltankPart> tanks{};
    Vector2 position{0.0f,0.0f};
    Vector2 velocity{0.0f,0.0f}; 
    float rotation_deg{0.0f};
    float angular_vel_deg{0.0f};
};

//--------------------------------------
// CONTAINERS
//--------------------------------------

std::vector<std::string> part_names = {
    "POD",
    "FUELTANK",
    "L ENGINE",
    "S ENGINE",
    "SEPARATOR"
};

std::vector<PartType> part_type = {
    PartType::NONE,
    PartType::TANK,
    PartType::ENGINE,
    PartType::ENGINE,
    PartType::NONE
};

std::vector<int> part_info_index = {
    -1,
    0,
    0,
    1,
    -1
};

std::vector<ArrayV2> parts_poly ={
    {{-0.2,-0.2},{-0.6,0.8},{0.6,0.8},{0.2,-0.2}},
    {{-0.6,-1},{-0.6,1},{0.6,1},{0.6,-1}},
    {{-0.2,-0.6},{-0.4,0.6},{0.4,0.6},{0.2,-0.6}},
    {{-0.2,-0.2},{-0.4,0.2},{0.4,0.2},{0.2,-0.2}},
    {{-0.6,-0.2},{-0.6,0.0},{0.6,0.0},{0.6,-0.2}}
};

std::vector<ArrayV2> linkages ={
    {{0.0,0.8}},
    {{0.0,1.0}, {0.0,-1.0}},
    {{0.0,-0.6},{0.0,0.6}},
    {{0.0,-0.2},{0.0,0.2}},
    {{0.0,0.0}, {0.0,-0.2}},
};

std::vector<EngineInfo> engine_info = {
    {2.0, 1.0, 0.8, 0.6},
    {1.0, 1.0, 0.8, 0.2},
};

std::vector<FueltankInfo> fuel_tank_info = {
    {100.0}
};

std::vector<Button> buttons;
std::vector<TextInput> text_inputs;


//--------------------------------------
// FUNCTIONS
//--------------------------------------

void update_gui(){
    Vector2 mouse_pos = GetMousePosition();
    bool pressed = IsMouseButtonPressed(0);
    int key = GetKeyPressed();

    for(Button & button: buttons){
        button.update(mouse_pos,pressed);
    }

    for(TextInput & input: text_inputs){
        input.update(mouse_pos,pressed,key);
    }
}

void draw_gui(){
    for(Button & button: buttons){
        button.draw();
    }
    
    for(TextInput & input: text_inputs){
        input.draw();
    }
}

void clear_gui_data(){
    buttons.clear();
    text_inputs.clear();
}

auto & get_part(PartID part){
    return parts_poly.at(part);
}

Vector2 operator-(Vector2 a, Vector2 b){
    return {a.x-b.x, a.y-b.y};
}

Vector2 operator-(Vector2 v){
    return {-v.x, -v.y};
}

Vector2 operator+(Vector2 a, Vector2 b){
    return {a.x+b.x, a.y+b.y};
}

Vector2 operator-=(Vector2 &a, Vector2 b){
    a = {a.x-b.x, a.y-b.y};
    return a;
}

Vector2 operator+=(Vector2 &a, Vector2 b){
    a = {a.x+b.x, a.y+b.y};
    return a;
}

Vector2 operator/(Vector2 a, float s){
    return {a.x/s, a.y/s};
}

Vector2 operator*(Vector2 a, float s){
    return {a.x*s, a.y*s};
}

Vector2 operator/=(Vector2 &a, float s){
    a = {a.x/s, a.y/s};
    return a;
}

Vector2 operator*=(Vector2 &a, float s){
    a = {a.x*s, a.y*s};
    return a;
}

bool IsFloatEqual(float a, float b){
    return( a-b < EPSILON && b-a < EPSILON );
}

bool operator==(Vector2 a, Vector2 b){
    return( IsFloatEqual(a.x,b.x) && IsFloatEqual(a.y,b.y) );
}

bool operator!=(Vector2 a, Vector2 b){
    return !(a==b);
}

float cross(Vector2 a, Vector2 b){
    return {a.x * b.y - a.y * b.x};
}

float dot(Vector2 a, Vector2 b){
    return {a.x * b.x + a.y * b.y};
}

float lensq(Vector2 v){
    return v.x*v.x + v.y*v.y;
}

float len(Vector2 v){
    return std::sqrt(lensq(v));
}

Vector2 normalized(Vector2 v){
    return v/len(v);
}

void normalize(Vector2 & v){
    v = normalized(v);
}

Vector2 get_normal(Vector2 v){
    return normalized({-v.y,v.x});
}

std::ostream & operator<<(std::ostream & os, Vector2 v){
    os<<"("<<v.x<<","<<v.y<<")";
    return os;
}

std::string to_string(Vector2 v){
    return ("("+std::to_string(v.x)+","+std::to_string(v.y)+")");
}

std::string to_string(Camera2D &camera){
    return( "Camera:{\n target: "+to_string(camera.target)+
        "\n offset: "+ to_string(camera.offset)+
        "\n zoom: "+ std::to_string(camera.zoom)+
        "\n rotation:"+ std::to_string(camera.rotation)+
        "\n}");
}

void draw_triangle_fan(ArrayV2 & vect, Vector2 pos = {0,0}, float rot = 0.0, float scale = 1.0, Color color = BLACK){
    size_t size = vect.size();
    if(size < 3) return;

    Vector2 a = vect.at(0),b,c;

    rlPushMatrix();
    rlTranslatef(pos.x,pos.y,0);
    rlRotatef(rot, 0,0,1);
    rlScalef(scale,scale,0);
    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r,color.g,color.b,color.a);
    for(size_t i=0; i<=size-3; ++i){
        b = vect.at(i+1);
        c = vect.at(i+2);
        rlVertex2f(a.x,a.y);
        rlVertex2f(b.x,b.y);
        rlVertex2f(c.x,c.y); 
    }
    rlEnd();
    rlPopMatrix();
}

void draw_ellipse(Vector2 center, Vector2 axses, float offset, float angle_deg, float thicc = 1.0, Color color = BLACK, int segments = 128){
    float angle = -DEG2RAD * angle_deg;
    Vector2 offset_vect = Vector2{-std::cos(angle), std::sin(angle)} * offset;
    rlPushMatrix();
    rlTranslatef(center.x + offset_vect.x, center.y + offset_vect.y, 0);
    rlRotatef(angle_deg, 0,0,1);    

    //DrawEllipseLines(0,0,axses.x, axses.y, BLACK);

    rlBegin(RL_TRIANGLES);
    rlColor4ub(color.r, color.g, color.b, color.a);
    float angle_per_seg = (2.0f*PI)/(float)segments;
    Vector2 p1, p2, p3, p4, n1, n2;

    p1 = {sinf(-angle_per_seg)*axses.x, cosf(-angle_per_seg)*axses.y};
    p2 = {sinf(0)*axses.x, cosf(0)*axses.y};
    n1 = get_normal(p2-p1);
    for (float a = 0; a < 2.0f*PI; a += angle_per_seg)
    {
        p1 = {sinf(a)*axses.x, cosf(a)*axses.y};
        p2 = {sinf(a+angle_per_seg)*axses.x, cosf(a+angle_per_seg)*axses.y};
        n2 = get_normal(p2-p1);
        p3 = p1 + n1 * thicc;
        p4 = p2 + n2 * thicc;

        n1 = n2;


        rlVertex2f(p2.x, p2.y);
        rlVertex2f(p1.x, p1.y);
        rlVertex2f(p3.x, p3.y);

        rlVertex2f(p3.x, p3.y);
        rlVertex2f(p4.x, p4.y);
        rlVertex2f(p2.x, p2.y);
    }
    rlEnd();

    rlPopMatrix();
}

int loops_observer = 0;
int test_loops_observer = 0;
double eccentric_anomaly(double mean_anomaly, double e, double tolerance = 1e-8) {
    double E = mean_anomaly; // Initial guess for E

    while (true && loops_observer < 10) {
        ++loops_observer;
        double f = E - e * sin(E) - mean_anomaly; // Function value
        double f_prime = 1.0 - e * cos(E); // Derivative value

        double delta = f / f_prime; // Update step

        E -= delta; // Update E

        if (std::abs(delta) < tolerance) {
            break; // Converged
        }
    }

    test_loops_observer = loops_observer;
    loops_observer = 0;

    return E;
}

void test_eccentric_anomaly(){
    double e = 0.4;
    for(double M = 0; M<=6.0; M+=0.5){
        double E = eccentric_anomaly(M,e);
        TraceLog(LOG_INFO, TextFormat("Keplers equation for M=%0.2f e=%0.2f -> E=%0.5f in %d loops", M, e, E, test_loops_observer));
        assert(test_loops_observer <= 6);
    }
}

Vector2 get_mouse_pos_world_space(Camera2D cam){
    Vector2 point = GetMousePosition();
    point = point - cam.offset; //screen center oriented
    point = point / cam.zoom;   //zoom adjusted
    point = point + cam.target; //back to world space
    return point;
}

bool has_point(const ArrayV2 & vect, Vector2 point){
    size_t size = vect.size();
    if(size < 3) return false;

    bool last_cross = false;
    bool cur_cross = false;
    Vector2 v, vp;
    for(size_t i=0; i<size; ++i){
        v = vect.at((i+1)%size) - vect.at(i);
        vp = point - vect.at(i);
        cur_cross = cross(v,vp) > 0;
        if(cur_cross != last_cross) return false;
        last_cross = cur_cross;
    }
    return true;
}

Vector2 rotate_point(Vector2 p, float rot){
    float s = std::sin(rot);
    float c = std::cos(rot);
    return {p.x*c - p.y*s,
            p.y*c + p.x*s};
}

Vector2 transform_point(Vector2 point, Vector2 pos, float rot, float scale = 1.0){
    point *= scale;
    point = rotate_point(point, rot);
    return point+pos;
}

float point_penetration(Vector2 a, Vector2 b, Vector2 point){
    Vector2 v = b-a;
    Vector2 n = get_normal(v);
    Vector2 p = point - a;
    return dot(p,n);
}

void test_point_penetration(){
    srand(time(NULL));
    Vector2 a = {(float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX};
    Vector2 b = {(float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX};
    Vector2 p = (a+b)*0.5;
    Vector2 n = get_normal(b-a);

    float result = point_penetration(a,b,p);
    //std::cout<<"a: "<<a<<" b: "<<b<<" p: "<<p<<"\n";
    //std::cout<<"result: "<<result<<"\n";
    assert(result < 0.001 && result > -0.001);

    result = point_penetration(a,b,p+n);
    assert(result < 1.001 && result > 0.999);

    result = point_penetration(a,b,p-n);
    assert(result < -0.999 && result > -1.001);
}

float deepest_line_penetration(Vector2 a, Vector2 b, const PartItem &p, size_t * index = nullptr){
    const ArrayV2 &v = get_part(p.id);
    size_t vsize = v.size();
    float deepest_penetration = std::numeric_limits<float>::lowest();
    Vector2 point;
    for(size_t i = 0; i<vsize; ++i){
        point = transform_point(v.at(i), p.position, 0.0, p.scale);
        float penetration = point_penetration(b,a, point);
        if(penetration > deepest_penetration){
            deepest_penetration = penetration;
            if(index) *index = i;
        }
    }
    return deepest_penetration;
}

// float least_penetration(const ArrayV2 &v1, const ArrayV2 &v2, size_t * index = nullptr){
//     size_t v1size = v1.size();
//     float least_penetration = -1.0;
//     Vector2 a,b;
//     size_t dp_index = -1;
//     float penetration = -1.0f;
//     for(size_t i=0; i<v1size; ++i){
//         a = v1.at(i);
//         b = v1.at((i+1)%v1size);
//         penetration = deepest_line_penetration(a,b, v2, &dp_index);
//         if(penetration < 0.0) return -1.0;
//         if(penetration < least_penetration){
//             least_penetration = penetration;
//             if(index) *index = dp_index;
//         }
//     }
//     return least_penetration;
// }

bool does_overlap(PartItem &p1, PartItem &p2){//const ArrayV2 &v1, const ArrayV2  &v2, Vector2 pos1, Vector2 pos2, float rot1 = 0.0f, float rot2 = 0.0f){
    const ArrayV2 &v1 = get_part(p1.id);
    const ArrayV2 &v2 = get_part(p2.id);
    size_t v1size = v1.size();
    size_t v2size = v2.size();
    Vector2 a,b;
    for(size_t i=0; i<v1size; ++i){
        a = v1.at(i);
        b = v1.at((i+1)%v1size);
        a = transform_point(a, p1.position, 0.0, p1.scale);
        b = transform_point(b, p1.position, 0.0, p1.scale);
        float penetration = deepest_line_penetration(a,b, p2);

        //std::cout<<" "<<penetration;
        if(penetration < 0.01f) return false;
    }

    for(size_t i=0; i<v2size; ++i){
        a = v2.at(i);
        b = v2.at((i+1)%v2size);
        a = transform_point(a, p2.position, 0.0, p2.scale);
        b = transform_point(b, p2.position, 0.0, p2.scale);
        float penetration = deepest_line_penetration(a,b, p1);

        //std::cout<<" "<<penetration;
        if(penetration < 0.01f) return false;
    }
    //std::cout<<"\n";
    
    return true;
}

struct CollisionInfo{
    bool is_colliding = {false};
    float penetration = {1.0f};
    Vector2 normal = {0.0f, 0.0f};
    Vector2 collision_point = {0.0f, 0.0f};
};

bool part_circle_collision(ArrayV2 &v, Vector2 position, float rotation_deg, float scale, Vector2 c_center, float c_radius, CollisionInfo & info){
    float rotation = DEG2RAD * rotation_deg;
    Vector2 t_center = c_center - position;         // reverse transform order into model space, tranlsate firstly
    t_center = rotate_point(t_center, -rotation);   // ten rotate
    t_center /= scale;                              // scale lastly 

    float t_radius = c_radius / scale;

    size_t v_size = v.size();
    size_t deepest_index = 0;
    float biggest_diff = 0.0;

    for(size_t i=0; i<v_size; ++i){
        Vector2 point = v.at(i);
        float distancesq = lensq(point-t_center);
        float diff = (t_radius * t_radius) - distancesq;

        if(diff > biggest_diff){
            biggest_diff = diff;
            deepest_index = i;
        }
    }

    info.is_colliding = biggest_diff > 0.0f;
    if(info.is_colliding){
        Vector2 point_offset = rotate_point( v.at(deepest_index) - t_center, rotation );
        info.collision_point = transform_point(v.at(deepest_index), {0,0}, rotation, 1.0);
        info.penetration = (len( point_offset ) - t_radius);
        info.normal =  normalized( point_offset );
    }
    return info.is_colliding;
}

bool vessel_circle_collision(Vessel &vessel, Vector2 c_center, float c_radius, CollisionInfo &info){
    CollisionInfo temp_info;
    info = temp_info;
    for(auto &part: vessel.parts){
        ArrayV2 &v = get_part(part.id);
        Vector2 position = rotate_point(part.position, DEG2RAD * vessel.rotation_deg) + vessel.position;
        if(part_circle_collision(v, position, vessel.rotation_deg, part.scale, c_center, c_radius, temp_info)){
            // info = temp_info;
            // info.collision_point += rotate_point(part.position, DEG2RAD * vessel.rotation_deg);
            // return true;
            if(temp_info.penetration < info.penetration){
                info = temp_info;
                info.collision_point += rotate_point(part.position, DEG2RAD * vessel.rotation_deg);
            }
        }
    }

    
    return( info.penetration < 0.0f);
}

void draw_vessel(std::vector<PartItem> &parts, Vector2 position, float rotation){
    for(PartItem & part : parts){
        ArrayV2 &vertices = get_part(part.id);
        Vector2 part_pos = rotate_point(part.position, DEG2RAD * rotation) + position;
        
        draw_triangle_fan(vertices, part_pos, rotation, part.scale, DARKBLUE);
        draw_triangle_fan(vertices, part_pos, rotation, part.scale*0.9, part.color);
    }
}

void DRAW_QUAD(Vector2 a1, Vector2 a2, Color ac, Vector2 b1, Vector2 b2, Color bc){
        rlColor4ub(ac.r,ac.g,ac.b,ac.a);
        rlVertex2f(a1.x,a1.y);
        rlColor4ub(bc.r,bc.g,bc.b,bc.a);
        rlVertex2f(b1.x,b1.y);
        rlVertex2f(b2.x,b2.y);
        rlColor4ub(ac.r,ac.g,ac.b,ac.a);
        rlVertex2f(a1.x,a1.y);
        rlColor4ub(bc.r,bc.g,bc.b,bc.a);
        rlVertex2f(b2.x,b2.y);
        rlColor4ub(ac.r,ac.g,ac.b,ac.a);
        rlVertex2f(a2.x,a2.y);
}

void get_point_pair(Vector2 &a, Vector2 &b, Vector2 pos, float width){
    a = pos;
    b = pos;
    a.x -= width * 0.5;
    b.x += width * 0.5;
}

void draw_exhaust(Vector2 pos, float rot, float width, float lenght, float cone = 1.0){
    int steps = 10;
    float step_lenght = lenght / float(steps);
    Color color = RED, prev_color = RED;

    Vector2 prev_a, prev_b;
    Vector2 a, b;
    float distortion = 0.0;
    static int distortion_offset = 0;
    distortion_offset = (distortion_offset+1)%1000;
    get_point_pair(prev_a, prev_b, pos, width);

    rlPushMatrix();
    //rlTranslatef(pos.x,pos.y,0);
    rlRotatef(rot, 0,0,1);
    //rlScalef(scale,scale,0);
    rlBegin(RL_TRIANGLES);
    
    for(int i=0; i<steps; ++i){
        color.a = 255 * ((float) (steps - i - 1) / (float) steps);
        distortion = std::sin((pos.y+(float)distortion_offset) * 123456789.0) * 0.05;
        pos.y += step_lenght;
        width *= cone;
        get_point_pair(a, b, pos, width);
        a.x -= distortion;
        b.x += distortion;

        DRAW_QUAD(prev_a, prev_b, prev_color, a, b, color);
        prev_a = a;
        prev_b = b;
        prev_color = color;
    }
    rlEnd();
    rlPopMatrix();
}

void test_has_point(){
    ArrayV2 & vect = get_part(POD);
    Vector2 point = {0,0};
    assert(has_point(vect,point));

    point = {10,0};
    assert(!has_point(vect,point));

    point = {-10,0};
    assert(!has_point(vect,point));
}

////////////////////////////////////////////////////
// PHYSICS
////////////////////////////////////////////////////

//offset and force are given in vessel space
void add_force(Vessel &vessel, Vector2 force, Vector2 offset){
    //get mass
    //add linear force
    //get inertia
    //add rotational force
}
#endif // COMMON_H