#ifndef COMMON_H
#define COMMON_H

#include "raylib.h"
#include "rlgl.h"
#include "gui.h"
#include <iostream>
#include <stdio.h>
#include <vector>
#include <functional>
#include <math.h>
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
    TANK,
    SEPARATE,
    LIFE_SUPPORT
};

struct PartItem
{
    PartID id;
    Vector2 position = {0.0,0.0};
    float rotation_degree = 0.0;
    float scale = 1.0;
    Color color = GRAY;//{GRAY.r,GRAY.g,GRAY.b, 128};
    int links_joined = {0};

    PartType type = {PartType::NONE};
    int data_index = {-1};
    int info_index = {-1};
};

bool operator==(const PartItem & a, const PartItem & b){
    return (&a == &b);
}

struct EngineInfo{
    float max_thrust;
    float consumptiom;
    float exhaust_width;
    float exhaust_offset;
    float gimbal_degrees;
};

struct EnginePart{
    int fuel_tank_index = {-1};
    bool active = {false};
    float gimbal_state = {0.0};
};

struct FueltankInfo{
    float capacity = {0};
};

struct FueltankPart{
    float fuel_amount = {1.0};
    int prev_tank_index = {-1};
};

struct OrbitalElemets{
    float apoapsis = 0.0;
    float periapsis = 0.0f;
    float eccentricity = 0.0;
    float orbital_period = 0.0;
    float peri_argument = 0.0;
    float sm_axis = 0.0;
    float anomaly = 0.0;
    float mean_anomaly = 0.0;
    float time_to_peri = 0.0;
    float time_to_apo = 0.0;
};

struct Vessel{
    std::vector<PartItem> parts{};

    Vector2 position{0.0f,0.0f};
    Vector2 velocity{0.0f,0.0f}; 
    float rotation_deg{0.0f};
    float angular_vel_deg{0.0f};

    float mass = 1.0;
    float inertia = 1.0;

    float throttle = 0.0;
    float broad_phase_radius = 0.0;
    OrbitalElemets elements{0};
};

struct CelestialBody{
    float radius = 5000.0f;
    float gravity = 50000000.0f;
};

struct CollisionInfo{
    bool is_colliding = {false};    //selfexplanatory
    float penetration = {1.0f};     //length of vector from a point overlaping another body to surface of this body
    Vector2 normal = {0.0f, 0.0f};  //surface normal vector from a body into part/vessel in world space
    Vector2 collision_point = {0.0f, 0.0f}; //vector from part/vessel origin to overlaping point in world space
    PartItem * colliding_part = nullptr;    //part of reference
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
    PartType::LIFE_SUPPORT,
    PartType::TANK,
    PartType::ENGINE,
    PartType::ENGINE,
    PartType::SEPARATE
};

std::vector<int> part_info_index = {
    -1,
    0,
    0,
    1,
    -1
};

std::vector<float> part_dry_mass = {
    0.3,
    1.0,
    0.2,
    0.2,
    0.01
};

std::vector<ArrayV2> parts_poly ={
    {{-0.2,-0.2},{-0.6,0.8},{0.6,0.8},{0.2,-0.2}},
    {{-0.6,-1},{-0.6,1},{0.6,1},{0.6,-1}},
    {{-0.2,-0.6},{-0.4,0.6},{0.4,0.6},{0.2,-0.6}},
    {{-0.2,-0.2},{-0.4,0.2},{0.4,0.2},{0.2,-0.2}},
    {{-0.6,-0.2},{-0.6,0.0},{0.6,0.0},{0.6,-0.2}}
};

std::vector<float> parts_broadphase_radius;

std::vector<ArrayV2> linkages ={
    {/*{0.0,0.8}*/},
    {/*{0.0,1.0},*/ {0.0,-1.1}},
    {/*{0.0,-0.6},*/{0.0,-0.7}},
    {/*{0.0,-0.2},*/{0.0,-0.3}},
    {/*{0.0,0.0},*/ {0.0,-0.3}},
};

std::vector<EngineInfo> engine_info = {
    {15.0, 1.0, 0.8, 0.6, 15.0},
    {5.0, 0.1, 0.8, 0.2, 15.0},
};

std::vector<FueltankInfo> fuel_tank_info = {
    {50.0}
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

auto & get_polygon(PartID part){
    try{
        return parts_poly.at(part);
    }catch(std::exception &e){
        std::cout<<("EXCEPTION at "+std::to_string(__LINE__)+" in "+__FILE__)<<std::endl;
        assert(false);
    }
}

float get_broadphase_radius(PartID part){
    try{
        return parts_broadphase_radius.at(part);
    }catch(std::exception &e){
        std::cout<<("EXCEPTION at "+std::to_string(__LINE__)+" in "+__FILE__)<<std::endl;
        assert(false);
    }
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

Vector2 tangent_vector(Vector2 v){
    return Vector2{-v.y, v.x};
}

Vector2 project_vector(Vector2 v1, Vector2 v2){
    Vector2 nv2 = normalized(v2);
    return nv2 * dot(v1, nv2);
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

std::string to_string(ArrayV2 &v){
    std::string s = "ArrayV2: [ ";
    for(Vector2 &vec: v)
        s += to_string(vec) + ", ";
    s += "]";
    return s;
}

float sign(float a){
    return (-1.0 * (a<0.0) + 1.0 * (a>0.0));
}

// float clamp(float number, float min, float max){
//     bool bigger = number > max;
//     bool lesser = number < min;
//     return ( bigger * max + lesser * min + (!bigger && ! lesser) * number);
// }

// void test_clamp(){
//     assert(clamp(0.123, -1, 1) == 0.123);
//     assert(clamp(-3.1, -1, 1) == -1.0);
//     assert(clamp(123, -1, 1) == 1);
// }

// float lerp(float from, float to, float amount){
//     amount = clamp(amount, 0.0, 1.0);
//     return (from + (to - from) * amount);
// }

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

Vector2 inv_transform_point(Vector2 point, Vector2 pos, float rot, float scale = 1.0){
    point -= pos;
    point = rotate_point(point, -rot);
    return point/scale;
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
    const ArrayV2 &v = get_polygon(p.id);
    size_t vsize = v.size();
    float deepest_penetration = std::numeric_limits<float>::lowest();
    Vector2 point;
    for(size_t i = 0; i<vsize; ++i){
        point = transform_point(v.at(i), p.position, p.rotation_degree * DEG2RAD, p.scale);
        float penetration = point_penetration(b,a, point);
        if(penetration > deepest_penetration){
            deepest_penetration = penetration;
            if(index) *index = i; 
        }
    }
    return deepest_penetration;
}

CollisionInfo collision_info_placeholder;
bool does_overlap(PartItem &p1, PartItem &p2, CollisionInfo &info = collision_info_placeholder){//const ArrayV2 &v1, const ArrayV2  &v2, Vector2 pos1, Vector2 pos2, float rot1 = 0.0f, float rot2 = 0.0f){
    info.is_colliding = false;
    info.penetration = std::numeric_limits<float>::max();
    const ArrayV2 &v1 = get_polygon(p1.id);
    const ArrayV2 &v2 = get_polygon(p2.id);
    size_t v1size = v1.size();
    size_t v2size = v2.size();
    Vector2 a,b;
    size_t index  = 0.0;
    for(size_t i=0; i<v1size; ++i){
        a = v1.at(i);
        b = v1.at((i+1)%v1size);
        a = transform_point(a, p1.position, p1.rotation_degree * DEG2RAD, p1.scale);
        b = transform_point(b, p1.position, p1.rotation_degree * DEG2RAD, p1.scale);

        float penetration = deepest_line_penetration(a,b, p2, &index);
        if(penetration < info.penetration){
            info.penetration = penetration;
            info.normal = tangent_vector(normalized(b - a));
            info.collision_point = v2.at(index);
            info.colliding_part = &p2; 
        }
        //std::cout<<" "<<penetration;
        if(penetration < 0.01f) return false;
    }

    for(size_t i=0; i<v2size; ++i){
        a = v2.at(i);
        b = v2.at((i+1)%v2size);
        a = transform_point(a, p2.position, p2.rotation_degree * DEG2RAD, p2.scale);
        b = transform_point(b, p2.position, p2.rotation_degree * DEG2RAD, p2.scale);

        float penetration = deepest_line_penetration(a,b, p1, &index);
        if(penetration < info.penetration){
            info.penetration = penetration;
            info.normal = tangent_vector(normalized(b - a));
            info.collision_point = v1.at(index);
            info.colliding_part = &p1; 
        }
        //std::cout<<" "<<penetration;
        if(penetration < 0.01f) return false;
    }
    //std::cout<<"\n";
    info.is_colliding = true;
    //info.collision_point = transform_point(info.collision_point, info.colliding_part->position, info.colliding_part->rotation_degree, info.colliding_part->scale);
    info.collision_point = rotate_point(info.collision_point, info.colliding_part->rotation_degree * DEG2RAD);
    return true;
}

void calculate_broadphase_radiuses(){
    for(ArrayV2 &poly: parts_poly){
        float max_radius = 0.0;
        for(Vector2 v: poly){
            float radius = len(v);
            if(max_radius < radius) max_radius = radius;
        }
        parts_broadphase_radius.push_back(max_radius);
    }
}

void calculate_vessel_broadphase_radius(Vessel &vessel){
    float max_radius = 0.0;
    for(PartItem &item: vessel.parts){
        float radius = len(item.position);
        radius += get_broadphase_radius(item.id);
        if(max_radius < radius) max_radius = radius;
    }
    vessel.broad_phase_radius = max_radius;
}

bool part_broadphase_collision(/*Vessel &vessel1, Vessel &vessel2,*/ PartItem &item1, PartItem &item2){
    // Vector2 part_pos1 = transform_point(item1.position, vessel1.position, vessel1.rotation_deg * DEG2RAD, 1.0);
    // Vector2 part_pos2 = transform_point(item2.position, vessel2.position, vessel2.rotation_deg * DEG2RAD, 1.0);
    //Vector2 dist_vect = part_pos2 - part_pos1;
    Vector2 dist_vect = item2.position - item1.position;
    float radius1 = get_broadphase_radius(item1.id);
    float radius2 = get_broadphase_radius(item2.id);
    float radius_sum = radius1 + radius2; 
    return (lensq(dist_vect) < radius_sum * radius_sum);
}

bool vessel_broadphase_collision(Vessel &vessel1, Vessel &vessel2){
    float radius_sum = vessel1.broad_phase_radius + vessel2.broad_phase_radius;
    return(lensq(vessel2.position - vessel1.position) < radius_sum * radius_sum);
}

bool part_circle_collision(ArrayV2 &v, Vector2 position, float rotation_deg, float scale, Vector2 c_center, float c_radius, CollisionInfo & info){
    float rotation = DEG2RAD * rotation_deg;
    Vector2 t_center = c_center - position;         // reverse transform order into model space, tranlsate firstly
    t_center = rotate_point(t_center, -rotation);   // then rotate
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
        Vector2 point_offset = rotate_point( v.at(deepest_index) - t_center, rotation ); // point from center of circle used to calculate normal and penetration in world space
        info.collision_point = transform_point(v.at(deepest_index), {0,0}, rotation, 1.0); // point from center of part in world space
        info.penetration = (len( point_offset ) - t_radius);
        info.normal =  normalized( point_offset );
    }
    return info.is_colliding;
}

bool vessel_circle_collision(Vessel &vessel, Vector2 c_center, float c_radius, CollisionInfo &info){
    CollisionInfo temp_info;
    info = temp_info;
    for(auto &part: vessel.parts){
        ArrayV2 &v = get_polygon(part.id);
        Vector2 position = rotate_point(part.position, DEG2RAD * vessel.rotation_deg) + vessel.position;
        if(part_circle_collision(v, position, vessel.rotation_deg + part.rotation_degree, part.scale, c_center, c_radius, temp_info)){
            if(temp_info.penetration < info.penetration){
                info = temp_info;
                info.collision_point += rotate_point(part.position, DEG2RAD * vessel.rotation_deg);
                info.colliding_part = &part;
            }
        }
    }
    return( info.penetration < 0.0f);
}

bool vessel_vessel_collision(Vessel &vessel1, Vessel &vessel2, CollisionInfo &info){
    info.is_colliding = false;
    PartItem p1, p2;
    for(PartItem &part1: vessel1.parts){
        p1 = part1;
        p1.position = transform_point(p1.position, vessel1.position, vessel1.rotation_deg * DEG2RAD);
        p1.rotation_degree += vessel1.rotation_deg;
        for(PartItem &part2: vessel2.parts){
            p2 = part2;
            p2.position = transform_point(p2.position, vessel2.position, vessel2.rotation_deg * DEG2RAD);
            p2.rotation_degree += vessel2.rotation_deg;
            if(part_broadphase_collision(p1,p2) && 
            does_overlap(p1, p2, info)){
                info.is_colliding = true;
                if(info.colliding_part == &p2){
                    info.normal = -info.normal;
                    info.collision_point += (p1.position - p2.position );
                }

                info.collision_point += p1.position - vessel1.position;
                return true;
            }
        }
    }

    return false;
}

void collision_position_correction(Vessel &vessel, CollisionInfo &info){
    vessel.position -= info.normal * info.penetration;
}

void collision_position_correction(Vessel &vessel, Vessel &other, CollisionInfo &info){
    float v_vel = len(vessel.velocity);
    float o_vel = len(other.velocity);
    float v_q = v_vel / (v_vel + o_vel);
    float o_q = o_vel / (v_vel + o_vel);
    vessel.position += info.normal * info.penetration * v_q;
    other.position  -= info.normal * info.penetration * o_q;
}

Vector2 get_point_velocity(Vessel &vessel, Vector2 &point){ // point relative to vessel in world space, result also in world space
    float point_lin_vel = vessel.angular_vel_deg * DEG2RAD * len(point);
    Vector2 point_vel_normal = normalized(tangent_vector(point));
    Vector2 point_velocity = vessel.velocity + (point_vel_normal * point_lin_vel);
    return point_velocity;
}

void draw_vessel(std::vector<PartItem> &parts, Vector2 position, float rotation){
    for(PartItem & part : parts){
        ArrayV2 &vertices = get_polygon(part.id);
        Vector2 part_pos = rotate_point(part.position, DEG2RAD * rotation) + position;
        
        //DrawCircleV(part_pos, get_broadphase_radius(part.id), Color(0,255,0,128));
        draw_triangle_fan(vertices, part_pos, rotation + part.rotation_degree, part.scale, DARKBLUE);
        draw_triangle_fan(vertices, part_pos, rotation + part.rotation_degree, part.scale*0.9, part.color);
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

float sinus_random = 5.454;
float& get_sin_random(){ std::cout<<"sin rand: "<<sinus_random<< "\n"; return sinus_random; }

void draw_exhaust(Vector2 pos, float rot, float width, float lenght, float cone = 1.0){
    int steps = 10;
    float step_lenght = lenght / float(steps);
    Color color = RED, prev_color = RED;

    Vector2 prev_a, prev_b;
    Vector2 a, b;
    float distortion = 0.0;
    static int distortion_offset = 0;
    distortion_offset = (distortion_offset+1)%1000;
    get_point_pair(prev_a, prev_b, {0,0}, width);

    rlPushMatrix();
    rlTranslatef(pos.x,pos.y,0);
    rlRotatef(rot, 0,0,1);
    rlBegin(RL_TRIANGLES);
    
    float y = 0.0;
    for(int i=0; i<steps; ++i){
        color.a = 255 * ((float) (steps - i - 1) / (float) steps);
        distortion = std::sin((y+(float)distortion_offset) * sinus_random) * 0.05;
        y += step_lenght;
        width *= cone;
        get_point_pair(a, b, {0.0,y}, width);
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
    ArrayV2 & vect = get_polygon(POD);
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

Vector2 get_centroid(Vector2 a, Vector2 b, Vector2 c){
    return (a + b + c) / 3.0f;
}

float triangle_area(Vector2 a, Vector2 b, Vector2 c){
    return std::abs(0.5 * (((b.x-a.x)*(c.y-a.y))-((c.x-a.x)*(b.y-a.y))));
}

float polygon_area(ArrayV2 &v){
    float area = 0.0;
    size_t prev = v.size() - 1;
    for(size_t i =0; i < v.size(); ++i){
        area += triangle_area({0.0,0.0}, v.at(prev), v.at(i));
        prev = i;
    }

    return area;
}

void test_polygon_area(){
    ArrayV2 v = {{-1,-1}, {-1, 1}, {1,1}, {1, -1}};
    assert(polygon_area(v) == 4.0f);
}

Vector2 get_centroid(ArrayV2 &v){
    assert(v.size() >= 3);
    Vector2 a{v.at(0)}, b{}, c{}, sum{}, tri_cen{};
    float area = 0.0, tri_area = 0.0;
    for(size_t i = 1; i < v.size()-1; ++i){
        b = v.at(i);
        c = v.at(i+1);
        std::cout<<"a: "<<a<<" b: "<<b<<" c: "<<c<<std::endl;
        tri_area = triangle_area(a,b,c);
        tri_cen = get_centroid(a,b,c);
        std::cout<<"centroid: "<<tri_cen<<" area: "<<tri_area<<std::endl;
        area += tri_area;
        sum += tri_cen * tri_area;
    }

    return sum / area;
}

void test_get_centroid(){
    ArrayV2 v = {{-1,-1}, {-1, 1}, {1,1}, {1, -1}};
    Vector2 centroid = get_centroid(v), expected{};
    std::cout<<"Test centroid: "<<centroid<< " expected: "<< expected<<std::endl;
    assert(centroid == expected);
}
 
const double prec = 0.1f;
void adjust_parts_centroids(){
    for(size_t i = 0; i < parts_poly.size(); ++i){
        ArrayV2 & a = get_polygon((PartID)i);
        Vector2 c = get_centroid(a);
        float bottom = 0.0;
        for(Vector2 & v: a){
            v -= c;
            v = Vector2{float(round(v.x/prec)), float(round(v.y/prec))} * prec;
            if(v.y < bottom) bottom = v.y;
        }
        ArrayV2 & links = linkages.at(i);
        for(Vector2 & l: links){
            l -= c;
            l = Vector2{float(round(l.x/prec)), float(round(l.y/prec))} * prec;
        }
        for(EngineInfo & info: engine_info){
            info.exhaust_offset = -bottom;
            //info.exhaust_offset = float(round(info.exhaust_offset/prec)) * prec;
        }
    }
}

void calculate_vessel_mass(Vessel &vessel){
    vessel.mass = 0.0;
    for(PartItem &item: vessel.parts){
        vessel.mass += part_dry_mass[item.id];
    }

    if(vessel.mass < 0.0001 || vessel.mass > 1e9){
        TraceLog(LOG_INFO, "MASS CALCULATION RESULT ABNORMAL");
        assert(false);
    }
}

float get_inertia(ArrayV2 &v, float mass){
    float inertia = 0.0;
    float area = polygon_area(v);
    float tri_mass = 0.0, tri_area = 0.0, tri_inertia = 0.0;

    size_t prev = v.size()-1;
    for(size_t i = 0; i< v.size(); ++i){
        tri_area = triangle_area({0,0}, v.at(prev), v.at(i));
        tri_mass = mass * tri_area / area;
        tri_inertia = tri_mass * (lensq(v.at(prev)) + lensq(v.at(i)) + dot(v.at(prev), v.at(i)));
        prev = i;
        inertia += tri_inertia;
    }
    return inertia;
}

void calculate_vessel_inertia(Vessel &vessel){
    vessel.inertia = 0.0;
    for(PartItem & item: vessel.parts){
        float inertia = get_inertia(get_polygon(item.id), part_dry_mass.at(item.id));
        Vector2 offset = item.position;
        float distsq = lensq(offset);
        vessel.inertia += inertia + vessel.mass*distsq; // Parallel axis theorem
    }

    if(vessel.inertia < 0.0001 || vessel.inertia > 1e9){
        TraceLog(LOG_INFO, "INERTIA CALCULATION RESULT ABNORMAL");
        assert(false);
    }
}

bool check_linkage(PartItem &a, PartItem &b){
    ArrayV2 &la = linkages.at(a.id);
    ArrayV2 &lb = linkages.at(b.id);
    ArrayV2 &va = get_polygon(a.id);
    ArrayV2 &vb = get_polygon(b.id);

    for(size_t i=0; i<la.size(); ++i){
        Vector2 lp = la.at(i);
        lp = transform_point(lp, a.position, a.rotation_degree * DEG2RAD, a.scale);
        //lp = transform_point(lp, -b.position, -b.rotation_degree * DEG2RAD, b.scale);
        lp -= b.position;
        lp = rotate_point(lp, -b.rotation_degree * DEG2RAD);
        lp *= 1.0f/b.scale;

        if(has_point(vb,lp)){
            a.links_joined = a.links_joined | ((1<<i));
            return true;
        }
    }

    for(size_t i=0; i<lb.size(); ++i){
        Vector2 lp = lb.at(i);
        lp = transform_point(lp, b.position, b.rotation_degree * DEG2RAD, b.scale);
        //lp = transform_point(lp, -a.position, -a.rotation_degree * DEG2RAD, a.scale);
        lp -= a.position;
        lp = rotate_point(lp, -a.rotation_degree * DEG2RAD);
        lp *= 1.0f/a.scale;
        std::cout<<to_string(lp)<<"\n"<<to_string(la)<<std::endl;
        if(has_point(va,lp)){
            b.links_joined = b.links_joined | ((1<<i));
            return true;
        }
    }

    return false;
}

void update_vessel_origin(Vessel & vessel){
    Vector2 position_sum = {0,0};
    float mass_sum = 0;
    for(PartItem &item: vessel.parts){
        size_t index = static_cast<size_t>(item.id);
        float part_mass = part_dry_mass.at(index);
        position_sum += item.position * part_mass;
        mass_sum += part_dry_mass.at(item.id);
    }
    Vector2 origin = position_sum/mass_sum;
    for(PartItem &item: vessel.parts){
        item.position -= origin;
    }
}

//offset and force are given in world space with origin at vessel position
void add_force(Vessel &vessel, Vector2 force, Vector2 offset){
    assert(!isnan(vessel.angular_vel_deg));

    if(lensq(force) < 0.0000001) return;
    //add linear force
    vessel.velocity += force / vessel.mass;

    //add rotational force
    Vector2 force_normal = normalized(force);
    Vector2 projected_offset = project_vector(offset, force);
    Vector2 tangent_offset = offset - projected_offset;
    float offcenter = len(tangent_offset) * sign(cross(tangent_offset, force_normal));

    vessel.angular_vel_deg += RAD2DEG * (offcenter*len(force))/vessel.inertia;

    if(isnan(vessel.angular_vel_deg)){
        std::cout<<"AAAAAAAAAAAA"<<std::endl;
    }
}
#endif // COMMON_H