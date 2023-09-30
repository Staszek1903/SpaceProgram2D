#ifndef FLIGHT_SCREEN_H
#define FLIGHT_SCREEN_H

#include "screen.h"
#include "common.h"

class FlightScreen: public Screen{
    Vessel current_vessel;
    Camera2D camera = {{screenWidth*0.5, screenHeight*0.5}, {0.0, 0.0}, 0, 40.0f};
    
    float earth_radius = 100.0f;
    float gravity_magnitude = 10000.0f;
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

    float throttle = 0.0;

    CollisionInfo collision_info;

    void update_engines(Vessel & vessel){
        for(EnginePart & engine: vessel.engines){
            PartItem & item = vessel.parts.at(engine.part_item_index);
            Vector2 offset = rotate_point(item.position, vessel.rotation_deg * DEG2RAD);
            Vector2 force_normal = rotate_point( Vector2(0.0,-1.0), vessel.rotation_deg * DEG2RAD );

            add_force(vessel, force_normal * 0.03f * throttle, offset);
        }
    }

    void update_vessel(Vessel & vessel){
        float dt = GetFrameTime();

        if(IsKeyDown(KEY_D)) current_vessel.angular_vel_deg += 10.0/current_vessel.inertia;
        if(IsKeyDown(KEY_A)) current_vessel.angular_vel_deg -= 10.0/current_vessel.inertia;
        
        update_engines(vessel);
        // current_vessel.velocity += Vector2{
        //     std::cos( DEG2RAD * (current_vessel.rotation_deg - 90.0f) ),
        //     std::sin( DEG2RAD * (current_vessel.rotation_deg - 90.0f) )
        //     } * 0.03f * throttle;

        float radiussq = lensq(vessel.position);
        Vector2 gravity = -normalized(vessel.position) * gravity_magnitude / radiussq;
        vessel.velocity += gravity*dt;
        vessel.position += vessel.velocity*dt;
        vessel.rotation_deg += vessel.angular_vel_deg*dt;
        vessel.angular_vel_deg *= (1.0f-dt);

        //COLLISTION
        if( vessel_circle_collision(current_vessel, {0,0}, earth_radius, collision_info) )
        // (get_part(current_vessel.parts.at(0).id), 
        //     current_vessel.position, current_vessel.rotation_deg,
        //     1.0, {0,0}, earth_radius, collision_info))
        {
                vessel.position -= collision_info.normal * collision_info.penetration;
                vessel.velocity -= collision_info.normal * dot(vessel.velocity, collision_info.normal);
                vessel.velocity *= 0.99;
        }
    }

    void update_telemetry(){
        float velsq = lensq(current_vessel.velocity);
        float radius = len(current_vessel.position);       
        //float angular_momentum = cross(current_vessel.position, current_vessel.velocity);

        Vector2 eccentricity_vect = (current_vessel.position * (velsq - gravity_magnitude/radius) - 
            current_vessel.velocity*dot(current_vessel.velocity,current_vessel.position)) /
            gravity_magnitude;

        eccentricity = len(eccentricity_vect);

 

        float E = velsq/2.0f - gravity_magnitude/radius;
        sm_axis = -gravity_magnitude/(2*E);

        apoapsis = sm_axis*(1+eccentricity);
        periapsis = sm_axis*(1-eccentricity);

        orbital_period = 2*PI* sqrt(sm_axis*sm_axis*sm_axis / gravity_magnitude);
        peri_argument = std::atan2(eccentricity_vect.y, eccentricity_vect.x);

        float vessel_argument = std::atan2(current_vessel.position.y, current_vessel.position.x);
        anomaly = vessel_argument - peri_argument;
        anomaly += (anomaly < 0.0) * 2*PI; // normalize anomaly [0, 2PI]
        mean_anomaly = anomaly - eccentricity * std::sin(anomaly);
        time_to_peri = orbital_period - orbital_period * (mean_anomaly/(2*PI));

        float oposite_anomaly = mean_anomaly - PI;
        oposite_anomaly += 2*PI * (oposite_anomaly<0.0f);
        time_to_apo = orbital_period - orbital_period * (oposite_anomaly/(2*PI));
    }

    void update_throttle(){
        if(throttle < 1.0f && IsKeyDown(KEY_LEFT_SHIFT)) throttle += GetFrameTime();
        if(throttle > 0.0f && IsKeyDown(KEY_LEFT_CONTROL)) throttle -= GetFrameTime();
        if(throttle > 1.0f) throttle = 1.0f;
        if(throttle < 0.0f) throttle = 0.0f;
        if(IsKeyPressed(KEY_Z)) throttle = 1.0f;
        if(IsKeyPressed(KEY_X)) throttle = 0.0f;
    }

    void draw_throttle(){
        float h = 200;
        float w = 10;
        float top = screenHeight - 250;
        float bottom = screenHeight - 250 + h;
        float x_pos = screenWidth - 50;
        DrawRectangle(x_pos - w/2.0f, top , w, h, {0,0,0,100});

        float guage_h = bottom - ((bottom - top) * throttle);
        DrawRectangle(x_pos - w, guage_h -5, w*2, 10, {0,0,0,100});
    }

    void draw_ship(){
        draw_vessel(current_vessel.parts, {0,0}, current_vessel.rotation_deg);

        for(FueltankPart & part: current_vessel.tanks){
            PartItem item = current_vessel.parts[part.part_item_index];
            Vector2 position = transform_point(item.position, {0,0}, DEG2RAD * current_vessel.rotation_deg);
            DrawCircleSector(position, 0.35, 0, 360.0f * part.fuel_amount , 16, GREEN);
            part.fuel_amount -= 0.001;
        }

        for(EnginePart & part: current_vessel.engines){
            if(!part.active) continue;
            PartItem item = current_vessel.parts[part.part_item_index];
            EngineInfo info = engine_info[item.info_index];
            Vector2 exhaust_pos = item.position;
            exhaust_pos.y += info.exhaust_offset;
            draw_exhaust(exhaust_pos, current_vessel.rotation_deg, info.exhaust_width, 2.0 * throttle, 1.1 );
        }
    }

    void draw_planet(){
        DrawCircleSector(-current_vessel.position, earth_radius, 0, 360,64, BROWN);
    }

    void draw_trajectory(){
        Vector2 peri_vect{std::cos(peri_argument), std::sin(peri_argument)};
        DrawLineV(-current_vessel.position, -current_vessel.position + peri_vect * periapsis, RED);
        DrawLineV(-current_vessel.position, -current_vessel.position - peri_vect * apoapsis, BLUE);
        //DrawEllipseLines(-current_vessel.position.x, current_vessel.position.y, , BLACK);
        draw_ellipse(-current_vessel.position, 
            {sm_axis, sm_axis*(float)sqrt(1-eccentricity*eccentricity)}, 
            sm_axis-periapsis, RAD2DEG * peri_argument,
            2.5f/camera.zoom, {0,0,0,64});

        //DRAW COLLSION
        DrawCircleV(collision_info.collision_point, 4.0f/camera.zoom, PURPLE);
        DrawLineEx(collision_info.collision_point, collision_info.collision_point + collision_info.normal * 20.0f/camera.zoom,
        3.0f/camera.zoom, GREEN);
        DrawLineEx(collision_info.collision_point, collision_info.collision_point + collision_info.normal * collision_info.penetration,
        2.0f/camera.zoom, RED);
    }

    void draw_telemetry(){
        DrawText(TextFormat("V: %0.2f\nH: %0.2f\nApo: %0.2f\nPer: %0.2f\nEcc: %0.2f\nT: %0.2f\nTp: %0.2f\nTa: %0.2f\nCollision: %d p: %0.4f\nAngle: %0.2f",
            len(current_vessel.velocity),
            len(current_vessel.position),
            apoapsis,
            periapsis,
            eccentricity,
            orbital_period,
            time_to_peri,
            time_to_apo,
            collision_info.is_colliding,
            collision_info.penetration,
            current_vessel.rotation_deg),
            10,50,24,
            BLACK
            );

        DrawText(TextFormat("DEBUG INFO:\nparts: %d\nengines: %d\nfueltanks: %d\n pos: %s",
            current_vessel.parts.size(),
            current_vessel.engines.size(),
            current_vessel.tanks.size(),
            to_string(current_vessel.position).c_str()),
            200,50,8,BLACK);
    }

public:
    virtual void init(void* data = nullptr) override {
        buttons.push_back(Button{"EDIT",
            [](){ change_screen_handler(); },
            screenWidth - 110, screenHeight-40, 100, 30});

        
        // std::string a = typeid(data).name();
        // std::string b = typeid(Vessel*).name();
        // TraceLog(LOG_INFO, TextFormat("Incoming data type: \"%s\" expected type is: \"%s\"", a.c_str(), b.c_str()));
        current_vessel = *(Vessel*)(data);
        current_vessel.position = {0, -(earth_radius*1.1f)};
        current_vessel.velocity = {std::sqrt(gravity_magnitude / (earth_radius*1.1f)), 0.0f};

        TraceLog(LOG_INFO, "FLIGHT SCREEN INIT");
    }

    virtual void update() override {
        float wheel = GetMouseWheelMove();
        camera.zoom *= 1.0 + 0.1 * wheel;

        update_vessel(current_vessel);
        update_telemetry();
        update_gui();
        update_throttle();

        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode2D(camera);
            draw_planet();
            draw_ship();
            draw_trajectory();
        EndMode2D();
            draw_gui();
            draw_telemetry();
            draw_throttle();
            DrawFPS(10, 10);
        EndDrawing();
    }

    virtual void* release() override{
        clear_gui_data();
        current_vessel.parts.clear();
        TraceLog(LOG_INFO, "FLIGHT SCREEN RELEASE");
        return nullptr;
    }
};

#endif // FLIGHT_SCREEN_H