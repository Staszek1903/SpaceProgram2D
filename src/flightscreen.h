#ifndef FLIGHT_SCREEN_H
#define FLIGHT_SCREEN_H

#include "screen.h"
#include "common.h"

struct VesselsPartsData{
    std::vector<EnginePart> engines{};
    std::vector<FueltankPart> tanks{};
};

class FlightScreen: public Screen{
    Vessel current_vessel;
    Camera2D camera = {{screenWidth*0.5, screenHeight*0.5}, {0.0, 0.0}, 0, 40.0f};
    VesselsPartsData vessels_parts_data;
    
    float earth_radius = 5000.0f;
    float gravity_magnitude = 50000000.0f;
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

    void update_engines(Vessel & vessel, float dt){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::ENGINE) continue;
            EnginePart &engine = vessels_parts_data.engines.at(item.data_index);
            EngineInfo &info = engine_info[item.info_index];
            if(engine.fuel_tank_index < 0) continue;;

            PartItem &tank_part = vessel.parts.at(engine.fuel_tank_index);
            FueltankPart &tank_data = vessels_parts_data.tanks.at(tank_part.data_index);
            FueltankInfo &tank_info = fuel_tank_info.at(tank_part.info_index);

            if(tank_data.fuel_amount < 0.0 || !engine.active){
                engine.active = false;
                continue;
            }

            float quotent_of_used_fuel = (info.consumptiom * dt * throttle) / tank_info.capacity;
            tank_data.fuel_amount -= quotent_of_used_fuel;
            
            Vector2 offset = rotate_point(item.position, vessel.rotation_deg * DEG2RAD);
            Vector2 force_normal = rotate_point( Vector2{0.0,-1.0}, vessel.rotation_deg * DEG2RAD + item.rotation_degree * DEG2RAD);

            if(vessel.mass < 0.0001){
                TraceLog(LOG_INFO, "VESSEL HAS NO MASS");
                assert(false);
            }

            add_force(vessel, force_normal * info.max_thrust * throttle * dt / vessel.mass, offset);
        }
    }

    void engines_gimbal(Vessel &vessel, float input, float dt){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::ENGINE) continue;
            EnginePart &engine = vessels_parts_data.engines.at(item.data_index);
            EngineInfo &info = engine_info[item.info_index];

            float prev_gimbal_state = engine.gimbal_state;
            engine.gimbal_state = std::lerp(engine.gimbal_state, input,  dt*2.0);
            float d_state = engine.gimbal_state - prev_gimbal_state;

            item.rotation_degree += info.gimbal_degrees * d_state;
        }
    }

    void update_tanks(Vessel &vessel){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::TANK) continue;
            FueltankPart &part = vessels_parts_data.tanks.at(item.data_index);
            FueltankInfo &info = fuel_tank_info.at(item.info_index);
            
            if(part.prev_tank_index < 0) continue;
            PartItem &item2 = vessel.parts.at(part.prev_tank_index);
            FueltankPart &part2 = vessels_parts_data.tanks.at(item2.data_index);
            FueltankInfo &info2 = fuel_tank_info.at(item.info_index);

            float missing_fuel = (1.0 - part.fuel_amount) * info.capacity;
            float available_fuel = part2.fuel_amount * info2.capacity;

            float fuel_transfer = (available_fuel > missing_fuel) ? missing_fuel : available_fuel;

            part.fuel_amount += fuel_transfer / info.capacity;
            part2.fuel_amount -= fuel_transfer / info2.capacity;
        }
    }

    void update_vessel(Vessel & vessel){
        float dt = GetFrameTime();

        if(vessel.inertia < 0.0001){
                TraceLog(LOG_INFO, "VESSEL HAS NO INERTIA");
                assert(false);
            }

        if(IsKeyDown(KEY_I)) get_sin_random() += dt;
        if(IsKeyDown(KEY_K)) get_sin_random() -= dt;
        float gimbal = 0.0;
        if(IsKeyDown(KEY_D)){
            //current_vessel.angular_vel_deg += 10.0/current_vessel.inertia;
            gimbal = -1.0;
        }
        if(IsKeyDown(KEY_A)){
            //current_vessel.angular_vel_deg -= 10.0/current_vessel.inertia;
            gimbal = 1.0;
        }
        if(IsKeyPressed(KEY_SPACE)){
            for(EnginePart &engine: vessels_parts_data.engines)
                engine.active = true;
        }
        engines_gimbal(current_vessel, gimbal, dt);

        //std::cout<<current_vessel.inertia<<std::endl;
        
        update_engines(vessel, dt);
        update_tanks(vessel);
        // current_vessel.velocity += Vector2{
        //     std::cos( DEG2RAD * (current_vessel.rotation_deg - 90.0f) ),
        //     std::sin( DEG2RAD * (current_vessel.rotation_deg - 90.0f) )
        //     } * 0.03f * throttle;

        float radiussq = lensq(vessel.position);
        Vector2 gravity = -normalized(vessel.position) * gravity_magnitude / radiussq;
        vessel.velocity += gravity*dt;
        vessel.position += vessel.velocity*dt;
        vessel.rotation_deg += vessel.angular_vel_deg*dt;
       // vessel.angular_vel_deg *= (1.0f-dt); //DAMPING

        //COLLISTION
        float bouncyness = 0.0;
        float friction = 0.1;
        if( vessel_circle_collision(current_vessel, {0,0}, earth_radius, collision_info) )
        // (get_part(current_vessel.parts.at(0).id), 
        //     current_vessel.position, current_vessel.rotation_deg,
        //     1.0, {0,0}, earth_radius, collision_info))
        {
                vessel.position -= collision_info.normal * collision_info.penetration;
                Vector2 normal_tangent = tangent_vector(collision_info.normal);

                float point_lin_vel = vessel.angular_vel_deg * DEG2RAD * len(collision_info.collision_point);
                Vector2 point_vel_normal = normalized(tangent_vector(collision_info.collision_point));
                Vector2 point_velocity = vessel.velocity + (point_vel_normal * point_lin_vel);

                float normal_velocity = dot(point_velocity, collision_info.normal);
                float tangent_velocity = dot(point_velocity, normal_tangent);

                add_force(current_vessel, collision_info.normal * -normal_velocity * vessel.mass * (1.0+bouncyness), collision_info.collision_point);
                add_force(current_vessel, normal_tangent * -tangent_velocity * vessel.mass * friction, collision_info.collision_point);
                //vessel.velocity -= collision_info.normal * dot(vessel.velocity, collision_info.normal);
                //vessel.velocity *= 0.99;
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
        //std::cout<<current_vessel.rotation_deg<<std::endl;

        for(PartItem &item : current_vessel.parts){
            switch (item.type)
            {
            case PartType::TANK:
            {
                FueltankPart &tank = vessels_parts_data.tanks.at(item.data_index);
                Vector2 position = transform_point(item.position, {0,0}, DEG2RAD * current_vessel.rotation_deg);
                DrawCircleSector(position, 0.35, 0, 360.0f * tank.fuel_amount , 16, GREEN);
                //tank.fuel_amount -= 0.001;
                
            }
            break;
            
            case  PartType::ENGINE:
            {
                EnginePart &engine = vessels_parts_data.engines.at(item.data_index);
                if(!engine.active) continue;
                EngineInfo info = engine_info.at(item.info_index);
                Vector2 exhaust_pos = rotate_point(item.position, current_vessel.rotation_deg * DEG2RAD);
                exhaust_pos += rotate_point({0, info.exhaust_offset}, (current_vessel.rotation_deg + item.rotation_degree) * DEG2RAD);
                draw_exhaust(exhaust_pos, current_vessel.rotation_deg + item.rotation_degree, info.exhaust_width, 2.0 * throttle, 1.1 );
            }
            break;

            case PartType::NONE: break;

            default:
                TraceLog(LOG_INFO, TextFormat("ERROR: Unknown part type: %s", std::to_string((int)item.type).c_str()));
                assert(false);
                break;
            }
        }

        // for(FueltankPart & part: vessels_parts_data.tanks){
        //     PartItem item = current_vessel.parts[part.part_item_index];
        //     Vector2 position = transform_point(item.position, {0,0}, DEG2RAD * current_vessel.rotation_deg);
        //     DrawCircleSector(position, 0.35, 0, 360.0f * part.fuel_amount , 16, GREEN);
        //     part.fuel_amount -= 0.001;
        // }

        // for(EnginePart & part: current_vessel.engines){
        //     if(!part.active) continue;
        //     PartItem item = current_vessel.parts[part.f];
        //     EngineInfo info = engine_info[item.info_index];
        //     Vector2 exhaust_pos = item.position;
        //     exhaust_pos.y += info.exhaust_offset;
        //     draw_exhaust(exhaust_pos, current_vessel.rotation_deg, info.exhaust_width, 2.0 * throttle, 1.1 );
        // }
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
            vessels_parts_data.engines.size(),
            vessels_parts_data.tanks.size(),
            to_string(current_vessel.position).c_str()),
            200,50,8,BLACK);
    }

    void init_new_vessel(){

    }

    int find_atached_tank(Vessel vessel, PartItem &part){
        for(int i=0; i<vessel.parts.size(); ++i){
            PartItem &tank = vessel.parts.at(i);
            if(&tank == &part) continue;
            if(tank.type != TANK) continue;
            if(!check_linkage(part,tank)) continue;
            if(tank.position.y > part.position.y) continue;
            return i;
        }
        return -1;
    }

public:
    virtual void init(void* data = nullptr) override {
        buttons.push_back(Button{"EDIT",
            [](){ change_screen_handler(); },
            screenWidth - 110, screenHeight-40, 100, 30});

        current_vessel = *(Vessel*)(data);
        current_vessel.position = {0, -earth_radius};//-(earth_radius*1.01f)};
        //current_vessel.velocity = {std::sqrt(gravity_magnitude / (earth_radius*1.01f)), 0.0f};

        for(PartItem &part: current_vessel.parts){
            std::cout<<"part_type: "<<part.type<<std::endl;
        }

        for(int i = 0; i < current_vessel.parts.size(); ++i){
            PartItem & part = current_vessel.parts[i];
            switch (part.type){
            case PartType::ENGINE:
            { 
                vessels_parts_data.engines.push_back(EnginePart());
                part.data_index = vessels_parts_data.engines.size() - 1;
                EnginePart &engine = vessels_parts_data.engines[part.data_index];
                engine.fuel_tank_index = find_atached_tank(current_vessel, part);
                if(engine.fuel_tank_index < 0) assert(false);
            }
                break;
            case PartType::TANK:
            {
                vessels_parts_data.tanks.push_back(FueltankPart());
                part.data_index = vessels_parts_data.tanks.size() - 1;
                FueltankPart &tank = vessels_parts_data.tanks[part.data_index];
                tank.prev_tank_index = find_atached_tank(current_vessel, part);
                if(tank.prev_tank_index < 0) TraceLog(LOG_INFO, "Prev_tank NOT found");
                else TraceLog(LOG_INFO, "Prev_tamk found");
            }
                break;
            case PartType::NONE: break;
            default:
                TraceLog(LOG_INFO, TextFormat("ERROR: Unknown part type: %s", std::to_string((int)part.type).c_str()));
                assert(false);
                break;
            }
        }

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