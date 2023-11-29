#ifndef FLIGHT_SCREEN_H
#define FLIGHT_SCREEN_H

#include "screen.h"
#include "common.h"

struct VesselsPartsData{
    std::vector<EnginePart> engines{};
    std::vector<FueltankPart> tanks{};

    EnginePart &get_engine(size_t index){
        assert(index >=0 && index < engines.size());
        return engines.at(index);
    }

    FueltankPart &get_tank(size_t index){
        assert(index >=0 && index < tanks.size());
        return tanks.at(index);
    }
};

class FlightScreen: public Screen{
    std::vector<Vessel> vessels;
    int active_vessel_index = 0;
    Camera2D camera = {{screenWidth*0.5, screenHeight*0.5}, {0.0, 0.0}, 0, 40.0f};
    VesselsPartsData parts_data;
    CelestialBody earth {5000.0f, 50000000.0f};

    float throttle = 0.0;

    CollisionInfo collision_info;
    int warp = 0;

    // DO NOT ERASE DATA BCUS IT SEVERS ALL INDEXES IN PART ITEM
    // MARK AS EMPT SOMEHOW AND OVERWRITE WHEN CERATING NEXT PART INSTEAD
    // void remove_part_data(PartItem &item){
    //     switch (item.type)
    //     {
    //     case ENGINE:
    //         parts_data.engines.erase(parts_data.engines.begin() + item.data_index);
    //         break;
    //     case TANK:
    //         parts_data.tanks.erase(parts_data.tanks.begin() + item.data_index);
    //     default:
    //         break;
    //     }
    // }

    Vessel& get_active_vessel(){
        assert(active_vessel_index >= 0 && active_vessel_index < vessels.size());
        return vessels.at(active_vessel_index);
    }

    void update_engines(Vessel & vessel, float dt){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::ENGINE) continue;
            EnginePart &engine = parts_data.get_engine(item.data_index);
            EngineInfo &info = engine_info.at(item.info_index);
            if(engine.fuel_tank_index < 0) continue;;

            PartItem &tank_part = vessel.parts.at(engine.fuel_tank_index);
            FueltankPart &tank_data = parts_data.get_tank(tank_part.data_index);
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

            add_force(vessel, force_normal * info.max_thrust * vessel.throttle * dt / vessel.mass, offset);
        }
    }

    void engines_gimbal(Vessel &vessel, float input, float dt){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::ENGINE) continue;
            EnginePart &engine = parts_data.get_engine(item.data_index);
            EngineInfo &info = engine_info.at(item.info_index);

            float prev_gimbal_state = engine.gimbal_state;
            engine.gimbal_state = std::lerp(engine.gimbal_state, input,  dt*2.0);
            float d_state = engine.gimbal_state - prev_gimbal_state;

            item.rotation_degree += info.gimbal_degrees * d_state;
        }
    }

    void update_tanks(Vessel &vessel){
        for(PartItem &item: vessel.parts){
            if(item.type != PartType::TANK) continue;
            FueltankPart &part = parts_data.get_tank(item.data_index);
            FueltankInfo &info = fuel_tank_info.at(item.info_index);
            
            if(part.prev_tank_index < 0) continue;
            PartItem &item2 = vessel.parts.at(part.prev_tank_index);
            FueltankPart &part2 = parts_data.get_tank(item2.data_index);
            FueltankInfo &info2 = fuel_tank_info.at(item.info_index);

            float missing_fuel = (1.0 - part.fuel_amount) * info.capacity;
            float available_fuel = part2.fuel_amount * info2.capacity;

            float fuel_transfer = (available_fuel > missing_fuel) ? missing_fuel : available_fuel;

            part.fuel_amount += fuel_transfer / info.capacity;
            part2.fuel_amount -= fuel_transfer / info2.capacity;
        }
    }

    const float REACTION_WHEELS_FORCE = 200.0f;
    void vessel_input(Vessel &vessel, float dt){
        if(IsKeyDown(KEY_I)) get_sin_random() += dt;
        if(IsKeyDown(KEY_K)) get_sin_random() -= dt;
        float gimbal = 0.0;
        if(IsKeyDown(KEY_D)){
            gimbal = -1.0;
        }
        if(IsKeyDown(KEY_A)){

            gimbal = 1.0;
        }
        if(IsKeyPressed(KEY_SPACE)){
            for(EnginePart &engine: parts_data.engines)
                engine.active = true;
        }
        engines_gimbal(vessel, gimbal, dt);
        vessel.angular_vel_deg -= gimbal * REACTION_WHEELS_FORCE * dt / vessel.inertia;
    }

    void update_vessel(Vessel & vessel, float dt){
        if(vessel.inertia < 0.0001){
                TraceLog(LOG_INFO, "VESSEL HAS NO INERTIA");
                assert(false);
            }
        
        update_engines(vessel, dt);
        update_tanks(vessel);

        float radiussq = lensq(vessel.position);
        Vector2 gravity = -normalized(vessel.position) * earth.gravity / radiussq;
        vessel.velocity += gravity*dt;
        vessel.position += vessel.velocity*dt;
        vessel.rotation_deg += vessel.angular_vel_deg*dt;
       // vessel.angular_vel_deg *= (1.0f-dt); //DAMPING

    }

    void planet_collision(Vessel &vessel, CelestialBody &planet){
        CollisionInfo info;
        Vector2 planet_center = {0,0};
        float bouncyness = 0.0;
        float friction = 0.1;
        if( vessel_circle_collision(vessel, planet_center, earth.radius, info) )
        {
                collision_position_correction(vessel, info);
                Vector2 normal_tangent = tangent_vector(info.normal);
                Vector2 point_velocity = get_point_velocity(vessel, info.collision_point);
                float normal_velocity = dot(point_velocity, info.normal);
                float tangent_velocity = dot(point_velocity, normal_tangent);

                add_force(vessel, info.normal * -normal_velocity * vessel.mass * (1.0+bouncyness), info.collision_point);
                add_force(vessel, normal_tangent * -tangent_velocity * vessel.mass * friction, info.collision_point);
                //vessel.velocity -= info.normal * dot(vessel.velocity, info.normal);
                //vessel.velocity *= 0.99;
        }
    }

    void vessel2vessel_collision(Vessel& vessel, Vessel &other){
        //CollisionInfo info;
        float bouncyness = 0.0;
        float friction = 0.1;
        if(vessel_broadphase_collision(vessel, other) && vessel_vessel_collision(vessel, other, collision_info)){
            collision_position_correction(vessel, other, collision_info);

            Vector2 v_point = collision_info.collision_point;
            Vector2 o_point = collision_info.collision_point - (other.position - vessel.position);
            Vector2 v_point_velocity = get_point_velocity(vessel, v_point);
            Vector2 o_point_velocity = get_point_velocity(other, o_point);
            Vector2 collision_velocity = v_point_velocity - o_point_velocity;
            float normal_velocity = dot(collision_velocity, collision_info.normal) * 0.5;
            Vector2 normal_tangent = tangent_vector(collision_info.normal);
            float tangent_velocity = dot(collision_velocity, normal_tangent);

            float force_aplied = normal_velocity * vessel.mass;
            TraceLog(LOG_DEBUG, "collision force: %f", force_aplied);
            if(force_aplied > 2.0 || force_aplied < -2.0) 
                TraceLog(LOG_DEBUG, "collision_velocity: %s", to_string(collision_velocity).c_str());

            add_force(vessel, -collision_info.normal * force_aplied * (1.0+bouncyness), v_point);
            add_force(vessel, -normal_tangent * tangent_velocity * vessel.mass * friction, v_point);
        }
    }

    void update_telemetry(Vessel &vessel, CelestialBody &body, OrbitalElemets &elements){
        float velsq = lensq(vessel.velocity);
        float radius = len(vessel.position);       
        //float angular_momentum = cross(vessel.position, vessel.velocity);

        Vector2 eccentricity_vect = (vessel.position * (velsq - earth.gravity/radius) - 
            vessel.velocity*dot(vessel.velocity,vessel.position)) /
            earth.gravity;

        elements.eccentricity = len(eccentricity_vect);

 

        float E = velsq/2.0f - earth.gravity/radius;
        elements.sm_axis = -earth.gravity/(2*E);

        elements.apoapsis = elements.sm_axis*(1+elements.eccentricity);
        elements.periapsis = elements.sm_axis*(1-elements.eccentricity);

        elements.orbital_period = 2*PI* sqrt(elements.sm_axis*elements.sm_axis*elements.sm_axis / earth.gravity);
        elements.peri_argument = std::atan2(eccentricity_vect.y, eccentricity_vect.x);

        float vessel_argument = std::atan2(vessel.position.y, vessel.position.x);
        elements.anomaly = vessel_argument - elements.peri_argument;
        elements.anomaly += (elements.anomaly < 0.0) * 2*PI; // normalize anomaly [0, 2PI]
        elements.mean_anomaly = elements.anomaly - elements.eccentricity * std::sin(elements.anomaly);
        elements.time_to_peri = elements.orbital_period - elements.orbital_period * (elements.mean_anomaly/(2*PI));

        float oposite_anomaly = elements.mean_anomaly - PI;
        oposite_anomaly += 2*PI * (oposite_anomaly<0.0f);
        elements.time_to_apo = elements.orbital_period - elements.orbital_period * (oposite_anomaly/(2*PI));
    }

    void update_throttle(){
        if(throttle < 1.0f && IsKeyDown(KEY_LEFT_SHIFT)) throttle += GetFrameTime();
        if(throttle > 0.0f && IsKeyDown(KEY_LEFT_CONTROL)) throttle -= GetFrameTime();
        if(throttle > 1.0f) throttle = 1.0f;
        if(throttle < 0.0f) throttle = 0.0f;
        if(IsKeyPressed(KEY_Z)) throttle = 1.0f;
        if(IsKeyPressed(KEY_X)) throttle = 0.0f;
    }

    void get_next_stage_indices(std::vector<bool> &next_stage, Vessel &vessel, PartItem &item){
        std::vector<size_t> to_check;
        to_check.reserve(vessel.parts.size());
        next_stage.clear();
        next_stage.resize(vessel.parts.size());

        int stage_first_index = find_atached_part(vessel, item);
        to_check.push_back(stage_first_index); 
        next_stage.at(stage_first_index) = true; 

        const int guard = 10000;
        for(int guard_counter = 0; to_check.size() > 0 && guard_counter < guard; ++guard_counter)
        {
            size_t check_item_index = to_check.back();
            PartItem &check_item = vessel.parts.at(check_item_index);
            to_check.pop_back();
            for(size_t i = 0; i<vessel.parts.size(); ++i){
                PartItem &next_item = vessel.parts.at(i);
                if(&item != &next_item &&
                    check_linkage(check_item, next_item) &&
                    !next_stage.at(i)){
                        to_check.push_back(i);
                        next_stage.at(i) = true;
                }
            }
        }
    }

    void stage_separation(Vessel &vessel, PartItem &item){
        if(item.type != PartType::SEPARATE) return;
        std::vector<bool> next_stage;
        get_next_stage_indices(next_stage, vessel, item);

        // vessels.push_back(Vessel());
        // Vessel & debris = vessels.back();
        // debris = vessel;

        std::vector <PartItem> vessel_parts;
        //std::vector <PartItem> debris_parts;
        vessel_parts.reserve(next_stage.size());
       // debris_parts.reserve(next_stage.size());
        for(size_t i =0; i<next_stage.size(); ++i)
            if(next_stage.at(i))vessel_parts.push_back(vessel.parts.at(i));
           // else debris_parts.push_back(vessel.parts.at(i));
        vessel.parts = std::move(vessel_parts);
        //debris.parts = std::move(debris_parts);
        
        Vector2 old_pos = vessel.parts.at(0).position;
        update_vessel_origin(vessel);
        Vector2 new_pos = vessel.parts.at(0).position;
        vessel.position += (old_pos - new_pos);
        calculate_vessel_mass(vessel);
        calculate_vessel_inertia(vessel);
        calculate_vessel_broadphase_radius(vessel);

        // old_pos = debris.parts.at(0).position;
        // update_vessel_origin(debris);
        // new_pos = debris.parts.at(0).position;
        // debris.position += (new_pos - old_pos);
        // calculate_vessel_mass(debris);
        // calculate_vessel_inertia(debris);
        // calculate_vessel_broadphase_radius(debris);
    }

    void update_click_action(Vessel &vessel){
        if(!IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) return;
        Vector2 mouse_pos = get_mouse_pos_world_space(camera);
        //TraceLog(LOG_INFO, TextFormat("CLICKED %s", to_string(mouse_pos).c_str()));
        mouse_pos = rotate_point(mouse_pos, -vessel.rotation_deg * DEG2RAD);

        for(PartItem &item: vessel.parts){
            Vector2 trans_mouse = mouse_pos - item.position;
            trans_mouse = rotate_point(trans_mouse, -item.rotation_degree * DEG2RAD);
            trans_mouse /= item.scale;
            if(has_point(get_polygon(item.id), trans_mouse)){
                TraceLog(LOG_INFO, TextFormat("clicked type: %d", item.type));
                switch (item.type)
                {
                case PartType::ENGINE:
                    {
                        EnginePart &engine = parts_data.get_engine(item.data_index);
                        engine.active = !engine.active;
                    }
                    break;
                
                case PartType::SEPARATE:
                    stage_separation(vessel, item);
                    return;
                
                default:
                    break;
                }
            }
        }
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

    void draw_exhausts(Vessel &vessel, Vector2 active_pos){
        for(PartItem &item : vessel.parts){
            if(item.type != PartType::ENGINE) continue;
            EnginePart &engine = parts_data.get_engine(item.data_index);
            if(!engine.active) continue;
            EngineInfo info = engine_info.at(item.info_index);
            Vector2 exhaust_pos = rotate_point(item.position, vessel.rotation_deg * DEG2RAD);
            exhaust_pos += rotate_point({0, info.exhaust_offset}, (vessel.rotation_deg + item.rotation_degree) * DEG2RAD);
            exhaust_pos += vessel.position - active_pos;
            draw_exhaust(exhaust_pos, vessel.rotation_deg + item.rotation_degree, info.exhaust_width, 2.0 * vessel.throttle, 1.0 );
        }
    }

    void draw_tank_guages(Vessel &vessel){
        for(PartItem &item : vessel.parts){
            if(item.type != PartType::TANK) continue;
            FueltankPart &tank = parts_data.get_tank(item.data_index);
            Vector2 position = transform_point(item.position, {0,0}, DEG2RAD * vessel.rotation_deg);
            DrawCircleSector(position, 0.35, 0, 360.0f * tank.fuel_amount , 16, GREEN);
            //tank.fuel_amount -= 0.001;

            // case PartType::NONE: break;
            // case PartType::SEPARATE: break;
            // case PartType::LIFE_SUPPORT: break;

            // default:
            //     TraceLog(LOG_INFO, TextFormat("ERROR: Unknown part type: %s", std::to_string((int)item.type).c_str()));
            //     assert(false);
            //     break;
            
        }
    }

    void draw_planet(Vessel &vessel, CelestialBody &body){
        DrawCircleSector(-vessel.position, body.radius, 0, 360, 0, BROWN);
    }

    void draw_trajectory(Vector2 pos, OrbitalElemets &elements){
        Vector2 peri_vect{std::cos(elements.peri_argument), std::sin(elements.peri_argument)};
        DrawLineV(pos, pos + peri_vect * elements.periapsis, RED);
        DrawLineV(pos, pos - peri_vect * elements.apoapsis, BLUE);
        //DrawEllipseLines(pos.x, vessel.position.y, , BLACK);
        draw_ellipse(pos, 
            {elements.sm_axis, elements.sm_axis*(float)sqrt(1-elements.eccentricity*elements.eccentricity)}, 
            elements.sm_axis-elements.periapsis, RAD2DEG * elements.peri_argument,
            2.5f/camera.zoom, {0,0,0,64});

        //DRAW COLLSION
        DrawCircleV(collision_info.collision_point, 4.0f/camera.zoom, PURPLE);
        DrawLineEx(collision_info.collision_point, collision_info.collision_point + collision_info.normal * 20.0f/camera.zoom,
        3.0f/camera.zoom, GREEN);
        DrawLineEx(collision_info.collision_point, collision_info.collision_point + collision_info.normal * collision_info.penetration,
        2.0f/camera.zoom, RED);
    }

    void draw_telemetry(Vessel &vessel,  OrbitalElemets &elements){
        DrawText(TextFormat("V: %0.2f\nH: %0.2f\nApo: %0.2f\nPer: %0.2f\nEcc: %0.2f\nT: %0.2f\nTp: %0.2f\nTa: %0.2f\nCollision: %d p: %0.4f\nAngle: %0.2f",
            len(vessel.velocity),
            len(vessel.position) - earth.radius,
            elements.apoapsis - earth.radius,
            elements.periapsis - earth.radius,
            elements.eccentricity,
            elements.orbital_period,
            elements.time_to_peri,
            elements.time_to_apo,
            collision_info.is_colliding,
            collision_info.penetration,
            vessel.rotation_deg),
            10,50,24,
            BLACK
            );

        DrawText(TextFormat("DEBUG INFO:\nparts: %d\nengines: %d\nfueltanks: %d\n pos: %s\n camera:%s\n warp: %d",
            vessel.parts.size(),
            parts_data.engines.size(),
            parts_data.tanks.size(),
            to_string(vessel.position).c_str(),
            to_string(camera).c_str(),
            warp),
            200,50,8,BLACK);            
    }

    int find_atached_tank(Vessel &vessel, PartItem &part){
        // for(int i=0; i<vessel.parts.size(); ++i){
        //     PartItem &tank = vessel.parts.at(i);
        //     if(&tank == &part) continue;
        //     if(tank.type != TANK) continue;
        //     if(!check_linkage(part,tank)) continue;
        //     if(tank.position.y > part.position.y) continue;
        //     return i;
        // }
        int index = find_atached_part(vessel, part);
        return ((index >= 0 && vessel.parts.at(index).type == PartType::TANK)? index : -1);
    }

    int find_atached_part(Vessel &vessel, PartItem &item){
        Vector2 linkage_pos = linkages.at(item.id).at(0);
        linkage_pos = transform_point(linkage_pos, item.position, item.rotation_degree * DEG2RAD, item.scale);
        for(size_t i = 0; i<vessel.parts.size(); ++i){
            PartItem &item_check = vessel.parts.at(i);
            Vector2 link_pos_check = inv_transform_point(linkage_pos, item_check.position, item_check.rotation_degree * DEG2RAD, item_check.scale);
            if(has_point(get_polygon(item_check.id), link_pos_check)) return i;
        }
        return -1;
    }

public:
    virtual void init(void* data = nullptr) override {
        buttons.push_back(Button{"EDIT",
            [](){ change_screen_handler(); },
            screenWidth - 110, screenHeight-40, 100, 30});

        
        vessels.push_back(*(Vessel*)(data));
        active_vessel_index = vessels.size() - 1;
        Vessel &new_vessel = get_active_vessel();
        new_vessel.position = {0, -earth.radius};//-(earth.radius*1.01f)};
        //new_vessel.velocity = {std::sqrt(earth.gravity / (earth.radius*1.01f)), 0.0f};

        for(PartItem &part: new_vessel.parts){
            std::cout<<"part_type: "<<part.type<<std::endl;
        }

        for(int i = 0; i < new_vessel.parts.size(); ++i){
            PartItem & part = new_vessel.parts.at(i);
            switch (part.type){
            case PartType::ENGINE:
            { 
                parts_data.engines.push_back(EnginePart());
                part.data_index = parts_data.engines.size() - 1;
                EnginePart &engine = parts_data.get_engine(part.data_index);
                engine.fuel_tank_index = find_atached_tank(new_vessel, part);
                if(engine.fuel_tank_index < 0) assert(false);
            }
                break;
            case PartType::TANK:
            {
                parts_data.tanks.push_back(FueltankPart());
                part.data_index = parts_data.tanks.size() - 1;
                FueltankPart &tank = parts_data.get_tank(part.data_index);
                tank.prev_tank_index = find_atached_tank(new_vessel, part);
                if(tank.prev_tank_index < 0) TraceLog(LOG_INFO, "Prev_tank NOT found");
                else TraceLog(LOG_INFO, "Prev_tamk found");
            }
                break;
            case PartType::NONE: break;
            case PartType::SEPARATE: break;
            case PartType::LIFE_SUPPORT: break;
            default:
                TraceLog(LOG_INFO, TextFormat("ERROR: Unknown part type: %s", std::to_string((int)part.type).c_str()));
                assert(false);
                break;
            }

            // for(PartItem &part: new_vessel.parts){
            //     if(part.data_index == -1 || part.info_index == -1){
            //         assert(false);
            //     }
            // }
        }

        TraceLog(LOG_INFO, "FLIGHT SCREEN INIT");
    }

    
    virtual void update() override {
        float dt = GetFrameTime();
        float wheel = GetMouseWheelMove();
        camera.zoom *= 1.0 + 0.1 * wheel;

        dt *= std::pow(2.0, (double)warp);
        if(IsKeyPressed(KEY_COMMA) && warp > 0) --warp;
        if(IsKeyPressed(KEY_PERIOD)) ++warp;

        Vessel & active_vessel = get_active_vessel();

        vessel_input(active_vessel, dt);
        for(Vessel& vessel: vessels){
            update_vessel(vessel, dt);
            planet_collision(vessel,earth);
            update_telemetry(vessel, earth, vessel.elements);
            for(Vessel &other_vessel: vessels){
                if(&vessel == &other_vessel) continue;
                vessel2vessel_collision(vessel, other_vessel);
            }
        }
        update_gui();
        update_throttle();
        active_vessel.throttle = throttle;
        update_click_action(active_vessel);

        BeginDrawing();
        ClearBackground(RAYWHITE);
        BeginMode2D(camera);
            draw_planet(active_vessel, earth);
            for(Vessel& vessel: vessels){
                Vector2 pos = vessel.position - active_vessel.position;
                //DrawCircleV(pos, vessel.broad_phase_radius, Color{0,0,255,128});
                if(camera.zoom > 5.0) draw_vessel(vessel.parts, pos, vessel.rotation_deg);
                else DrawCircleSector(pos, 5.0/camera.zoom, 0,360, 4, GRAY);
                draw_trajectory(-active_vessel.position, vessel.elements);
                draw_exhausts(vessel, active_vessel.position);
            }
            draw_tank_guages(active_vessel);
        EndMode2D();
            draw_gui();
            draw_telemetry(active_vessel, active_vessel.elements);
            draw_throttle();
            DrawFPS(10, 10);
        EndDrawing();
    }

    virtual void* release() override{
        clear_gui_data();
        // current_vessel.parts.clear();
        // parts_data.engines.clear();
        // parts_data.tanks.clear();
        TraceLog(LOG_INFO, "FLIGHT SCREEN RELEASE");
        return nullptr;
    }
};

#endif // FLIGHT_SCREEN_H