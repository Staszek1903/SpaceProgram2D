#ifndef EDITOR_SCREEN_H
#define EDITOR_SCREEN_H

#include "screen.h"
#include "common.h"

std::vector<PartItem> editable_parts;

class EditorScreen: public Screen{
    class Action{
    public:
        virtual void do_action() = 0;
        virtual void undo_action() = 0;
        virtual void print_str() { std::cout<<"UNKNOWN_ACTION"<< std::endl; };
    };
    
    class SpawnAction: public Action{
        PartID m_part_id;
        size_t m_index;
    public:
        virtual void print_str() override{
            std::cout<<"SPAWN_ACTION part_id:"<<m_part_id<<" index:"<<m_index<<std::endl;
        } 
        SpawnAction(PartID part_id, size_t index) : m_part_id(part_id), m_index(index) {}
        virtual void do_action() override;
        virtual void undo_action() override;
    };

    class MoveAction: public Action{
        size_t m_index;
        Vector2 m_src_pos;
        Vector2 m_dst_pos;
    public:
        MoveAction(size_t index, Vector2 src, Vector2 dst)
            :m_index(index), m_src_pos(src), m_dst_pos(dst) {}
        virtual void print_str() override{
            std::cout<<"MOVE_ACTION: index:"<<m_index<<" src_pos: "<<m_src_pos<<" dst_pos: "<<m_dst_pos<<std::endl;;
        }
        virtual void do_action() override;
        virtual void undo_action() override;
    };

    class DeleteAction: public Action{
        PartID m_part_id;
        Vector2 m_src_pos;
        size_t m_index;
    public:
        DeleteAction(PartID part_id, Vector2 src_pos, size_t index) 
            :m_part_id(part_id), m_src_pos(src_pos), m_index(index) {}
        virtual void print_str() override{
            std::cout<<"DELETE_ACTION: index:"<<m_index<<" src_pos: "<<m_src_pos<<" part_id: "<<m_part_id<<std::endl;
        }
        virtual void do_action() override;
        virtual void undo_action() override;
    };

    const Rectangle editor_area = {-5, -10, 10, 20};
    Camera2D camera = {{screenWidth*0.5, screenHeight*0.5}, {0.0, 0.0}, 0, 40.0f};

    Vector2 prev_camera_pos = {0,0};
    Vector2 prev_pan_mouse_pos = {0,0};

    PartItem * part_grabbed = nullptr;
    Vector2 prev_mouse_pos = {0,0};
    Vector2 prev_part_pos = {0,0};

    Image editor_grid;
    Rectangle editor_grid_rect = { 0.0f, 0.0f, editor_area.width*5.0f, editor_area.height*5.0f};
    Texture2D editor_grid_texture;

    std::vector<std::unique_ptr<Action>> actions;
    size_t actions_top = -1;

    bool rotate_part = true;

    std::string error_text = "";

    Vessel vessel_to_launch;

    void push_action(std::unique_ptr<Action> & action){
        ++actions_top;
        TraceLog(LOG_INFO, TextFormat("action queue top: %d", actions_top));
        action->print_str();

        if(actions.size() <= actions_top) actions.resize(actions.size() + 10);
        action.swap(actions.at(actions_top));
    }

    void undo_action(){
        TraceLog(LOG_INFO, TextFormat("action queue top: %d", actions_top));
        
        if(actions_top == (size_t)(-1)) return;
        actions.at(actions_top)->print_str();
        actions.at(actions_top)->undo_action();
        --actions_top;
    }

    void update_camera(){
        if(IsMouseButtonPressed(0)){
            prev_camera_pos = camera.target;
            prev_pan_mouse_pos = GetMousePosition();
        }
        if(!part_grabbed && IsMouseButtonDown(0)){
            Vector2 mouse_pos = GetMousePosition();
            Vector2 dp = (mouse_pos - prev_pan_mouse_pos);
            camera.target = prev_camera_pos - dp/camera.zoom; 
        }

        float wheel = GetMouseWheelMove() * float(!rotate_part);
        camera.zoom *= 1.0 + 0.1 * wheel;

        // if(IsKeyDown(KEY_RIGHT)) camera.target.x += 2;
        // if(IsKeyDown(KEY_LEFT)) camera.target.x -= 2;
        // if(IsKeyDown(KEY_DOWN)) camera.target.y += 2;
        // if(IsKeyDown(KEY_UP)) camera.target.y -= 2;
    }


    void update_linkages(){
        for(PartItem &parta: editable_parts){
            parta.links_joined = 0;
        }
        for(PartItem &parta: editable_parts){
            for(PartItem &partb: editable_parts){
                if(parta == partb) continue;
                check_linkage(parta, partb);
                std::cout<<"linkages: "<<parta.links_joined<<std::endl;
                std::cout<<"linkages: "<<partb.links_joined<<std::endl;
            }
        }

        // for(PartItem &parta: editable_parts){
        //     std::cout<<"linkages: "<<parta.links_joined<<std::endl;
        // }
    }

    //void check_part_connected

    bool are_all_parts_connected(){
        if(editable_parts.size() == 0) return false;

        std::vector<bool> connected_array;
        std::vector<size_t> index_to_check;
        connected_array.resize(editable_parts.size(), false);
        index_to_check.reserve(editable_parts.size());
        index_to_check.push_back(0);
        connected_array.at(0) = true;
        
        while(index_to_check.size() > 0){
            size_t cur_index = index_to_check.back();
            index_to_check.pop_back();
            PartItem & part_a = editable_parts.at(cur_index);
            for(size_t i = 0; i < editable_parts.size(); ++i){
                if( cur_index == i ) continue;
                if( connected_array.at(i) ) continue;
                if( check_linkage(part_a, editable_parts.at(i)) ){
                    connected_array.at(i) = true;
                    index_to_check.push_back(i);
                }
            }
        }

        bool result = true;
        for(bool c : connected_array) result = result && c;
        return result;
    }

    bool check_for_overlap(){
        bool result = false;
        for(PartItem & p: editable_parts)
            p.color.a = 255;

        for(PartItem &part : editable_parts){
            for(PartItem & p: editable_parts){
                if(&p == &part) continue;
                if(does_overlap(p,part)){
                    //TraceLog(LOG_INFO,"OVERLAP");
                    p.color.a = 128;
                    part.color.a = 128;
                    result = true;
                }
            }
        }

        return result;
    }

    size_t spawn_part(PartID p){
        PartItem item {p};
        item.type = part_type[p];
        item.info_index = part_info_index[p];

        editable_parts.push_back(item);
        size_t index = editable_parts.size()-1;
        std::unique_ptr<Action> spawn_action = std::make_unique<SpawnAction>(p, index);
        push_action(spawn_action);

        check_for_overlap();
        return index;
    }

    void generate_part_buttons(){
        int x = 5, y = 200;
        for(int i = 0; i<PartID::COUNT; ++i){
            buttons.push_back(Button{part_names.at(i), [&, i](){spawn_part((PartID)i);}, x,y+35*i,100,30});
        }
    }

    void draw_editable_parts(){
        for(auto &part: editable_parts){
            ArrayV2 &points = get_part(part.id);
            ArrayV2 &linkage = linkages.at(part.id);
            draw_triangle_fan(points, part.position, part.rotation_degree, part.scale, part.color);
            for(size_t i = 0; i<linkage.size(); ++i){
                Vector2 point = linkage.at(i);
                point = transform_point(point, part.position, part.rotation_degree * DEG2RAD, part.scale);
                DrawCircleV({point.x, point.y}, 5.0f / camera.zoom, ((part.links_joined>>i)&1)? GREEN : BLUE );
            }
            DrawCircleV(part.position, 5.0f / camera.zoom, RED );
            //DrawText(TextFormat("%d",part.links_joined), part.position.x, part.position.y, 8, BLACK);
            //draw_triangle_fan(points, part.position, 0.0, part.scale-2, part.color);
        }
    }

    bool is_in_editor_area(PartItem &part){
        ArrayV2 &p = get_part(part.id);
        for(Vector2 v: p){
            v = transform_point(v, part.position, 0.0, part.scale);
            if(!CheckCollisionPointRec(v,editor_area)) return false;
        }
        return true;
    }

    bool is_outside_editor_area(PartItem &part){
        ArrayV2 &p = get_part(part.id);
        for(Vector2 v: p){
            v = transform_point(v, part.position, 0.0, part.scale);
            if(CheckCollisionPointRec(v,editor_area)) return false;
        }
        return true;
    }

    void editor_area_check(PartItem &part, bool &inside, bool &outside){
        ArrayV2 &p = get_part(part.id);
        inside = true;
        outside = true;
        for(Vector2 v: p){
            v = transform_point(v, part.position, 0.0, part.scale);
            if(CheckCollisionPointRec(v,editor_area)) outside = false;
            else inside = false;
            if(!inside && !outside) return;
        }
    }
    
    void update_shortcuts(){
        if(!IsKeyDown(KEY_LEFT_CONTROL)) return;

        if(IsKeyPressed(KEY_Z)) undo_action();
    }


    float grid_snap = 0.1f;
    void update_parts(){

        if(rotate_part && part_grabbed){
            part_grabbed->rotation_degree += GetMouseWheelMove() * 5.0;
            std::cout<<part_grabbed->rotation_degree<<std::endl;
        }

        // CHECK PART DRAG
        Vector2 point = {0,0}, mouse_pos = get_mouse_pos_world_space(camera);
        bool inside = false, outside = false;
        if(part_grabbed) editor_area_check(*part_grabbed, inside, outside);
        if(part_grabbed && IsMouseButtonDown(0)){
            auto dp = (mouse_pos - prev_mouse_pos);
            part_grabbed->position = (prev_part_pos + dp);
            part_grabbed->position.x = 
                floor(part_grabbed->position.x/grid_snap)*grid_snap;
            part_grabbed->position.y = 
                floor(part_grabbed->position.y/grid_snap)*grid_snap;

            part_grabbed->color.a = inside?255:128;
            check_for_overlap();
        }else if(part_grabbed){
            PartItem part = *part_grabbed;
            auto iter = std::find(editable_parts.begin(), editable_parts.end(), *part_grabbed);
            size_t index = iter - editable_parts.begin();
            if(outside){
                editable_parts.erase( iter );
                std::unique_ptr<Action> del_action = std::make_unique<DeleteAction>(part.id, prev_part_pos, index);
                push_action(del_action);
            }else if(prev_part_pos != part.position) {
                std::unique_ptr<Action> move_action = std::make_unique<MoveAction>(index, prev_part_pos, part.position);
                push_action(move_action);
            }
            update_linkages();
            part_grabbed = nullptr;
        }

        // CHECK PART GRAB
        if(!IsMouseButtonPressed(0)) return;
        part_grabbed = nullptr;
        for(auto & part : editable_parts){
            point = mouse_pos - part.position;
            point = rotate_point(point, -part.rotation_degree * DEG2RAD);
            point = point / part.scale;

            if(has_point(get_part(part.id), point)){
                part_grabbed = &part;
                prev_part_pos = part_grabbed->position;
            }
        }
        prev_mouse_pos = mouse_pos;
    }

    void make_vessel(){
        vessel_to_launch.parts = {};
        // vessel_to_launch.engines = {};
        // vessel_to_launch.tanks = {};

        Vector2 origin = {0,0};
        for(PartItem & part : editable_parts) origin += part.position;
        origin /= editable_parts.size();
        
        vessel_to_launch.parts = editable_parts;
        for(PartItem & part : vessel_to_launch.parts) part.position -= origin;
        TraceLog(LOG_INFO, TextFormat("Origin: %s", to_string(origin).c_str()));

        // for(int i = 0; i < vessel_to_launch.parts.size(); ++i){
        //     PartItem & part = vessel_to_launch.parts[i];
        //     switch (part.type){
        //     case PartType::ENGINE: vessel_to_launch.engines.push_back(EnginePart{i}); break;
        //     case PartType::TANK: vessel_to_launch.tanks.push_back(FueltankPart{i}); break;
        //     }
        // }

        calculate_vessel_mass(vessel_to_launch);
        calculate_vessel_inertia(vessel_to_launch);
    }

public:
    EditorScreen(){
        editor_grid = GenImageChecked(editor_grid_rect.width, editor_grid_rect.height, 1, 1, LIGHTGRAY, RAYWHITE);
    }

    ~EditorScreen(){
        UnloadImage(editor_grid);
    }

    virtual void init(void* data = nullptr) override{
        editor_grid_texture = LoadTextureFromImage(editor_grid);
        generate_part_buttons();
        buttons.push_back(Button{"LAUNCH",
            [&](){
                if(editable_parts.size() == 0){
                    error_text = "CANT LAUNCH EMPTY VESSEL";
                    TraceLog(LOG_INFO, error_text.c_str());
                    return;
                }
                if(!are_all_parts_connected()){
                    error_text = "DETECTED UNATACHED PARTS";
                    TraceLog(LOG_INFO, error_text.c_str());
                    return;
                }

                if(check_for_overlap()){
                    error_text = "DETECTED OVERLAPING PARTS";
                    TraceLog(LOG_INFO, error_text.c_str());
                    return;
                }
    
                make_vessel();
                change_screen_handler(); 

            },
            screenWidth - 110, screenHeight-40, 100, 30});

        TraceLog(LOG_INFO, "EDITOR SCREEN INIT");
    }

    virtual void update() override{
        if(IsWindowResized()){
            TraceLog(LOG_INFO, "Window resized to (%d , %d)", GetScreenWidth(), GetScreenHeight());
            camera.offset = {(float)GetScreenWidth() * 0.5f, (float)GetScreenHeight() * 0.5f};
        }
        if(IsMouseButtonPressed(0)) error_text = "";

        rotate_part =  part_grabbed && IsKeyDown(KEY_R);

        update_gui();
        update_parts();
        update_camera();
        update_shortcuts();

        BeginDrawing();
            ClearBackground(RAYWHITE); 
            BeginMode2D(camera);
                DrawTexturePro(editor_grid_texture, editor_grid_rect, editor_area, {0,0}, 0.0, WHITE);
                //DrawTexture(editor_grid_texture, editor_area.x, editor_area.y, WHITE);      
                draw_editable_parts();
            EndMode2D();
            draw_gui();
            DrawText(TextFormat("has: %d\n%s\n%s",
                part_grabbed,
                to_string(camera).c_str(),
                to_string(get_mouse_pos_world_space(camera)).c_str()), 10, 40, 8, DARKGRAY);

            DrawText(error_text.c_str(), 100,100, 32, RED);
            
            DrawFPS(10, 10);
        EndDrawing();
    }

    virtual void* release() override{
        UnloadTexture(editor_grid_texture);
        clear_gui_data();
        TraceLog(LOG_INFO, "EDITOR SCREEN RELEASE");
        return &vessel_to_launch;
    }

};

void EditorScreen::SpawnAction::do_action(){
    editable_parts.push_back(PartItem{m_part_id, {200,100}});
    m_index = editable_parts.size()-1;
}
void EditorScreen::SpawnAction::undo_action(){
    if(m_index >= 0 && m_index < editable_parts.size())
        editable_parts.erase(editable_parts.begin()+m_index);
    
    m_index = -1;
}

void EditorScreen::DeleteAction::do_action(){
    if(m_index >= 0 && m_index < editable_parts.size())
        editable_parts.erase(editable_parts.begin()+m_index);
}
void EditorScreen::DeleteAction::undo_action(){
    editable_parts.emplace(editable_parts.begin() + m_index, PartItem{m_part_id, {200,100}});
    editable_parts.at(m_index).position = m_src_pos;
}

void EditorScreen::MoveAction::do_action(){
    if(m_index >= 0 && m_index < editable_parts.size())
        editable_parts.at(m_index).position = m_dst_pos;
}

void EditorScreen::MoveAction::undo_action(){
    if(m_index >= 0 && m_index < editable_parts.size())
        editable_parts.at(m_index).position = m_src_pos;
}

#endif // EDITOR_SCREEN_H