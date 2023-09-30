#include "editorscreen.h"
#include "flightscreen.h"


void tests(){
    test_has_point();
    test_point_penetration();
    test_eccentric_anomaly();
    test_get_centroid();
    test_polygon_area();
}



int main()
{
    EditorScreen editor_screen;
    FlightScreen flight_screen;
    std::vector<Screen*> screens = {
        &editor_screen,
        &flight_screen
    };

    size_t current_screen = 0;
    Screen::change_screen_handler = [&](){
        void *transfer_data = screens[current_screen]->release();
        current_screen = (current_screen + 1)%screens.size();
        screens[current_screen]->init(transfer_data);
    };


    bool pause = false;




    tests();

    adjust_parts_centroids();

    SetConfigFlags(FLAG_MSAA_4X_HINT | FLAG_WINDOW_RESIZABLE);
    InitWindow(screenWidth, screenHeight, "KSP");        
    SetTargetFPS(60);


    // TraceLog(LOG_INFO, TextFormat("int: %d long: %d long long: %d long long int: %d",
    //                         sizeof(int), sizeof(long), sizeof(long long), sizeof(long long int)));


    screens[current_screen]->init(); 
    while (!WindowShouldClose()){
        pause = IsKeyPressed(KEY_P) ^ pause;
        //TraceLog(LOG_INFO, TextFormat("pause: %d P:%d", pause, IsKeyPressed(KEY_P)));
        if(!pause) screens[current_screen]->update();
        else{ 
            BeginDrawing(); 
            DrawText("PAUSED", 100,100, 32, RED);
            EndDrawing();
        }
    }
    screens[current_screen]->release();

    CloseWindow();                 
    return 0;
}






















