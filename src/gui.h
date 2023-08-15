#ifndef GUI_H
#define GUI_H

#include "raylib.h"
#include <functional>
#include <string>
#include <iostream>

#define GUI_COLOR LIGHTGRAY
#define GUI_HOOVER_COLOR SKYBLUE
#define GUI_CLICKED_COLOR BLUE
#define GUI_BORDER_COLOR DARKGRAY
#define GUI_BORDER_WIDTH 3
#define GUI_TEXT_SIZE 16

class Button{
public:
    std::string text{"BUTTON"};
    std::function<void ()> handler;
    int x{100},y{100},w{100},h{30};

    enum{
        NORMAL,
        HOOVER,
        CLICKED,
    };

    uint8_t state{NORMAL};
    void update(Vector2 mouse_pos, bool mouse_button)
    {
        int prev_state = state;
        state = has_mouse(mouse_pos)? ((mouse_button && state == HOOVER) ? CLICKED : HOOVER) : NORMAL;
        if(handler && prev_state == CLICKED && state != CLICKED)
            handler();
    }

    void draw()
    {
        DrawRectangle(x,y,w,h, GUI_BORDER_COLOR);
        Color c[] = {GUI_COLOR, GUI_HOOVER_COLOR, GUI_CLICKED_COLOR};
        Color curent_color = c[state];
        DrawRectangle(x+GUI_BORDER_WIDTH,y+GUI_BORDER_WIDTH,
            w - 2*GUI_BORDER_WIDTH, h - 2* GUI_BORDER_WIDTH,
            curent_color);

        Font font = GetFontDefault();
        Vector2 text_size = MeasureTextEx(font, text.c_str(), GUI_TEXT_SIZE, 1);

        DrawText(text.c_str(), (x+w/2)-text_size.x/2, (y+h/2)-text_size.y/2, GUI_TEXT_SIZE, DARKGRAY);
    }

protected:
    bool has_mouse(Vector2 mouse_pos)
    {
        return mouse_pos.x > x && mouse_pos.x < (x+w) && mouse_pos.y > y && mouse_pos.y < (y+h);
    }
};

class TextInput : public Button
{
public:
    bool active = false;

    void update(Vector2 mouse_pos, bool mouse_button, int key)
    {
        //std::cout<<"aa";
        int prev_state = state;
        state = has_mouse(mouse_pos)? mouse_button? CLICKED : HOOVER : NORMAL;
        if(prev_state == CLICKED && state != CLICKED){
            active = !active;
            //std::cout<<active<<std::endl;
        }

        if(!active || !key) return;
        if(key == KEY_BACKSPACE) text = text.substr(0,text.size()-1);
        else if((key >= (int)'0' && key <= (int)'9') ||
            //key >= (int)'a' && key <= (int)'z' ||
            ( key >= (int)'A' && key <= (int)'Z' )||
            key == (int)'-' || key == (int)'_' )
        {
            text += (char)key;
        }
    }

    void draw(){
        Button::draw();
        if(!active) return;
        double time = GetTime();
        if(time - (double)(int)time > 0.5) return;
        DrawRectangleV({(float)(x+3),(float)(y+3)}, {10,20}, BLACK);
    }
};

#endif // GUI_H