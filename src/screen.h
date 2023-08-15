#ifndef SCREEN_H
#define SCREEN_H

#include "common.h"

class Screen{
public:
    static std::function<void()> change_screen_handler;
    virtual void init(void* data = nullptr) = 0;
    virtual void update() = 0;
    virtual void* release() {return nullptr;};
};

std::function<void()> Screen::change_screen_handler;

#endif // SCREEN_H