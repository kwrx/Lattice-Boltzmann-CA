 
//
// MIT License

// Copyright (c) 2020 Antonino Natale

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
 


#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <allegro5/allegro.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_primitives.h>




#define WINDOW_WIDTH                640
#define WINDOW_HEIGHT               240
#define WINDOW_FPS                  60


#define VIEWPORT_BLOCKSIZE          (4)
#define VIEWPORT_WIDTH              (WINDOW_WIDTH  / VIEWPORT_BLOCKSIZE)
#define VIEWPORT_HEIGHT             (WINDOW_HEIGHT / VIEWPORT_BLOCKSIZE)


#define WIND_SPEED                  0.20
#define WIND_VISCOSITY              1.40




class v2d {

    public:

        v2d()
            : pX(0.0), pY(0.0) { } 

        v2d(double xy)
            : pX(xy), pY(xy) { }

        v2d(double x, double y)
            : pX(x), pY(y) { }




        const double& x() const { return pX; }
        const double& y() const { return pY; }

        double& x() { return pX; }
        double& y() { return pY; }




        double len() const {
            return sqrt(pX * pX + pY * pY);
        }

        double len2() const {
            return len() * len();
        }


        static double dot(const v2d& v1, const v2d& v2) {
            return (v1.x() * v2.x()) + (v1.y() * v2.y());
        }

        static double dot2(const v2d& v1, const v2d& v2) {
            return dot(v1, v2) * dot(v1, v2);
        }



    private:

        double pX;
        double pY;

};






constexpr double wZero = 4.0 / 9.0;
constexpr double wCard = 1.0 / 9.0;
constexpr double wDiag = 1.0 / 36.0;


const double W[] = {
    wZero,
    wCard,
    wCard,
    wCard,
    wCard,
    wDiag,
    wDiag,
    wDiag,
    wDiag
};

const v2d E[] = {
    {  0,  0 },
    {  1,  0 },
    {  0,  1 },
    { -1,  0 },
    {  0, -1 },
    {  1,  1 },
    { -1,  1 },
    { -1, -1 },
    {  1, -1 },
};




struct unit {


    union {
        
        struct {
            double n0;
            double nE;
            double nN;
            double nW;
            double nS;
            double nNE;
            double nNW;
            double nSW;
            double nSE;
        };

        double n[9];

    };


    bool barrier;

    v2d u;
    double rho;
    double curl;



    void zero() {

        for(auto i = 0; i < 9; i++)
            n[i] = 0.0;

        rho = 0;
        curl = 0;

    }

    void eq(const double w, const double rho) {

        this->rho = rho;    
   
        for(auto i = 0; i < 9; i++)
            n[i] += w * (rho * W[i] * (1 + 3 * v2d::dot(E[i], u) + 4.5 * v2d::dot2(E[i], u) - 1.5 * u.len2()) - n[i]); 
        
    }



    const double new_rho() const {

        return n[0] + n[1] + n[2]
             + n[3] + n[4] + n[5]
             + n[6] + n[7] + n[8];

    }


};




class WIND {

    public:


        WIND() 
            : pUnits(VIEWPORT_WIDTH, std::vector<unit>(VIEWPORT_HEIGHT)), pFlowSpeed(0.0), pFlowViscosity(WIND_VISCOSITY) {
                
                setDirection( 1.0, 0.0);
                reset();

            }



        void setBlock(uint16_t x, uint16_t y) {

            pUnits[x][y].barrier = true;
            pUnits[x][y].rho = 0;
            pUnits[x][y].u = v2d();

            for(auto j = 0; j < 9; j++)
                pUnits[x][y].n[0];

        }


        void setDirection(double x, double y) {

            pFlowSpeed.x() = WIND_SPEED * x;
            pFlowSpeed.y() = WIND_SPEED * y;

            reset();

        }


        void clear() {

            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                    pUnits[x][y].barrier = false;
                    pUnits[x][y].zero();

                }
            }


            reset();


        }



        void step() {

            
            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                    if(!pUnits[x][y].barrier) {

                        auto& i = pUnits[x][y];
                        auto rho = i.new_rho();


                        if(rho > 0.0) {
                            i.u.x() = ((i.nE + i.nNE + i.nSE - i.nW - i.nNW - i.nSW) / rho);
                            i.u.y() = ((i.nN + i.nNE + i.nNW - i.nS - i.nSE - i.nSW) / rho);
                        } else
                            i.u = v2d();
                        

                        i.eq(pFlowViscosity, rho);

                    }

                }
            }



            for(auto x = 0; x < VIEWPORT_WIDTH - 1; x++) {
                for(auto y = VIEWPORT_HEIGHT - 1; y > 0; y--) {

                    pUnits[x][y].nN  = pUnits[x][y - 1].nN;
                    pUnits[x][y].nNW = pUnits[x + 1][y - 1].nNW;

                }
            }

            for(auto x = VIEWPORT_WIDTH - 1; x > 0; x--) {
                for(auto y = VIEWPORT_HEIGHT - 1; y > 0; y--) {

                    pUnits[x][y].nE  = pUnits[x - 1][y].nE;
                    pUnits[x][y].nNE = pUnits[x - 1][y - 1].nNE;

                }
            }

            for(auto x = VIEWPORT_WIDTH - 1; x > 0; x--) {
                for(auto y = 0; y < VIEWPORT_HEIGHT - 1; y++) {

                    pUnits[x][y].nS  = pUnits[x][y + 1].nS;
                    pUnits[x][y].nSE = pUnits[x - 1][y + 1].nSE;

                }
            }

            for(auto x = 0; x < VIEWPORT_WIDTH - 1; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT - 1; y++) {

                    pUnits[x][y].nW  = pUnits[x + 1][y].nW;
                    pUnits[x][y].nSW = pUnits[x + 1][y + 1].nSW;

                }
            }


            for(auto y = 0; y < VIEWPORT_HEIGHT - 1; y++)
                pUnits[0][y].nS = pUnits[0][y + 1].nS;
            
            for(auto y = VIEWPORT_HEIGHT - 1; y > 0; y--)
                pUnits[VIEWPORT_WIDTH - 1][y].nN = pUnits[VIEWPORT_WIDTH - 1][y - 1].nN;



            for(auto y = 0; y < VIEWPORT_HEIGHT - 1; y++) {
             
                if(!pUnits[0][y].barrier) {

                    const auto& u = pFlowSpeed;

                    pUnits[0][y].nE  = W[1] * (1 + 3 * v2d::dot(E[1], u) + 4.5 * v2d::dot2(E[1], u) - 1.5 * u.len2());
                    pUnits[0][y].nNE = W[5] * (1 + 3 * v2d::dot(E[5], u) + 4.5 * v2d::dot2(E[5], u) - 1.5 * u.len2());
                    pUnits[0][y].nSE = W[8] * (1 + 3 * v2d::dot(E[8], u) + 4.5 * v2d::dot2(E[8], u) - 1.5 * u.len2());

                }

                if(!pUnits[VIEWPORT_WIDTH - 1][y].barrier) {

                    const auto& u = pFlowSpeed;

                    pUnits[VIEWPORT_WIDTH - 1][y].nW  = W[3] * (1 + 3 * v2d::dot(E[3], u) + 4.5 * v2d::dot2(E[3], u) - 1.5 * u.len2());
                    pUnits[VIEWPORT_WIDTH - 1][y].nNW = W[6] * (1 + 3 * v2d::dot(E[6], u) + 4.5 * v2d::dot2(E[6], u) - 1.5 * u.len2());
                    pUnits[VIEWPORT_WIDTH - 1][y].nSW = W[7] * (1 + 3 * v2d::dot(E[7], u) + 4.5 * v2d::dot2(E[7], u) - 1.5 * u.len2());

                }

            }



            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {

                pUnits[x][0].zero();
                pUnits[x][0].eq(1, 1);

                pUnits[x][VIEWPORT_HEIGHT - 1].zero();
                pUnits[x][VIEWPORT_HEIGHT - 1].eq(1, 1);

            }
            


            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                    if(pUnits[x][y].barrier) {

                            pUnits[x][y - 1].nS += pUnits[x][y].nN;
                                                   pUnits[x][y].nN = 0;

                            pUnits[x][y + 1].nN += pUnits[x][y].nS;
                                                   pUnits[x][y].nS = 0;

                            pUnits[x - 1][y].nW += pUnits[x][y].nE;
                                                   pUnits[x][y].nE = 0;

                            pUnits[x + 1][y].nE += pUnits[x][y].nW;
                                                   pUnits[x][y].nW = 0;



                            pUnits[x + 1][y - 1].nSE += pUnits[x][y].nNW;
                                                        pUnits[x][y].nNW = 0;

                            pUnits[x - 1][y - 1].nSW += pUnits[x][y].nNE;
                                                        pUnits[x][y].nNE = 0;

                            pUnits[x + 1][y + 1].nNE += pUnits[x][y].nSW;
                                                        pUnits[x][y].nSW = 0;

                            pUnits[x - 1][y + 1].nNW += pUnits[x][y].nSE;
                                                        pUnits[x][y].nSE = 0;

            

                    }

                }
            }



            for(auto x = 1; x < VIEWPORT_WIDTH - 1; x++) {
                for(auto y = 1; y < VIEWPORT_HEIGHT - 1; y++) {

                    pUnits[x][y].curl = (pUnits[x + 1][y].u.y() - pUnits[x - 1][y].u.y()) 
                                      - (pUnits[x][y + 1].u.x() - pUnits[x][y - 1].u.x());

                }
            }

            for(auto y = 1; y < VIEWPORT_HEIGHT - 1; y++) {

                pUnits[0][y].curl = 2 * (pUnits[1][y].u.y()     - pUnits[0][y].u.y())
                                      - (pUnits[0][y - 1].u.x() - pUnits[0][y + 1].u.x());


                pUnits[VIEWPORT_WIDTH - 1][y].curl = 2 * (pUnits[VIEWPORT_WIDTH - 1][y].u.y()     - pUnits[VIEWPORT_WIDTH - 2][y].u.y())
                                                       - (pUnits[VIEWPORT_WIDTH - 1][y - 1].u.x() - pUnits[VIEWPORT_WIDTH - 1][y + 1].u.x());

            }

        }






        void reset() {

            
            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                    auto& i = pUnits[x][y];

                    i.zero();


                    if(i.barrier) {

                        i.rho = 0;
                        i.u = v2d();

                    } else {

                        i.rho = 1;
                        i.u = pFlowSpeed;
                        i.eq(1.0, 1.0);

                    }
                    
                }
            }



        }



        const auto& units() const {
            return pUnits;
        }

        const auto& flowSpeed() const {
            return pFlowSpeed;
        }

        const auto& flowViscosity() const {
            return pFlowViscosity;
        }

    

    private:

        std::vector<std::vector<unit>> pUnits;
        v2d pFlowSpeed;
        double pFlowViscosity;


};









class Application {

    public:

        Application(ALLEGRO_FONT* font) 
            : font(font), running(true) {

                
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 0);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 1);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 2);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 3);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 4);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 5);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 6);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 7);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 8);
                // wind.setBlock(VIEWPORT_WIDTH / 2, VIEWPORT_HEIGHT / 2 - 5 + 9);

                mode = 0;

                wind.reset();
 
                
            }
        


        void redraw(ALLEGRO_EVENT* e) {

            al_clear_to_color(al_map_rgb(0, 0, 0));


            for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                    auto& i = wind.units()[x][y];

                    const double dx = x * VIEWPORT_BLOCKSIZE;
                    const double dy = y * VIEWPORT_BLOCKSIZE;
                    const double dw = x * VIEWPORT_BLOCKSIZE + VIEWPORT_BLOCKSIZE;
                    const double dh = y * VIEWPORT_BLOCKSIZE + VIEWPORT_BLOCKSIZE;


                    if(i.barrier) {

                        al_draw_filled_rectangle(dx, dy, dw, dh, al_map_rgb(0, 0, 0));


                    } else {


                        double value = 0.0;

                        switch(mode) {

                            case 0:
                                value = i.curl;
                                break;

                            case 1:
                                value = i.new_rho() / 16;
                                break;

                            case 2:
                                value = i.u.len();
                                break;

                            case 3:
                                value = i.u.x();
                                break;

                            case 4:
                                value = i.u.y();
                                break;

                        }
                    

                        int c = 1024 * (value * 4 + 0.5);

                        if(c > 1024)
                            c = 1024;

                        if(c < 0)
                            c = 0;


                        int out[3];
                        HSVtoRGB((c / 1024.0) * (2.0 / 3.0) * 360.0 + 90, 0.75, 1, out);



                        char r, g, b;

                        r = out[0];
                        g = out[1];
                        b = out[2];
                    
                        al_draw_filled_rectangle(dx, dy, dw, dh, al_map_rgb(r, g, b));


                    }
                       

                }
            }

        }


        void step() {
            wind.step();
        }


        void update(ALLEGRO_EVENT* e) {


            switch(e->type) {

                case ALLEGRO_EVENT_KEY_DOWN:

                    
                    if(e->keyboard.keycode == ALLEGRO_KEY_ESCAPE)
                        running = false;

                    else if(e->keyboard.keycode == ALLEGRO_KEY_C)
                        wind.clear();


                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_7)
                        wind.setDirection(-0.5, -0.5);
                    
                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_8)
                        wind.setDirection(   0, -1.0);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_9)
                        wind.setDirection( 0.5, -0.5);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_4)
                        wind.setDirection(-1.0,  0.0);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_5)
                        wind.setDirection( 0.0,  0.0);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_6)
                        wind.setDirection( 1.0,  0.0);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_1)
                        wind.setDirection(-0.5,  0.5);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_2)
                        wind.setDirection( 0.0,  1.0);

                    else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_3)
                        wind.setDirection( 0.5,  0.5);


                    else if(e->keyboard.keycode == ALLEGRO_KEY_1)
                        mode = 0;

                    else if(e->keyboard.keycode == ALLEGRO_KEY_2)
                        mode = 1;

                    else if(e->keyboard.keycode == ALLEGRO_KEY_3)
                        mode = 2;

                    else if(e->keyboard.keycode == ALLEGRO_KEY_4)
                        mode = 3;

                    else if(e->keyboard.keycode == ALLEGRO_KEY_5)
                        mode = 4;

                    break;


                case ALLEGRO_EVENT_MOUSE_AXES:

                    if(e->mouse.pressure > 0) {

                        wind.setBlock(e->mouse.x / VIEWPORT_BLOCKSIZE, e->mouse.y / VIEWPORT_BLOCKSIZE);
                        wind.reset();

                    }

                    break;

            }

        

        }

        void exit() {

        }

        inline bool isRunning() { return running; }



    private:

        ALLEGRO_FONT* font;
        WIND wind;

        bool running;
        uint8_t mode;




        void HSVtoRGB(int H, double S, double V, int output[3]) {


            double C = S * V;
            double X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
            double m = V - C;
            double Rs, Gs, Bs;


            if(H >= 0 && H < 60) {
                Rs = C;
                Gs = X;
                Bs = 0;	
            }

            else if(H >= 60 && H < 120) {	
                Rs = X;
                Gs = C;
                Bs = 0;	
            }

            else if(H >= 120 && H < 180) {
                Rs = 0;
                Gs = C;
                Bs = X;	
            }

            else if(H >= 180 && H < 240) {
                Rs = 0;
                Gs = X;
                Bs = C;	
            }

            else if(H >= 240 && H < 300) {
                Rs = X;
                Gs = 0;
                Bs = C;	
            }

            else {
                Rs = C;
                Gs = 0;
                Bs = X;	
            }
            
            output[0] = (Rs + m) * 255;
            output[1] = (Gs + m) * 255;
            output[2] = (Bs + m) * 255;


        }



};




int main(int argc, char** argv) {




    al_init();
    al_install_keyboard();
    al_install_mouse();

    al_set_app_name("Game of Life");


    ALLEGRO_EVENT_QUEUE* queue;
    if((queue = al_create_event_queue()) == NULL)
        return std::cerr << "al_create_event_queue() failed!" << std::endl, 1;

    ALLEGRO_DISPLAY* disp;
    if((disp = al_create_display(WINDOW_WIDTH, WINDOW_HEIGHT)) == NULL)
        return std::cerr << "al_create_display() failed!" << std::endl, 1;

    ALLEGRO_TIMER* timer;
    if((timer = al_create_timer(1.0 / WINDOW_FPS)) == NULL)
        return std::cerr << "al_create_timer() failed!" << std::endl, 1;

    ALLEGRO_FONT* font;
    if((font = al_create_builtin_font()) == NULL)
        return std::cerr << "al_create_builtin_font() failed!" << std::endl, 1;



    al_set_window_title(disp, "WIND - Parallel Algorithm and Data Structures - Exam");


    al_register_event_source(queue, al_get_keyboard_event_source());
    al_register_event_source(queue, al_get_mouse_event_source());
    al_register_event_source(queue, al_get_display_event_source(disp));
    al_register_event_source(queue, al_get_timer_event_source(timer));

    al_start_timer(timer);





    Application app(font);

    do {



        ALLEGRO_EVENT e;
        if(al_get_next_event(queue, &e)) {

            
            switch(e.type) {

                case ALLEGRO_EVENT_TIMER:
                    app.redraw(&e);
                    al_flip_display();
                    break;

                case ALLEGRO_EVENT_KEY_DOWN:
                case ALLEGRO_EVENT_KEY_UP:
                case ALLEGRO_EVENT_MOUSE_BUTTON_DOWN:
                case ALLEGRO_EVENT_MOUSE_BUTTON_UP:
                case ALLEGRO_EVENT_MOUSE_AXES:
                case ALLEGRO_EVENT_MOUSE_ENTER_DISPLAY:
                case ALLEGRO_EVENT_MOUSE_LEAVE_DISPLAY:
                    app.update(&e);
                    break;

                case ALLEGRO_EVENT_DISPLAY_CLOSE:
                    app.exit();
                    break;

            }

        }



        app.step();



#if !defined(BENCH)
        usleep(1000);
#endif    

    } while(app.isRunning());


    return 0;

}
