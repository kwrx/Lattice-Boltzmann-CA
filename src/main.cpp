 
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
#include <sstream>
#include <cmath>
#include <cassert>

#include <allegro5/allegro.h>
#include <allegro5/allegro_font.h>
#include <allegro5/allegro_primitives.h>

#include <mpi.h>




#define PRIMARY                     0


#define WINDOW_WIDTH                640
#define WINDOW_HEIGHT               240
#define WINDOW_FPS                  30


#define VIEWPORT_BLOCKSIZE          (4)
#define VIEWPORT_WIDTH              (WINDOW_WIDTH  / VIEWPORT_BLOCKSIZE)
#define VIEWPORT_HEIGHT             (WINDOW_HEIGHT / VIEWPORT_BLOCKSIZE)


#define WIND_SPEED                  0.20
#define WIND_VISCOSITY              1.40


#define XY(x, y, w)     \
    (((y) * (w)) + (x))





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

        rho = 0.0;
        curl = 0.0;

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







static ALLEGRO_EVENT_QUEUE* queue;
static ALLEGRO_DISPLAY* disp;
static ALLEGRO_TIMER* timer;
static ALLEGRO_FONT* font;


static MPI_Comm MPI_COMM_LOCAL;
static MPI_Win MPI_LOCAL_WINDOW;
static MPI_Datatype MPI_TYPE_V2D;
static MPI_Datatype MPI_TYPE_UNIT;


static int world_rank;
static int world_num_procs;

static int local_rank;
static int local_num_procs;


static unit* up_units       = nullptr;
static unit* bottom_units   = nullptr;
static unit* units          = nullptr;
static unit* frame          = nullptr;
static unit* current_unit   = nullptr;

static uint16_t current_unit_x = 0;
static uint16_t current_unit_y = 0;




static v2d flow_speed = v2d(WIND_SPEED, 0.0);
static double flow_viscosity = WIND_VISCOSITY;

static bool running = true;
static bool resetting = true;
static bool paused = false;
static bool draw_groups = false;
static bool draw_nodes = false;
static uint8_t draw_mode = 0;







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



void reset() {
    resetting = true;
}


void clear() {

    for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
        for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

            frame[XY(x, y, VIEWPORT_WIDTH)].barrier = false;
            frame[XY(x, y, VIEWPORT_WIDTH)].zero();

        }
    }

    reset();

}


void setDirection(double x, double y) {

    flow_speed.x() = WIND_SPEED * x;
    flow_speed.y() = WIND_SPEED * y;

    reset();

}


void setBarrier(int x, int y) {

    if(x <= 0 || y <= 0 || x >= VIEWPORT_WIDTH - 1 || y >= VIEWPORT_HEIGHT - 1)
        return;

    if(frame[XY(x, y, VIEWPORT_WIDTH)].barrier)
        return;


    frame[XY(x, y, VIEWPORT_WIDTH)].barrier = true;
    frame[XY(x, y, VIEWPORT_WIDTH)].rho = 0.0;
    frame[XY(x, y, VIEWPORT_WIDTH)].u = v2d();
    frame[XY(x, y, VIEWPORT_WIDTH)].zero();

    reset();

}


void setCurrentUnit(int x, int y) {

    if(x < 0 || y < 0 || x > VIEWPORT_WIDTH - 1 || y > VIEWPORT_HEIGHT - 1)
        return;

    if(frame[XY(x, y, VIEWPORT_WIDTH)].barrier)
        return;

    current_unit = &frame[XY(x, y, VIEWPORT_WIDTH)];
    current_unit_x = x;
    current_unit_x = y;

}



void redraw(ALLEGRO_EVENT* e) {


    al_clear_to_color(al_map_rgb(0, 0, 0));


    for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
        for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {


            auto& i = frame[XY(x, y, VIEWPORT_WIDTH)];

            const double dx = x * VIEWPORT_BLOCKSIZE;
            const double dy = y * VIEWPORT_BLOCKSIZE;
            const double dw = x * VIEWPORT_BLOCKSIZE + VIEWPORT_BLOCKSIZE;
            const double dh = y * VIEWPORT_BLOCKSIZE + VIEWPORT_BLOCKSIZE;


            if(i.barrier) {

                al_draw_filled_rectangle(dx, dy, dw, dh, al_map_rgb(0, 0, 0));


            } else {


                double value = 0.0;

                switch(draw_mode) {

                    case 0: {
                        
                            value = i.curl;


                            constexpr int steps = 8;
                            constexpr int averg = 3;

                            if(y > steps && y < VIEWPORT_HEIGHT - steps) {

                                if(((y + averg) % (VIEWPORT_HEIGHT / world_num_procs)) < (averg << 1)) {

                                    for(auto n = 0; n < (steps >> 1); n++)
                                        value += frame[XY(x, y - n, VIEWPORT_WIDTH)].curl;

                                    for(auto n = 0; n < (steps >> 1); n++)
                                        value += frame[XY(x, y + n, VIEWPORT_WIDTH)].curl;

                                    value /= steps;

                                }


                            }


                        } break;

                    case 1:
                        value = i.new_rho() / 16.0;
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



    if(draw_groups) {


        static double group_colors[8][3] = {
            { 1.0, 0.2, 0.2 },
            { 0.2, 1.0, 0.2 },
            { 0.2, 0.2, 1.0 },
            { 0.5, 0.5, 0.5 },
            { 0.8, 0.2, 0.6 },
            { 0.6, 0.8, 0.2 },
            { 0.2, 0.8, 0.6 },
            { 0.8, 0.8, 0.6 },
        };


        const double n = world_num_procs / local_num_procs;
        const double h = VIEWPORT_HEIGHT / (world_num_procs / local_num_procs);

        for(auto i = 0; i < n; i++) {

            const double dx = 10;
            const double dw = WINDOW_WIDTH - 10;
            const double dy = (i * h) * VIEWPORT_BLOCKSIZE;
            const double dh = (dy + (h * VIEWPORT_BLOCKSIZE));

            std::stringstream ss;
            ss << i;

            al_draw_filled_rounded_rectangle(dx, dy, dw, dh, 10, 10, al_map_rgba_f(group_colors[i][0], group_colors[i][1], group_colors[i][2], 0.25));
            al_draw_text(font, al_map_rgb(25, 25, 25), dx + 10, dy + 10, 0, ss.str().c_str());

        }

    }


    if(draw_nodes) {


        static double nodes_colors[8][3] = {
            { 1.0, 0.2, 0.2 },
            { 0.2, 1.0, 0.2 },
            { 0.2, 0.2, 1.0 },
            { 0.5, 0.5, 0.5 },
            { 0.8, 0.2, 0.6 },
            { 0.6, 0.8, 0.2 },
            { 0.2, 0.8, 0.6 },
            { 0.8, 0.8, 0.6 },
        };


        const double n = world_num_procs;
        const double h = VIEWPORT_HEIGHT / world_num_procs;

        for(auto i = 0; i < n; i++) {

            const double dx = 40;
            const double dw = WINDOW_WIDTH - 20;
            const double dy = (i * h) * VIEWPORT_BLOCKSIZE + 5;
            const double dh = (dy + (h * VIEWPORT_BLOCKSIZE)) - 10;

            
            std::stringstream ss1, ss2;
            ss1 << i % world_num_procs;
            ss2 << i % local_num_procs;

            al_draw_filled_rounded_rectangle(dx, dy, dw, dh, 10, 10, al_map_rgba_f(nodes_colors[i][0], nodes_colors[i][1], nodes_colors[i][2], 0.75));

            al_draw_text(font, al_map_rgb(25, 25, 25), dx + 10, dy + 10, 0, ss1.str().c_str());
            al_draw_text(font, al_map_rgb(25, 25, 25), dw - 20, dy + 10, 0, ss2.str().c_str());

        }

    }


    if(current_unit) {

        std::stringstream ss;
        ss << "Unit(" << current_unit_x << ", " << current_unit_y << ") "
           << "Density: " << current_unit->new_rho() << ", Curl: " << current_unit->curl << ", "
           << "Speed: X(" << current_unit->u.x() << "," << current_unit->u.y() << ") " << current_unit->u.len();

        al_draw_text(font, al_map_rgb(25, 25, 25), 10, WINDOW_HEIGHT - 15, 0, ss.str().c_str());

    }




}





void update(ALLEGRO_EVENT* e) {


    switch(e->type) {

        case ALLEGRO_EVENT_KEY_DOWN:

            
            if(e->keyboard.keycode == ALLEGRO_KEY_ESCAPE) 
                running = false;

            else if(e->keyboard.keycode == ALLEGRO_KEY_C)
                clear();

            else if(e->keyboard.keycode == ALLEGRO_KEY_P)
                paused = !paused;

            else if(e->keyboard.keycode == ALLEGRO_KEY_I)
                draw_groups = !draw_groups;
            
            else if(e->keyboard.keycode == ALLEGRO_KEY_U)
                draw_nodes = !draw_nodes;


            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_7)
                setDirection(-0.5, -0.5);
            
            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_8)
                setDirection(   0, -1.0);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_9)
                setDirection( 0.5, -0.5);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_4)
                setDirection(-1.0,  0.0);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_5)
                setDirection( 0.0,  0.0);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_6)
                setDirection( 1.0,  0.0);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_1)
                setDirection(-0.5,  0.5);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_2)
                setDirection( 0.0,  1.0);

            else if(e->keyboard.keycode == ALLEGRO_KEY_PAD_3)
                setDirection( 0.5,  0.5);


            else if(e->keyboard.keycode == ALLEGRO_KEY_1)
                draw_mode = 0;

            else if(e->keyboard.keycode == ALLEGRO_KEY_2)
                draw_mode = 1;

            else if(e->keyboard.keycode == ALLEGRO_KEY_3)
                draw_mode = 2;

            else if(e->keyboard.keycode == ALLEGRO_KEY_4)
                draw_mode = 3;

            else if(e->keyboard.keycode == ALLEGRO_KEY_5)
                draw_mode = 4;

            break;


        case ALLEGRO_EVENT_MOUSE_BUTTON_DOWN:
        case ALLEGRO_EVENT_MOUSE_AXES:

            if(e->mouse.pressure > 0) {

                setBarrier(e->mouse.x / VIEWPORT_BLOCKSIZE, e->mouse.y / VIEWPORT_BLOCKSIZE);

            }


            setCurrentUnit(e->mouse.x / VIEWPORT_BLOCKSIZE, e->mouse.y / VIEWPORT_BLOCKSIZE);

            break;

    }



}





int main(int argc, char** argv) {


    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_num_procs);



    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, world_rank, MPI_INFO_NULL, &MPI_COMM_LOCAL);

    /**
     * Divisione in gruppi a scopo di test
     *
     * MPI_Comm_split(MPI_COMM_WORLD, world_rank, world_rank, &MPI_COMM_LOCAL);                             // 1 nodo per ogni gruppo
     * MPI_Comm_split(MPI_COMM_WORLD, world_rank < (world_num_procs / 2), world_rank, &MPI_COMM_LOCAL);     // n/2 nodi per ogni gruppo
     */




    MPI_Comm_rank(MPI_COMM_LOCAL, &local_rank);
    MPI_Comm_size(MPI_COMM_LOCAL, &local_num_procs);


    assert(world_num_procs > 0);
    assert(local_num_procs > 0);




#if !defined(BENCH)

    std::cout << "Running Node " << world_rank << " of " << world_num_procs 
              << " (" << local_rank << " of " << local_num_procs << ")" << std::endl;

#endif








    MPI_Datatype v2d_types[] = { MPI_DOUBLE, MPI_DOUBLE };
    MPI_Aint v2d_offsets[]   = { 0, sizeof(double) };
    int v2d_blocks[]         = { 1, 1 };

    MPI_Type_create_struct(2, v2d_blocks, v2d_offsets, v2d_types, &MPI_TYPE_V2D);
    MPI_Type_commit(&MPI_TYPE_V2D);





    MPI_Datatype unit_types[] =
        { MPI_DOUBLE, MPI_CXX_BOOL, MPI_TYPE_V2D, MPI_DOUBLE, MPI_DOUBLE };
        
    MPI_Aint unit_offsets[] = {
        offsetof(unit, n),
        offsetof(unit, barrier),
        offsetof(unit, u),
        offsetof(unit, rho),
        offsetof(unit, curl),
    };

    int unit_blocks[] = { 9, 1, 1, 1, 1 };


    MPI_Type_create_struct(5, unit_blocks, unit_offsets, unit_types, &MPI_TYPE_UNIT);
    MPI_Type_commit(&MPI_TYPE_UNIT);





#if !defined(BENCH)

    if(world_rank == PRIMARY) {


        al_init();
        al_install_keyboard();
        al_install_mouse();

        al_set_app_name("WIND - APSD");

    
        if((queue = al_create_event_queue()) == NULL)
            return std::cerr << "al_create_event_queue() failed!" << std::endl, 1;

        if((disp = al_create_display(WINDOW_WIDTH, WINDOW_HEIGHT)) == NULL)
            return std::cerr << "al_create_display() failed!" << std::endl, 1;

        if((timer = al_create_timer(1.0 / WINDOW_FPS)) == NULL)
            return std::cerr << "al_create_timer() failed!" << std::endl, 1;

        if((font = al_create_builtin_font()) == NULL)
            return std::cerr << "al_create_builtin_font() failed!" << std::endl, 1;



        al_set_window_title(disp, "WIND - Parallel Algorithm and Data Structures - Exam");


        al_register_event_source(queue, al_get_keyboard_event_source());
        al_register_event_source(queue, al_get_mouse_event_source());
        al_register_event_source(queue, al_get_display_event_source(disp));
        al_register_event_source(queue, al_get_timer_event_source(timer));

        al_start_timer(timer);

    }


#endif



    const size_t unit_width  = VIEWPORT_WIDTH;
    const size_t unit_height = VIEWPORT_HEIGHT / world_num_procs;
    const size_t unit_size   = unit_width * unit_height;

    if(MPI_Win_allocate_shared (unit_size * sizeof(struct unit), sizeof(struct unit), MPI_INFO_NULL, MPI_COMM_LOCAL, &units, &MPI_LOCAL_WINDOW) != MPI_SUCCESS)
        MPI_Abort(MPI_COMM_WORLD, __LINE__);



    if(world_rank == PRIMARY) {

        frame = new unit[VIEWPORT_WIDTH * VIEWPORT_HEIGHT];

        if(!frame)
            MPI_Abort(MPI_COMM_WORLD, __LINE__);


    }




    up_units = new unit[unit_width];
    bottom_units = new unit[unit_width];


    if(!up_units || !bottom_units)
        MPI_Abort(MPI_COMM_WORLD, __LINE__);


    




#if defined(BENCH)

    double bench_start = MPI_Wtime();

#endif


    do {


#if !defined(BENCH)

        if(world_rank == PRIMARY) {

            ALLEGRO_EVENT e;
            if(al_get_next_event(queue, &e)) {

                
                switch(e.type) {

                    case ALLEGRO_EVENT_TIMER:
                        redraw(&e);
                        al_flip_display();
                        break;

                    case ALLEGRO_EVENT_KEY_DOWN:
                    case ALLEGRO_EVENT_KEY_UP:
                    case ALLEGRO_EVENT_MOUSE_BUTTON_DOWN:
                    case ALLEGRO_EVENT_MOUSE_BUTTON_UP:
                    case ALLEGRO_EVENT_MOUSE_AXES:
                    case ALLEGRO_EVENT_MOUSE_ENTER_DISPLAY:
                    case ALLEGRO_EVENT_MOUSE_LEAVE_DISPLAY:
                        update(&e);
                        break;

                    case ALLEGRO_EVENT_DISPLAY_CLOSE:
                        running = false;
                        break;

                }

            }

        }

#endif



        MPI_Bcast(&resetting,       1, MPI_CXX_BOOL, PRIMARY, MPI_COMM_WORLD);
        MPI_Bcast(&running,         1, MPI_CXX_BOOL, PRIMARY, MPI_COMM_WORLD);
        MPI_Bcast(&paused,          1, MPI_CXX_BOOL, PRIMARY, MPI_COMM_WORLD);
        MPI_Bcast(&flow_viscosity,  1, MPI_DOUBLE,   PRIMARY, MPI_COMM_WORLD);
        MPI_Bcast(&flow_speed,      1, MPI_TYPE_V2D, PRIMARY, MPI_COMM_WORLD);
        MPI_Bcast(&draw_mode,       1, MPI_INT,      PRIMARY, MPI_COMM_WORLD);



#if defined(BENCH)
        
        static uint32_t iterations = 0;

        if(++iterations == ITERATIONS)
            running = false;

#endif




        MPI_Barrier(MPI_COMM_WORLD);

        

        if(__sync_bool_compare_and_swap(&resetting, true, false)) {


            if(world_rank != PRIMARY) {
                MPI_Recv(units, unit_size, MPI_TYPE_UNIT, PRIMARY, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


            } else {


                
                for(auto x = 0; x < VIEWPORT_WIDTH; x++) {
                    for(auto y = 0; y < VIEWPORT_HEIGHT; y++) {

                        auto& i = frame[XY(x, y, VIEWPORT_WIDTH)];

                        i.zero();


                        if(i.barrier) {

                            i.rho = 0.0;
                            i.u = v2d();

                        } else {

                            i.rho = 1.0;
                            i.u = flow_speed;
                            i.eq(1.0, 1.0);

                        }

                    }
                }


                for(auto i = PRIMARY; i < world_num_procs; i++) {


                    if(i == PRIMARY)
                        memcpy(units, &frame[XY(0, i * unit_height, VIEWPORT_WIDTH)], unit_size * sizeof(unit));

                    else
                    
                        MPI_Send(&frame[XY(0, i * unit_height, VIEWPORT_WIDTH)], unit_size, MPI_TYPE_UNIT, i, 0, MPI_COMM_WORLD);
                    

                }





            }

        }



        if(paused)
            continue;



        for(auto x = 0; x < unit_width; x++) {
            for(auto y = 0; y < unit_height; y++) {

                if(!units[XY(x, y, unit_width)].barrier) {

                    auto& i = units[XY(x, y, unit_width)];
                    auto rho = i.new_rho();


                    if(rho > 0.0) {

                        i.u.x() = ((i.nE + i.nNE + i.nSE - i.nW - i.nNW - i.nSW) / rho);
                        i.u.y() = ((i.nN + i.nNE + i.nNW - i.nS - i.nSE - i.nSW) / rho);

                    } else

                        i.u = v2d();
                    

                    i.eq(flow_viscosity, rho);


                }

            }
        }


        for(auto y = 0; y < unit_height; y++) {
    
            const auto& u = flow_speed;

            units[XY(0, y, unit_width)].nE  = W[1] * (1 + 3 * v2d::dot(E[1], u) + 4.5 * v2d::dot2(E[1], u) - 1.5 * u.len2());
            units[XY(0, y, unit_width)].nNE = W[5] * (1 + 3 * v2d::dot(E[5], u) + 4.5 * v2d::dot2(E[5], u) - 1.5 * u.len2());
            units[XY(0, y, unit_width)].nSE = W[8] * (1 + 3 * v2d::dot(E[8], u) + 4.5 * v2d::dot2(E[8], u) - 1.5 * u.len2());


            units[XY(unit_width - 1, y, unit_width)].nW  = W[3] * (1 + 3 * v2d::dot(E[3], u) + 4.5 * v2d::dot2(E[3], u) - 1.5 * u.len2());
            units[XY(unit_width - 1, y, unit_width)].nNW = W[6] * (1 + 3 * v2d::dot(E[6], u) + 4.5 * v2d::dot2(E[6], u) - 1.5 * u.len2());
            units[XY(unit_width - 1, y, unit_width)].nSW = W[7] * (1 + 3 * v2d::dot(E[7], u) + 4.5 * v2d::dot2(E[7], u) - 1.5 * u.len2());


        }
            








        
        MPI_Win_fence(0, MPI_LOCAL_WINDOW);


        #define LOCAL_WIDTH     (unit_width)
        #define LOCAL_HEIGHT    (unit_height)




        for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {
            for(auto y = LOCAL_HEIGHT - 1; y > 0; y--) {

                units[XY(x, y, LOCAL_WIDTH)].nN  = units[XY(x + 0, y - 1, LOCAL_WIDTH)].nN;
                units[XY(x, y, LOCAL_WIDTH)].nNW = units[XY(x + 1, y - 1, LOCAL_WIDTH)].nNW;

            }
        }

        for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {
            for(auto y = LOCAL_HEIGHT - 1; y > 0; y--) {

                units[XY(x, y, LOCAL_WIDTH)].nE  = units[XY(x - 1, y, LOCAL_WIDTH)].nE;
                units[XY(x, y, LOCAL_WIDTH)].nNE = units[XY(x - 1, y - 1, LOCAL_WIDTH)].nNE;

            }
        }

        
        for(auto y = LOCAL_HEIGHT - 1; y > 0; y--)
            units[XY(LOCAL_WIDTH - 1, y, LOCAL_WIDTH)].nN = units[XY(LOCAL_WIDTH - 1, y - 1, LOCAL_WIDTH)].nN;




        if(world_num_procs > 1) {


            if(world_rank != PRIMARY) {

                if(local_rank == PRIMARY) {

                    MPI_Sendrecv (
                        &units[0],    unit_width, MPI_TYPE_UNIT, world_rank - 1, 0,
                        &up_units[0], unit_width, MPI_TYPE_UNIT, world_rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
                    );

                } else {

                    memcpy(&up_units[0], (void*) ((uintptr_t) &units[0] - (LOCAL_WIDTH * sizeof(unit))), LOCAL_WIDTH * sizeof(unit));

                }


                for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, 0, LOCAL_WIDTH)].nN  = up_units[XY(x + 0, 0, LOCAL_WIDTH)].nN;
                    units[XY(x, 0, LOCAL_WIDTH)].nNW = up_units[XY(x + 1, 0, LOCAL_WIDTH)].nNW;
                    
                }

                for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {

                    units[XY(x, 0, LOCAL_WIDTH)].nNE = up_units[XY(x - 1, 0, LOCAL_WIDTH)].nNE;
                    
                }


                units[XY(LOCAL_WIDTH - 1, 0, LOCAL_WIDTH)].nN = up_units[XY(LOCAL_WIDTH - 1, 0, LOCAL_WIDTH)].nN;




                for(auto x = 1; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, 0, LOCAL_WIDTH)].curl = (units[XY(x + 1, 0, LOCAL_WIDTH)].u.y() - units[XY(x - 1, 0, LOCAL_WIDTH)].u.y())
                                                      - (units[XY(x, 1, LOCAL_WIDTH)].u.x() - up_units[XY(x, 0, LOCAL_WIDTH)].u.x());
                    
                }



                for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, 0, LOCAL_WIDTH)].nW  = units[XY(x + 1, 0, LOCAL_WIDTH)].nW;
                    units[XY(x, 0, LOCAL_WIDTH)].nSW = units[XY(x + 1, 1, LOCAL_WIDTH)].nSW;


                }


                for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {

                    units[XY(x, 0, LOCAL_WIDTH)].nS  = units[XY(x, 1, LOCAL_WIDTH)].nS;
                    units[XY(x, 0, LOCAL_WIDTH)].nSE = units[XY(x - 1, 1, LOCAL_WIDTH)].nSE;
                    units[XY(x, 0, LOCAL_WIDTH)].nE  = units[XY(x - 1, 0, LOCAL_WIDTH)].nE;


                }


                units[XY(0, 0, LOCAL_WIDTH)].nS = units[XY(0, 1, LOCAL_WIDTH)].nS;


            }


        } else {

            for(auto x = 0; x < LOCAL_WIDTH; x++) {

                units[XY(x, 0, LOCAL_WIDTH)].zero();
                units[XY(x, 0, LOCAL_WIDTH)].eq(1, 1);

            }

        }

        







        for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {
            for(auto y = 0; y < LOCAL_HEIGHT - 1; y++) {

                units[XY(x, y, LOCAL_WIDTH)].nS  = units[XY(x, y + 1, LOCAL_WIDTH)].nS;
                units[XY(x, y, LOCAL_WIDTH)].nSE = units[XY(x - 1, y + 1, LOCAL_WIDTH)].nSE;

            }
        }

        for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {
            for(auto y = 0; y < LOCAL_HEIGHT - 1; y++) {

                units[XY(x, y, LOCAL_WIDTH)].nW  = units[XY(x + 1, y, LOCAL_WIDTH)].nW;
                units[XY(x, y, LOCAL_WIDTH)].nSW = units[XY(x + 1, y + 1, LOCAL_WIDTH)].nSW;

            }
        }   



        for(auto y = 0; y < LOCAL_HEIGHT - 1; y++)
            units[XY(0, y, LOCAL_WIDTH)].nS = units[XY(0, y + 1, LOCAL_WIDTH)].nS;






        if(world_num_procs > 1) {


            if(world_rank != (world_num_procs - 1)) {

                if(local_rank == local_num_procs - 1) {

                    MPI_Sendrecv (
                        &units[XY(0, LOCAL_HEIGHT - 1, LOCAL_WIDTH)], unit_width, MPI_TYPE_UNIT, world_rank + 1, 0,
                        &bottom_units[0],                             unit_width, MPI_TYPE_UNIT, world_rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE
                    );

                } else {

                    memcpy(&bottom_units[0], &units[XY(0, LOCAL_HEIGHT, LOCAL_WIDTH)], LOCAL_WIDTH * sizeof(unit));

                }


                for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nSW = bottom_units[XY(x + 1, 0, LOCAL_WIDTH)].nSW;
                    
                }

                for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {

                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nS  = bottom_units[XY(x + 0, 0, LOCAL_WIDTH)].nS;
                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nSE = bottom_units[XY(x - 1, 0, LOCAL_WIDTH)].nSE;
                
                }


                units[XY(0, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nS = bottom_units[XY(0, 0, LOCAL_WIDTH)].nS;




                for(int x = 1; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].curl = (units[XY(x + 1, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].u.y() - units[XY(x - 1, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].u.y())
                                                                     - (bottom_units[XY(x, 0, LOCAL_WIDTH)].u.x() - units[XY(x, LOCAL_HEIGHT - 2, LOCAL_WIDTH)].u.x());
                    
                }




                for(auto x = 0; x < LOCAL_WIDTH - 1; x++) {

                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nN  = units[XY(x + 0, LOCAL_HEIGHT - 2, LOCAL_WIDTH)].nN;
                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nNW = units[XY(x + 1, LOCAL_HEIGHT - 2, LOCAL_WIDTH)].nNW;
                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nW  = units[XY(x + 1, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nW;


                }


                for(auto x = LOCAL_WIDTH - 1; x > 0; x--) {

                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nE  = units[XY(x - 1, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nE;
                    units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nNE = units[XY(x - 1, LOCAL_HEIGHT - 2, LOCAL_WIDTH)].nNE;

                }


                units[XY(LOCAL_WIDTH - 1, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].nN = units[XY(LOCAL_WIDTH - 1, LOCAL_HEIGHT - 2, LOCAL_WIDTH)].nN;


            }
            

        } else {

            for(auto x = 0; x < LOCAL_WIDTH; x++) {

                units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].zero();
                units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].eq(1, 1);

            }

        }







        if(world_rank == (world_num_procs - 1)) {

            for(auto x = 0; x < LOCAL_WIDTH; x++) {

                units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].zero();
                units[XY(x, LOCAL_HEIGHT - 1, LOCAL_WIDTH)].eq(1, 1);

            }

        }

        if(world_rank == PRIMARY) {

            for(auto x = 0; x < LOCAL_WIDTH; x++) {

                units[XY(x, 0, LOCAL_WIDTH)].zero();
                units[XY(x, 0, LOCAL_WIDTH)].eq(1, 1);

            }


        } 






        for(auto x = 1; x < LOCAL_WIDTH - 1; x++) {
            for(auto y = 1; y < LOCAL_HEIGHT - 1; y++) {

                units[XY(x, y, LOCAL_WIDTH)].curl = (units[XY(x + 1, y, LOCAL_WIDTH)].u.y() - units[XY(x - 1, y, LOCAL_WIDTH)].u.y()) 
                                                    - (units[XY(x, y + 1, LOCAL_WIDTH)].u.x() - units[XY(x, y - 1, LOCAL_WIDTH)].u.x());

            }
        }

        for(auto y = 1; y < LOCAL_HEIGHT - 1; y++) {

            units[XY(0, y, LOCAL_WIDTH)].curl = (units[XY(1, y, LOCAL_WIDTH)].u.y()     - units[XY(0, y, LOCAL_WIDTH)].u.y())
                                                - (units[XY(0, y - 1, LOCAL_WIDTH)].u.x() - units[XY(0, y + 1, LOCAL_WIDTH)].u.x());


            units[XY(LOCAL_WIDTH - 1, y, LOCAL_WIDTH)].curl = (units[XY(LOCAL_WIDTH - 1, y, LOCAL_WIDTH)].u.y()     - units[XY(LOCAL_WIDTH - 2, y, LOCAL_WIDTH)].u.y())
                                                            - (units[XY(LOCAL_WIDTH - 1, y - 1, LOCAL_WIDTH)].u.x() - units[XY(LOCAL_WIDTH - 1, y + 1, LOCAL_WIDTH)].u.x());

        }






        for(auto x = 1; x < LOCAL_WIDTH - 1; x++) {
            for(auto y = 1; y < LOCAL_HEIGHT - 1; y++) {

                if(units[XY(x, y, LOCAL_WIDTH)].barrier) {


                        units[XY(x, y - 1, LOCAL_WIDTH)].nS += units[XY(x, y, LOCAL_WIDTH)].nN;
                                                                units[XY(x, y, LOCAL_WIDTH)].nN = 0;


                        units[XY(x, y + 1, LOCAL_WIDTH)].nN += units[XY(x, y, LOCAL_WIDTH)].nS;
                                                                units[XY(x, y, LOCAL_WIDTH)].nS = 0;


                        units[XY(x - 1, y, LOCAL_WIDTH)].nW += units[XY(x, y, LOCAL_WIDTH)].nE;
                                                                units[XY(x, y, LOCAL_WIDTH)].nE = 0;


                        units[XY(x + 1, y, LOCAL_WIDTH)].nE += units[XY(x, y, LOCAL_WIDTH)].nW;
                                                                units[XY(x, y, LOCAL_WIDTH)].nW = 0;



                        units[XY(x + 1, y - 1, LOCAL_WIDTH)].nSE += units[XY(x, y, LOCAL_WIDTH)].nNW;
                                                                    units[XY(x, y, LOCAL_WIDTH)].nNW = 0;


                        units[XY(x - 1, y - 1, LOCAL_WIDTH)].nSW += units[XY(x, y, LOCAL_WIDTH)].nNE;
                                                                    units[XY(x, y, LOCAL_WIDTH)].nNE = 0;


                        units[XY(x + 1, y + 1, LOCAL_WIDTH)].nNE += units[XY(x, y, LOCAL_WIDTH)].nSW;
                                                                    units[XY(x, y, LOCAL_WIDTH)].nSW = 0;


                        units[XY(x - 1, y + 1, LOCAL_WIDTH)].nNW += units[XY(x, y, LOCAL_WIDTH)].nSE;
                                                                    units[XY(x, y, LOCAL_WIDTH)].nSE = 0;

        

                }

            }
        }





        
        MPI_Gather(units, unit_size, MPI_TYPE_UNIT, frame, unit_size, MPI_TYPE_UNIT, PRIMARY, MPI_COMM_WORLD);


#if !defined(BENCH)
        usleep(1000);
#endif

    } while(running);



#if defined(BENCH)

    double bench_end = MPI_Wtime();

    if(world_rank == PRIMARY)
        std::cerr << std::fixed << (bench_end - bench_start) << std::endl; 

#endif


    return MPI_Finalize();

}