#include <easy3d/viewer/model.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/fileio/resources.h>
#include <easy3d/util/logging.h>
#include <easy3d/viewer/viewer.h>
#include <easy3d/viewer/camera.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/fileio/point_cloud_io.h>
#include <fstream>


using namespace easy3d;

class ImprovedViewer: public Viewer {
    std::vector<vec3> points;
public:
    vec3 obracun(int x, int y){
        float visina;
        float sirina;

        this->camera()->getOrthoWidthHeight(sirina, visina);
        int W = camera()->screenWidth();
        int H = camera()->screenHeight();

        float ax = (W - 1)/(2*sirina);
        float ay = -(H - 1)/(2*visina);

        x = -sirina + x/ax;
        y = visina + y/ay;

        return vec3(0,x,y);
    }

    bool mouse_press_event(int x, int y, int button, int modifiers) {
        std::vector<vec3> tacke;
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            PointsDrawable* pointsDrawable = new PointsDrawable("vertices");

            vec3 tacka = obracun(x, y);
            points.push_back(tacka);

            pointsDrawable->update_vertex_buffer(points);
            pointsDrawable->set_impostor_type(PointsDrawable::SPHERE);
            pointsDrawable->set_default_color(vec3(1.0f, 0, 0));
            pointsDrawable->set_point_size(6.0f);
            this->add_drawable(pointsDrawable);

            update();
        }
        return 1;
    }

    ImprovedViewer(const std::string& title) : Viewer(title) {
        float r = 1000;
        camera()->setType(Camera::ORTHOGRAPHIC);
        camera()->setSceneRadius(r);
        camera()->showEntireScene();
    }

    bool callback_event_keyboard(int key, int action, int modifiers) {
        PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
        LinesDrawable* linesDrawable = new LinesDrawable("line");

        if (action == GLFW_PRESS || action == GLFW_REPEAT) {
            if(key == GLFW_KEY_N) {

            }
            else if(key == GLFW_KEY_K) {
                std::ifstream input;
                PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
                LinesDrawable* linesDrawable = new LinesDrawable("line");


            }
            else if(key == GLFW_KEY_P) {
                MakeSimplePolygon();
            }
            else if(key == GLFW_KEY_I){
                std::vector<vec3> tacke;
                std::ifstream strim;
                strim.open("C:/Users/Korisnik/Desktop/RG/tackeLV2.txt");
                int x, y;
                int brojVrhova;
                if(strim.is_open()) {
                    strim >> brojVrhova;
                    while (strim >> x >> y) {

                    }
                }

            }
        }
        return 1;
    }
    bool Polygon(std::vector<vec3> points){
        LinesDrawable* linesDrawable = new LinesDrawable("linesDrawable");
        std::vector<unsigned int> veze;
        for(int i=0; i<points.size()-1;i++){
            veze.push_back(i);
            veze.push_back(i+1);
        }
        veze.push_back(points.size()-1);
        veze.push_back(0);

        linesDrawable->update_vertex_buffer(points);
        linesDrawable->update_index_buffer(veze);
        linesDrawable->set_line_width(5.0f);
        this->add_drawable(linesDrawable);

        update();
        return 1;
    }

    bool MakeSimplePolygon(){
        int b = 0;
        vec3 Q;
        for(int i=0; i<points.size();i++){
            if((points[i].y < points[b].y) || (points[i].y == points[b].y && points[i].z < points[b].z))
                b = i;
            vec3 Q = points[b];
        }
        auto Cmp = [Q](vec3 R, vec3 S) {
            float a = R.y - Q.y;
            float b = R.z - Q.z;
            float c = S.y - Q.y;
            float d = S.z - Q.z;
            float e = a*d - b*c;
            if(e!=0) return e>0;
            return (pow(a,2) + pow(b,2))<(pow(c,2) + pow(d,2));
        };

        std::sort(points.begin(), points.begin() + points.size(), Cmp);

        return Polygon(points);

    }
};






int main(int argc, char** argv) {
    logging::initialize();

    try {
        ImprovedViewer viewer("LV2");
        viewer.run();

        return EXIT_SUCCESS;

    } catch (const std::runtime_error &e) {
        LOG(ERROR) << "caught a fatal error: " + std::string(e.what());
        return EXIT_FAILURE;
    }

}


