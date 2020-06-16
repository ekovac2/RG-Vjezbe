/**
* Created by sbecirovic
* Laboratorijska vjezba 2
* Improved viewer
*
*/


#include "improved_viewer.h"
#include <easy3d/viewer/model.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/setting.h>
#include <easy3d/util/logging.h>


using namespace easy3d;

int main(int argc, char** argv) {
    logging::initialize(argv[0]);
    try {
        ImprovedViewer viewer("LV8ImprovedViewer");
        viewer.run();
        return EXIT_SUCCESS;
    } catch (const std::runtime_error &e) {
        LOG(ERROR) << "Caught a fatal error: " + std::string(e.what());
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;

}

