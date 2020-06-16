/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "improved_viewer.h"

#include <easy3d/fileio/resources.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/algo/point_cloud_normals.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>	// for the KEYs
#include <fstream>
#include <easy3d/util/timer.h>
#include <easy3d/core/random.h>

using namespace easy3d;

ImprovedViewer::ImprovedViewer(const std::string& title) : Viewer(title) {
    camera()->setType(easy3d::Camera::ORTHOGRAPHIC);
    camera()->setUpVector(vec3(0, 0, 1));
    camera()->setViewDirection(vec3(0, 0, 0));
    camera()->setPosition(vec3(0, 0, 0));
    camera()->setSceneRadius(2.0f);
    camera_->showEntireScene();
}

std::string ImprovedViewer::usage() const {
    return ("----------- ImprovedViewer usage------------ \n"
            "Right click generates vertex\n"
            "------------------------------------------------ \n");
}

bool ImprovedViewer::callback_event_timer(){
    update_ = true;
    update();
    return true;
}

void ImprovedViewer::draw() const {
    if(update_) {
        double q3 = std::sin(quarter_pi / 40);
        double q4 = std::cos(quarter_pi / 40);
        quat quat1 = quat(0,0, q3, q4);
        std::vector<vec3> tacke;
        for(int i = 0; i < points.size(); i++){
            tacke.push_back(quat1.rotate(points[i]));
        }
        points = tacke;
        surface->update_vertex_buffer(points);
        normals.clear();
        surface->update_normal_buffer(normals);
    }
    update();
    update_ = false;
    Viewer::draw();
}

bool ImprovedViewer::key_press_event(int key, int modifiers) {
    if (key == GLFW_KEY_A){
        const std::vector<vec3>& vertices = resource::bunny_vertices;
        const std::vector<int>& indices = resource::bunny_indices;

        for (std::size_t i=0; i<indices.size(); i+=3) {
            const vec3& a = vertices[ indices[i]	 ];		points.push_back(a);
            const vec3& b = vertices[ indices[i + 1] ]; 	points.push_back(b);
            const vec3& c = vertices[ indices[i + 2] ];		points.push_back(c);
            const vec3& n = geom::triangle_normal(a, b, c);
            normals.insert(normals.end(), 3, n);
        }
        surface->update_vertex_buffer(points);
        surface->update_normal_buffer(normals);

        add_drawable(surface);
        update();
        timer_.set_interval(50, &ImprovedViewer::callback_event_timer,this);
    }

return true;

}














