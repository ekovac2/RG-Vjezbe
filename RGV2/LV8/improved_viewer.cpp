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

//Tacke
std::vector<vec3> points = {};
std::vector<vec3> colors = {};


ImprovedViewer::ImprovedViewer(const std::string& title) : Viewer(title) {
    camera()->setType(Camera::ORTHOGRAPHIC);
    camera()->setUpVector(vec3(0, 0, 1));
    camera()->setViewDirection(vec3(0,1,0));
    camera()->setSceneRadius(camera()->screenHeight());
    camera()->setPosition(vec3(0, 0, 0));
}

std::string ImprovedViewer::usage() const {
    return ("----------- ImprovedViewer usage------------ \n"
            "Right click generates vertex\n"
            "Key press of n enables input\n"
            "Key press of k enables read\n"
            "Key press of p generates simple polygon\n"
            "Key press of c checks if point is in polygon\n"
            "Key press of h created Convex Hull using brute force\n"
            "------------------------------------------------ \n");
}

}












