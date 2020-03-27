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

#include <easy3d/core/point_cloud.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/algo/point_cloud_normals.h>
#include <3rd_party/glfw/include/GLFW/glfw3.h>	// for the KEYs
#include <fstream>

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
            "------------------------------------------------ \n");
}

void ImprovedViewer::readFilePoints(const std::string &file)
{

    std::cout<<file;
    PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
    std::ifstream inputFile;
    inputFile.open(file);
    float x,y;
    char comma;
    while(inputFile)
    {
        inputFile >> x >> comma >>  y;
        std::cout<<x;
        vec3 point = camera()->unprojectedCoordinatesOf(vec3(x,y,1));
        point.y=-1;
        points.push_back(point);
    }
    inputFile.close();

    pointsDrawable->update_vertex_buffer(points);
    pointsDrawable->set_default_color(vec4(0.0f,0.0f,0.0f,0.0f));
    pointsDrawable->set_impostor_type(PointsDrawable::SPHERE);
    pointsDrawable->set_point_size(4.0f);
    add_drawable(pointsDrawable);

    camera()->showEntireScene();
    update();
}

void ImprovedViewer::readFileLines(const std::string &fileLines, const std::string &fileIndices)
{

    LinesDrawable * linesDrawable = new LinesDrawable(fileLines);
    std::vector<vec3> linesPoints;
    std::vector<unsigned int> linesIndices;
    std::ifstream inputFile;
    inputFile.open(fileLines);
    float x,y;
    int index;
    char comma;
    while(inputFile)
    {
        inputFile >> x >> comma >>  y;
        vec3 point = camera()->unprojectedCoordinatesOf(vec3(x,y,1));
        point.y=-1;
        linesPoints.push_back(point);
    }
    inputFile.close();
    inputFile.open(fileIndices);
    while(inputFile)
    {
        inputFile >> index >> comma;
        linesIndices.push_back(index);
    }
    inputFile.close();
    linesDrawable->update_vertex_buffer(linesPoints);
    linesDrawable->update_index_buffer(linesIndices);
    linesDrawable->set_default_color(vec4(0.0f,0.0f,0.0f,0.0f));
    add_drawable(linesDrawable);

    camera()->showEntireScene();
    update();
}

void ImprovedViewer::readPolygonPoints(const std::string &file) {
    std::cout<<file;
    std::ifstream inputFile;
    inputFile.open(file);
    float x,y;
    char space;
    int numberOfPoints = 0;
    inputFile >> numberOfPoints;
    int broj = 0;
    while(broj++ != numberOfPoints)
    {
        inputFile >> x >>  y;
        std::cout<<x<< " "<< y<< std::endl;
        vec3 point = camera()->unprojectedCoordinatesOf(vec3(x,y,0));
        point.y=-1;
        points.push_back(point);
    }
    inputFile.close();

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
    add_drawable(linesDrawable);

    camera()->showEntireScene();
    update();

}

bool ImprovedViewer::Polygon(std::vector<vec3> pointys){
    for(int i=0;i<pointys.size();i++){
        pointys[i] = camera()->unprojectedCoordinatesOf(pointys[i]);
    }
    LinesDrawable* linesDrawable = new LinesDrawable("linesDrawable");
    std::vector<unsigned int> veze;
    for(int i=0; i<pointys.size()-1;i++){
        veze.push_back(i);
        veze.push_back(i+1);
    }
    veze.push_back(pointys.size()-1);
    veze.push_back(0);
    std::cout<<" "<<veze;
    linesDrawable->update_vertex_buffer(pointys);
    linesDrawable->update_index_buffer(veze);
    linesDrawable->set_line_width(5.0f);
    add_drawable(linesDrawable);

    camera()->showEntireScene();
    update();
    return 1;
}

std::vector<vec3> pointsOnClick = {};
bool ImprovedViewer::MakeSimplePolygon(){
    std::vector<vec3> pointsOn;
    for(int i=1;i<pointsOnClick.size();i++){
        pointsOn.push_back(camera()->projectedCoordinatesOf(pointsOnClick[i]));
    }
    int b = 0;
    for(int i=0; i<pointsOn.size();i++){
        if((pointsOn[i].x < pointsOn[b].x) || (pointsOn[i].x == pointsOn[b].x && pointsOn[i].y < pointsOn[b].y))
            b = i;
    }
    vec3 Q = pointsOn[b];
    auto Cmp = [Q](vec3 R, vec3 S) {
        float a = R.x - Q.x;
        float b = R.y - Q.y;
        float c = S.x - Q.x;
        float d = S.y - Q.y;
        float e = a*d - b*c;
        if(e!=0) return e>0;
        return (pow(a,2) + pow(b,2))<(pow(c,2) + pow(d,2));
    };

    std::sort(pointsOn.begin(), pointsOn.end(), Cmp);
    return Polygon(pointsOn);

}

double ImprovedViewer::AreaOfTriangle_Aux(vec3 P1, vec3 P2, vec3 P3){
    return (P2.x - P1.x)*(P3.y-P1.y) - (P3.x - P1.x)*(P2.y - P1.y);
}

bool ImprovedViewer::IsPointInsidePolygon(vec3 P, std::vector<vec3> PI) {
    bool f = false;
    PI.push_back(PI[0]);
    for(int i=0; i<PI.size()-1; i++){
        int j = i%PI.size()+1;
        if((PI[j].y <= P.y && P.y < PI[i].y && AreaOfTriangle_Aux(PI[j],PI[i],P) > 0)
            || (PI[i].y <= P.y && P.y < PI[j].y && AreaOfTriangle_Aux(PI[i],PI[j],P) > 0))
            f = !f;
    }
    return f;
}

bool ImprovedViewer::mouse_press_event(int x, int y, int button, int modifiers)  {
    if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
        vec3 point = camera()->unprojectedCoordinatesOf(vec3(x,y,0));

        point.y=-1;
        pointsOnClick.push_back(point);
        colors.push_back(vec3(0.0f,0.0f,0.0f));
        pointsDrawable->update_vertex_buffer(pointsOnClick);
        pointsDrawable->update_color_buffer(colors);
        //pointsDrawable->set_default_color(vec4(0.0f,0.0f,0.0f,0.0f));
        pointsDrawable->set_impostor_type(PointsDrawable::SPHERE);
        pointsDrawable->set_point_size(5.0f);
        add_drawable(pointsDrawable);
        camera()->showEntireScene();
        update();
    }
   else
        return Viewer::mouse_press_event(x,y,button, modifiers);
}

bool ImprovedViewer::key_press_event(int key, int modifiers) {
    bool pritisnutoI = false;
    if (key == GLFW_KEY_N) {
        PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
        std::cout << "\n\nPoint coordinates: \n";
        float x=0,y=0;

        std::cout << "\n x:";
        std::cin >> x;
        std::cout << "\n y:";
        std::cin >> y;
        vec3 point = camera()->unprojectedCoordinatesOf(vec3(x,y,1));
        point.y=-1;
        points.push_back(point);
        std::cout<<point;
        pointsDrawable->update_vertex_buffer(points);
        pointsDrawable->set_default_color(vec4(0.0f,0.0f,0.0f,0.0f));
        pointsDrawable->set_impostor_type(PointsDrawable::SPHERE);
        pointsDrawable->set_point_size(4.0f);

        add_drawable(pointsDrawable);
        camera()->showEntireScene();
        update();
    }
    else if (key == GLFW_KEY_K) {
        std::cout << "Type V to insert vertices, or L to insert lines:\n";
        char sign = ' ';
        std::cin>>sign;
        if(sign == 'V')
        {
            std::string file;
            std::cout<<"Vertices file name: \n";
            std::cin>>file;
            readFilePoints(file);
        }
        else if (sign == 'L'){

            std::string fileVertices, fileLines;
            std::cout<<"Lines vertices file name: \n";
            std::cin>>fileVertices;
            std::cout<<"Lines indices file name:\n";
            std::cin>>fileLines;
            readFileLines(fileVertices, fileLines);
        }
    }
    else if (key == GLFW_KEY_P) {
        MakeSimplePolygon();
    }
    else if (key == GLFW_KEY_I) {
        std::string file;
        std::cout<<"Polygon vertices file name: \n";
        std::cin>>file;
        readPolygonPoints(file);
    }
    else if (key == GLFW_KEY_C){
        vec3 point = camera()->projectedCoordinatesOf(pointsOnClick[pointsOnClick.size()-1]);
        std::vector<vec3> pointsOn = {};
        for(int i=0;i<points.size();i++){
            pointsOn.push_back(camera()->projectedCoordinatesOf(points[i]));
        }
        PointsDrawable* pointsDrawable = new PointsDrawable("vertices");
        bool unutar = IsPointInsidePolygon(point, pointsOn);
        if(unutar){
            std::cout<<"\n"<<point<<std::endl;
            std::cout<<"Da"<<std::endl;
            colors[colors.size()-1] = vec3(0.0f,1.0f,0.0f);
        }
        else{
            std::cout<<"\n"<<point<<std::endl;
            std::cout<<"Ne"<<std::endl;
            colors[colors.size()-1] = vec3(1.0f,0.0f,0.0f);
        }
        pointsDrawable->update_vertex_buffer(pointsOnClick);
        pointsDrawable->update_color_buffer(colors);
        pointsDrawable->set_impostor_type(PointsDrawable::SPHERE);
        pointsDrawable->set_point_size(8.0f);
        add_drawable(pointsDrawable);
        camera()->showEntireScene();
        update();
    }
    else
        return Viewer::key_press_event(key, modifiers);
}





