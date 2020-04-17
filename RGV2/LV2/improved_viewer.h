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

#ifndef RG_LV_IMPROVED_VIEWER_H
#define RG_LV_IMPROVED_VIEWER_H

#include <easy3d/viewer/viewer.h>
using namespace easy3d;

class ImprovedViewer : public easy3d::Viewer
{
public:
    struct Rectangle{
        double x1, x2, y1, y2;
        Rectangle(){}
        Rectangle(double x1, double x2, double y1, double y2){
            this->x1 = x1;
            this->x2 = x2;
            this->y1 = y1;
            this->y2 = y2;
        }
    };
    ImprovedViewer(const std::string& title = "s");
    void readFilePoints(const std::string &ime);
    void readFileLines(const std::string &ime, const std::string &imeIndices);
    bool MakeSimplePolygon();
    bool Polygon(std::vector<vec3> points);
    void readPolygonPoints(const std::string &ime);
    double AreaOfTriangle_Aux(vec3 P1, vec3 P2, vec3 P3);
    bool IsPointInsidePolygon(vec3 P, std::vector<vec3> PI);
    bool ConvexHull_BruteForce_Classic(std::vector<vec3> P);
    bool ConvexHull_BruteForce_VeryBad(std::vector<vec3> P);
    bool IsPointInsideTriangle(vec3 P, std::vector<vec3> T);
    bool ArePointsOnSameSideOfLine(std::vector<vec3> P, vec3 l);
    vec3 LineThroughTwoPoints(vec3 P1, vec3 P2);
    void generatePointsRect();
    std::vector<vec3> RectangularSearch_Naive(std::vector<vec3> P, Rectangle R);


protected:
    bool mouse_press_event(int x, int y, int button, int modifiers) override;
    bool key_press_event(int key, int modifiers) override;
    std::string usage() const override;

private:


};


#endif // RG_LV_IMPROVED_VIEWER_H
