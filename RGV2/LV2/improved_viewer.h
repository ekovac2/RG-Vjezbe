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
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/util/timer.h>
#include <easy3d/core/random.h>
#include <easy3d/viewer/drawable_triangles.h>

using namespace easy3d;

class ImprovedViewer : public easy3d::Viewer
{
public:
    struct Triangle{
        vec3 P1, P2, P3;
        Triangle(vec3 P1, vec3 P2, vec3 P3){
            this->P1 = P1;
            this->P2 = P2;
            this->P3 = P3;
        }
    };
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
    struct LineSegment{
        vec3 Q, R;
        LineSegment(vec3 Q, vec3 R){
            this->Q = Q;
            this->R = R;
        }
    };
    enum Events {HLineStart, HLineEnd, VLine};
    struct Event {
        Events t;
        double x;
        int i;

        Event(Events t, double x, int i) {
            this->t = t;
            this->x = x;
            this->i = i;
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
    std::vector<vec3> HVLineSegmentsIntersections(std::vector<LineSegment> s);
    std::vector<vec3> AllLineSegmentsIntersections_Naive(std::vector<LineSegment> s, std::vector<vec3>& colors);
    vec3 LineSegmentsIntersection(LineSegment s1, LineSegment s2);
    bool AreLineSegmentsIntersect(LineSegment s1, LineSegment s2);
    int GeneralizedOrientationOfThreePoints(vec3 P1, vec3 P2, vec3 P3);
    std::vector<Triangle> DelaunayTriangulation_VerySlow(std::vector<vec3> P);

    Timer timer_;
    mutable bool update_;
    mutable std::vector<vec3> points;
    mutable std::vector<vec3> normals;
    TrianglesDrawable* surface = new TrianglesDrawable("faces");
    PointsDrawable *pDrawable = new PointsDrawable("LV8");

    void draw() const override;
    bool callback_event_timer();

protected:
    bool mouse_press_event(int x, int y, int button, int modifiers) override;
    bool key_press_event(int key, int modifiers) override;
    std::string usage() const override;

};


#endif
