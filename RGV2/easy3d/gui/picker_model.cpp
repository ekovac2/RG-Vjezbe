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


#include <easy3d/gui/picker_model.h>
#include <easy3d/core/point_cloud.h>
#include <easy3d/core/graph.h>
#include <easy3d/core/surface_mesh.h>
#include <easy3d/viewer/opengl_error.h>
#include <easy3d/viewer/camera.h>
#include <easy3d/viewer/drawable_points.h>
#include <easy3d/viewer/drawable_lines.h>
#include <easy3d/viewer/drawable_triangles.h>
#include <easy3d/viewer/framebuffer_object.h>
#include <easy3d/util/logging.h>


namespace easy3d {

    ModelPicker::ModelPicker(const Camera *cam) : Picker(cam) {
        use_gpu_if_supported_ = true;
    }


    ModelPicker::~ModelPicker() {
        if (fbo_) {
            delete fbo_;
            fbo_ = nullptr;
        }
    }


    Model *ModelPicker::pick(const std::vector<Model *> &models, int x, int y) {
        if (models.empty())
            return nullptr;

        int viewport[4];
        glGetIntegerv(GL_VIEWPORT, viewport);
        int width = viewport[2];
        int height = viewport[3];
        setup_framebuffer(width, height);

        //--------------------------------------------------------------------------
        // render the 'scene' into the new FBO.

        // TODO: the performance can be improved. Since the 'scene' is static, we need to render it to the fbo only
        //       once. Then just query. Re-rendering is needed only when the scene is changed/manipulated, or the canvas
        //       size is changed.

        // Bind the offscreen fbo for drawing
        fbo_->bind(); easy3d_debug_log_gl_error; easy3d_debug_log_frame_buffer_error;

        float color[4];
        glGetFloatv(GL_COLOR_CLEAR_VALUE, color);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        states_.clear();
        draw(models);

        // --- Maybe this is not necessary ---------
        glFlush();
        glFinish();
        // -----------------------------------------

        int gl_x, gl_y;
        screen_to_opengl(x, y, gl_x, gl_y, width, height);

        unsigned char c[4];
        fbo_->read_color(c, gl_x, gl_y);

        // switch back to the previous fbo
        fbo_->release(); easy3d_debug_log_gl_error; easy3d_debug_log_frame_buffer_error;

        restore(models);

        // restore the clear color
        glClearColor(color[0], color[1], color[2], color[3]);
        easy3d_debug_log_gl_error;
        easy3d_debug_log_frame_buffer_error;

        //--------------------------------------------------------------------------

        // Convert the color back to an integer ID
        int id = rgb::rgba(c[0], c[1], c[2], c[3]);
        if (id >= 0 && id < models.size()) {
            //LOG(INFO) << "selected model " << models[id]->name();
            return models[id];
        }

        return nullptr;
    }


    // draw the scene
    void ModelPicker::draw(const std::vector<Model *> &models) {
        for (std::size_t i = 0; i < models.size(); ++i) {
            Model *model = models[i];
            if (!model->is_visible())
                continue;

            // the color used to render this model
            int r, g, b, a;
            rgb::encode(static_cast<int>(i), r, g, b, a);
            const vec4 color(r / 255.0f, g / 255.0f, b / 255.0f, a / 255.0f);

            if (dynamic_cast<SurfaceMesh *>(model))
                draw(model->triangles_drawable("faces"), color);
            else if (dynamic_cast<PointCloud *>(model))
                draw(model->points_drawable("vertices"), color);
            else if (dynamic_cast<Graph *>(model)) {
                draw(model->points_drawable("vertices"), color);
                draw(model->lines_drawable("edges"), color);
            }
        }
    }


    // draw a drawable
    void ModelPicker::draw(Drawable *drawable, const vec4 &color) {
        // record
        State state;
        state.lighting = drawable->lighting();
        state.use_texture = drawable->use_texture();
        state.per_vertex_color = drawable->per_vertex_color();
        state.default_color = drawable->default_color();
        states_[drawable] = state;

        // change
        drawable->set_lighting(false);
        drawable->set_use_texture(false);
        drawable->set_per_vertex_color(false);
        drawable->set_default_color(color);

        // render
        drawable->draw(camera(), false);
    }


    void ModelPicker::restore(const std::vector<Model *> &models) {
        for (std::size_t i = 0; i < models.size(); ++i) {
            Model *model = models[i];
            if (!model->is_visible())
                continue;
            if (dynamic_cast<SurfaceMesh *>(model))
                restore(model->triangles_drawable("faces"));
            else if (dynamic_cast<PointCloud *>(model))
                restore(model->points_drawable("vertices"));
            else if (dynamic_cast<Graph *>(model)) {
                restore(model->points_drawable("vertices"));
                restore(model->lines_drawable("edges"));
            }
        }
    }


    void ModelPicker::restore(Drawable *drawable) {
        const State &state = states_[drawable];
        drawable->set_lighting(state.lighting);
        drawable->set_use_texture(state.use_texture);
        drawable->set_per_vertex_color(state.per_vertex_color);
        drawable->set_default_color(state.default_color);
    }

}
