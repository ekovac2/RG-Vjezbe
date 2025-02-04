#ifndef PAINT_CANVAS_H
#define PAINT_CANVAS_H

#include <easy3d/viewer/opengl.h>
#include <easy3d/core/types.h>

#include <QOpenGLWidget>

#include "tools/canvas.h"


namespace easy3d {
    class Camera;
    class Model;
    class TrianglesDrawable;
    class AmbientOcclusion;
    class Shadow;
    class Transparency;
    class EyeDomeLighting;
    class OpenGLTimer;
}

class QWidget;
class QOpenGLFunctions;

class PaintCanvas : public QOpenGLWidget, public easy3d::Canvas
{
       Q_OBJECT
public:
    explicit PaintCanvas(QWidget* parent = nullptr);
    virtual ~PaintCanvas() override;

    virtual std::string usage() const;

	// the actual samples received
	int samples() const { return samples_; }

    // Scaling factor for high DPI devices
    double dpi_scaling() const { return dpi_scaling_; }

    const easy3d::vec4& backGroundColor() const { return background_color_; }
    void setBackgroundColor(const easy3d::vec4& c);

    void addModel(easy3d::Model* model, bool create_default_drawables = true);
	void deleteModel(easy3d::Model* model);

	const std::vector<easy3d::Model*>& models() const override { return models_; }
    easy3d::Model* currentModel();

    // the camera
    easy3d::Camera* camera() override { return camera_; }
    const easy3d::Camera* camera() const override { return camera_; }

	// moves the camera so that the 'model' is centered on the screen.
	// if 'model' is NULL, it centers the entire scene (all models).
	void fitScreen(const easy3d::Model* model = nullptr);

    // Returns the coordinates of the 3D point located at pixel (x,y) on screen.
    // x, y: screen point expressed in pixel units with an origin in the upper left corner.
    // found: indicates whether a point was found or not.
    // NOTE: This method assumes that a GL context is available, and that its
    //		 content was drawn using the Camera (i.e. using its projection and modelview
    //		 matrices). This method hence cannot be used for offscreen Camera computations.
    //		 Use cameraCoordinatesOf() and worldCoordinatesOf() to perform similar
    //		 operations in that case.
    //       The precision of the z-Buffer highly depends on how the zNear() and zFar()
    //       values are fitted to your scene. Loose boundaries will result in imprecision
    //		 along the viewing direction.
    easy3d::vec3 pointUnderPixel(const QPoint& p, bool &found) const;

	// Save snapshot. This function has no limit on the image size.
    // file_name: the image file name
	// w, h: the width and height of the require snapshot
    // bk_color: the background color
    // expand: expand the frustum to ensure the image aspect ratio
	bool saveSnapshot(int w, int h, int samples, const QString& file_name, bool bk_white = true, bool expand = true);

    easy3d::AmbientOcclusion* ssao();

    easy3d::Shadow* shadow();
    void enableShadow(bool b);

    easy3d::Transparency* transparency();
    void enableTransparency(bool b);

    easy3d::EyeDomeLighting* edl();
    void enableEyeDomeLighting(bool b);

public slots:
    void invertSelection();
    void deleteSelectedPrimitives();

    void copyCamera();
    void pasteCamera();

protected:

    /* Set up required OpenGL resources/state and then calls user-defined init().
     * This method is called once before the first call to paintGL() or resizeGL().
     * Note:
     *  - Overload init() instead of this method to modify specific OpenGL state;
     *  - The framebuffer is not yet available at this stage.
     */
    virtual void initializeGL() override;

    /* User-defined initialization method.
     * This method is called within initializeGL() and should be overloaded to
     * initialize OpenGL flags/resources, e.g.,
     *  - OpenGL state modification;
     *  - shader program/texture/VAOs creation;
     *  - camera initialization;
     *  - previous viewer state restoration;
     *  - ...
     * All OpenGL specific initializations must be done in this method.
     * OpenGL context is not yet available in your viewer constructor.
     * NOTE:
     *  - If you derive you own viewer from this class, don't forget to call
     *    Viewer::init() at the beginning of your inherited function.
     *  - Do not call updateGL() in this method (resulting in an infinite loop).
     */
    virtual void init();

    /* Sets up the OpenGL viewport, projection, etc. Gets called whenever the 
	 * widget has been resized (and also when it is shown for the first time 
	 * because all newly created widgets get a resize event automatically).
     * If you overload this method, first call the inherited method in which
     * the projection matrix is updated.
     */
    virtual void resizeGL(int width, int height) override;

    /* Renders the OpenGL scene. Gets called whenever the widget needs to 
	 * be updated. Internally, it calls the following methods in order:
     *  - preDraw(): places the camera in the world coordinate system;
     *  - draw(): main drawing method. Should be overloaded.
     *  - postDraw(): display of visual hints (world axis, FPS...)
     * Note: For normal rendering, i.e., drawing triggered by the
     *       paintEvent(), the clearing of the color and depth buffers is
     *       done by the widget before entering paintGL(). However, if you
     *       want to reuse the paintGL() method for offscreen rendering,
     *       you have to clear both buffers before calling paintGL().
     */
    virtual void paintGL() override;

    /* This function will be called before the main draw procedure.
     */
    virtual void preDraw();

    /* The core method of the viewer, that draws the scene.
     */
    virtual void draw();

    /* Called after draw() to draw viewer visual hints.
     * By default, it displays axis and visual hints if the respective flags are set.
     */
    virtual void postDraw();

    // OpenGL resources (e.g., shaders, textures, VAOs) must destroyed when
    // there exists a valid rendering context. It is (usually) a bad idea to
    // clean up OpenGL in a destructor because the OpenGL context may not exist
    // (e.g., destroyed already) or the visible one is not *current*. This
    // cleanup() function is to ensure you have a valid rendering context.
    // See also init().
    // NOTE: Don't forget to call Viewer::cleanup() at the end of your
    //		 inherited function.
    virtual void cleanup();

protected:
    virtual void mousePressEvent(QMouseEvent *) override;    // Mouse button press event handler
    virtual void mouseMoveEvent(QMouseEvent *) override;
    virtual void mouseReleaseEvent(QMouseEvent *) override;  // Mouse button release event handler
    virtual void mouseDoubleClickEvent(QMouseEvent *) override;
    virtual void wheelEvent(QWheelEvent *) override;         // Mouse scroll event handler
    virtual void keyPressEvent(QKeyEvent *) override;        // Keyboard press event handler.
    virtual void keyReleaseEvent(QKeyEvent *) override;      // Keyboard press event handler.
    virtual void timerEvent(QTimerEvent *) override;
    virtual void closeEvent(QCloseEvent *) override;

protected:
    // Create default drawables for rendering
    //  - for point clouds, it creates a PointsDrawable
    //      - per point color will be enabled if vertex property 'v:color' exists
    //      - normals (stored in vertex property 'v:normal') will be used for rendering if exists.
    //  - for mesh surfaces, it creates a TrianglesDrawable
    //      - per vertex color will be enabled if vertex property 'v:color' exists
    // TODO: move this function to Renderer module; enable per face color for surface meshes.
    void create_drawables(easy3d::Model* model);

    void drawCornerAxes();

signals:
    void currentModelChanged();

protected:
	// Actually I can inherit the viewer from QOpenGLFunctions (thus no such a member 
	// variable). Having it as a member can eliminate including the header file.
	QOpenGLFunctions* func_;
    easy3d::OpenGLTimer* gpu_timer_;
    double gpu_time_;

    double  dpi_scaling_;
    int     samples_;

    easy3d::Camera*	camera_;
    easy3d::vec4	background_color_;

    Qt::MouseButton pressed_button_;
    QPoint  mouse_previous_pos_;

    bool    show_pivot_point_;

    //----------------- viewer data -------------------

    // corner axes
    easy3d::TrianglesDrawable* drawable_axes_;

    // camera path
    bool	show_camera_path_;

    std::vector<easy3d::Model*> models_;
    int model_idx_;

    //----------------- filters -------------------

    easy3d::AmbientOcclusion* ssao_;

    easy3d::Transparency* transparency_;
    bool    transparency_enabled_;

    easy3d::Shadow* shadow_;
    bool    shadow_enabled_;

    easy3d::EyeDomeLighting* edl_;
    bool    edl_enabled_;
};


#endif // PAINT_CANVAS_H
