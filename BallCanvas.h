#ifndef BALLCANVAS_H
#define BALLCANVAS_H

#include <QWidget>
#include <QTimer>

class BallCanvas : public QWidget {
    Q_OBJECT

public:
    explicit BallCanvas(QWidget *parent = nullptr);

public slots:
    void resetAnimation();
    void setSlopeAngle(int angle);
    void setGravity(int g);
    void setFriction(int f);

protected:
    void paintEvent(QPaintEvent *event) override;

private:
    double ball_x;
    double ball_y;
    double velocity;
    double rotation_angle;
    double gravity;
    double slope_angle;
    double friction;
    QTimer *timer;
};

#endif
