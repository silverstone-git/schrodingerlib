#include "BallCanvas.h"
#include <QPainter>
#include <QtMath>

BallCanvas::BallCanvas(QWidget *parent) : QWidget(parent) {
    gravity = 0.15;
    slope_angle = 30.0;
    friction = 0.0;
    
    timer = new QTimer(this);
    connect(timer, &QTimer::timeout, this, [this]() {
        double angle_rad = qDegreesToRadians(slope_angle);
        const double radius = 10.0;
        
        // Physics: Acceleration with kinetic friction opposing motion
        double accel = gravity * qSin(angle_rad) - friction * gravity * qCos(angle_rad);
        
        // Prevent friction from pushing the ball upwards
        if (accel < 0 && velocity <= 0) {
            accel = 0;
        }
        
        velocity += accel;
        
        // Prevent going backwards if friction is very high
        if (velocity < 0) {
            velocity = 0;
        }
        
        double next_x = ball_x + velocity * qCos(angle_rad);
        double next_y = ball_y + velocity * qSin(angle_rad);
        
        // Stop at the bottom (assuming wedge base is at y=250)
        if (next_y + radius * qCos(angle_rad) > 250.0) {
            timer->stop();
        } else {
            ball_x = next_x;
            ball_y = next_y;
            // Update rotation based on distance traveled to visualize rolling
            rotation_angle += velocity / radius;
        }
        
        update(); // Trigger paintEvent() to redraw
    });
    
    resetAnimation();
}

void BallCanvas::resetAnimation() {
    const double radius = 10.0;
    const double start_x = 50.0;
    const double start_y = 50.0;
    double angle_rad = qDegreesToRadians(slope_angle);
    
    ball_x = start_x + radius * qSin(angle_rad);
    ball_y = start_y - radius * qCos(angle_rad);
    velocity = 0;
    rotation_angle = 0;
    timer->start(16);
    update();
}

void BallCanvas::setSlopeAngle(int angle) { slope_angle = angle; resetAnimation(); }
void BallCanvas::setGravity(int g) { gravity = g / 100.0; resetAnimation(); }
void BallCanvas::setFriction(int f) { friction = f / 100.0; resetAnimation(); }

void BallCanvas::paintEvent(QPaintEvent *) {
    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);

    // 1. Calculate Wedge points based on slope_angle
    double start_x = 50.0;
    double start_y = 50.0;
    double base_y = 250.0;
    double height = base_y - start_y;
    double width = height / qTan(qDegreesToRadians(slope_angle));
    double end_x = start_x + width;

    QPolygonF wedge;
    wedge << QPointF(start_x, start_y) 
          << QPointF(end_x, base_y) 
          << QPointF(start_x, base_y);
    
    painter.setBrush(Qt::lightGray);
    painter.drawPolygon(wedge);

    // 2. Draw the Ball
    painter.translate(ball_x, ball_y);
    painter.rotate(qRadiansToDegrees(rotation_angle));
    
    painter.setBrush(Qt::red);
    painter.setPen(Qt::NoPen);
    painter.drawEllipse(QPointF(0, 0), 10, 10);
    
    // Draw crosshairs to visualize rotation
    painter.setPen(QPen(Qt::black, 2));
    painter.drawLine(-10, 0, 10, 0);
    painter.drawLine(0, -10, 0, 10);
    
    painter.resetTransform();
}
