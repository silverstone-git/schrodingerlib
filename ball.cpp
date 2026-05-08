#include <QApplication>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QLabel>
#include <QSlider>
#include "BallCanvas.h"

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);

    QWidget window;
    QVBoxLayout *layout = new QVBoxLayout(&window);

    BallCanvas *canvas = new BallCanvas();
    layout->addWidget(canvas, 1); // Give canvas stretch priority

    // --- Controls ---
    // Angle Slider
    QHBoxLayout *thetaLayout = new QHBoxLayout();
    QLabel *thetaLabel = new QLabel("Angle (θ):");
    QSlider *thetaSlider = new QSlider(Qt::Horizontal);
    thetaSlider->setRange(10, 80);
    thetaSlider->setValue(30);
    thetaLayout->addWidget(thetaLabel);
    thetaLayout->addWidget(thetaSlider);
    
    // Gravity Slider
    QHBoxLayout *gLayout = new QHBoxLayout();
    QLabel *gLabel = new QLabel("Gravity (g):");
    QSlider *gSlider = new QSlider(Qt::Horizontal);
    gSlider->setRange(1, 100);
    gSlider->setValue(15);
    gLayout->addWidget(gLabel);
    gLayout->addWidget(gSlider);

    // Friction Slider
    QHBoxLayout *fLayout = new QHBoxLayout();
    QLabel *fLabel = new QLabel("Friction (μ):");
    QSlider *fSlider = new QSlider(Qt::Horizontal);
    fSlider->setRange(0, 100);
    fSlider->setValue(0);
    fLayout->addWidget(fLabel);
    fLayout->addWidget(fSlider);

    layout->addLayout(thetaLayout);
    layout->addLayout(gLayout);
    layout->addLayout(fLayout);

    QPushButton *btn = new QPushButton("Rerun Simulation");
    layout->addWidget(btn);

    // Connections
    QObject::connect(btn, &QPushButton::clicked, canvas, &BallCanvas::resetAnimation);
    QObject::connect(thetaSlider, &QSlider::valueChanged, canvas, &BallCanvas::setSlopeAngle);
    QObject::connect(gSlider, &QSlider::valueChanged, canvas, &BallCanvas::setGravity);
    QObject::connect(fSlider, &QSlider::valueChanged, canvas, &BallCanvas::setFriction);

    window.setWindowTitle("Rolling Ball Simulation");
    window.resize(500, 500);
    window.show();

    return app.exec();
}
