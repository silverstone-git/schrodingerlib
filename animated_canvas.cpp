#include <QWidget>
#include <QPainter>
#include <QTimer>
#include <QPushButton>
#include <QVBoxLayout>

class AnimatedCanvas : public QWidget {
    Q_OBJECT
    int x_pos = 0;
    QTimer *timer;

public:
    AnimatedCanvas(QWidget *parent = nullptr) : QWidget(parent) {
        timer = new QTimer(this);
        // Connect timer to our update function
        connect(timer, &QTimer::timeout, this, [=]() {
            x_pos += 2;      // Move 2 pixels per frame
            update();        // Tell Qt to call paintEvent()
        });
        timer->start(16);    // Run at ~60 frames per second
    }

    void resetAnimation() {
        x_pos = 0;           // Reset position
        timer->start(16);    // Restart timer if it was stopped
    }

protected:
    void paintEvent(QPaintEvent *) override {
        QPainter painter(this);
        painter.setBrush(Qt::blue);
        painter.drawEllipse(x_pos, 100, 50, 50); // Draw at current x_pos
    }
};

