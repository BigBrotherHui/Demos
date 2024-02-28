#ifndef PROGRESSWIDGET_H
#define PROGRESSWIDGET_H

#include <QWidget>

namespace Ui {
class ProgressWidget;
}

class ProgressWidget : public QWidget
{
    Q_OBJECT

public:
    explicit ProgressWidget(QWidget *parent = nullptr);
    ~ProgressWidget();
    /**
     * @brief 等待
     * @param pTitle 标题
     * @param pTips 提示
     * @param pCanQuit 是否允许退出
     */
    void Waiting(QString pTitle, QString pTips, bool pCanQuit);
    /**
     * @brief  启动进度条
     * @param pTitle 标题
     * @param pMax 最大进度值，默认100
     * @param pMin 最小进度值，默认0
     */
    void Start(QString pTitle, int pMax=100, int pMin=0);
    /**
     * @brief 设置进度提示
     * @param pTips 进度提示
     */
    void SetTips(QString pTips);
    /**
     * @brief 设置进度
     * @param pStep  进度值
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void Step(int pStep, QString pTips="");
    /**
     * @brief 结束进度条
     */
    void Finish();
    /**
     * @brief SetCanQuit 设置是否允许退出
     * @param canQuit 是否允许
     */
    void SetCanQuit(bool canQuit=true);

protected:
    virtual void keyPressEvent(QKeyEvent *event) override;
private:
    Ui::ProgressWidget *ui;
    /**
     * @brief 最小进度值
     */
    int mMin;
    /**
     * @brief 最大进度值
     */
    int mMax;
    /**
     * @brief 当前进度值
     */
    int mStep;
    /**
     * @brief 是否可以退出
     */
    bool mCanQuit;
};

#endif // PROGRESSWIDGET_H
