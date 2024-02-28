#ifndef PROGRESSUTIL_H
#define PROGRESSUTIL_H

#include "ProgressWidget.h"

#include <QObject>
#include <functional>


class ProgressUtil:public QObject
{
    Q_OBJECT
public:
    ProgressUtil(QObject* parent=nullptr);
    ~ProgressUtil();
    /**
     * @brief 等待
     * @param pTitle 标题
     * @param pTips 提示
     * @param pCanQuit 是否允许退出，默认不允许
     */
    void Waiting(QString pTitle, QString pTips, bool pCanQuit=false);
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
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void AutoStep(std::function<void()> pCallback, QString pTips);
    /**
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void AutoStep(std::function<void()> pCallback, int pStep, QString pTips);
    /**
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
     * @param pSeconds 预计耗时，单位秒，负数表示1秒
     * @param pInterval 进度更新间隔，单位毫秒，负数表示100毫秒
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void AutoStep(std::function<void()> pCallback, int pStep=-1, int pSeconds=-1, int pInterval=-1, QString pTips="");
    /**
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void* AutoStepWithReturn(std::function<void*()> pCallback, QString pTips);
    /**
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     */
    void* AutoStepWithReturn(std::function<void*()> pCallback, int pStep, QString pTips);
    /**
     * @brief 自动进度处理（使用多线程处理耗时较久的操作）
     * @param pCallback 要调用的操作
     * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
     * @param pSeconds 预计耗时，单位秒，负数表示1秒
     * @param pInterval 进度更新间隔，单位毫秒，负数表示100毫秒
     * @param pTips 进度提示，空表示不修改进度提示，默认为空
     * @return 处理结果
     */
    void* AutoStepWithReturn(std::function<void*()> pCallback, int pStep=-1, int pSeconds=-1, int pInterval=-1, QString pTips="");
    /**
    * @brief 自动进度处理（带图像处理，暂不使用多线程）
    * @param pCallback 要调用的操作
    * @param pTips 进度提示，空表示不修改进度提示，默认为空
    */
    void AutoStepWithImage(std::function<void()> pCallback, QString pTips);
    /**
    * @brief 自动进度处理（带图像处理，暂不使用多线程）
    * @param pCallback 要调用的操作
    * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
    * @param pTips 进度提示，空表示不修改进度提示，默认为空
    */
    void AutoStepWithImage(std::function<void()> pCallback, int pStep, QString pTips);
    /**
    * @brief 自动进度处理（带图像处理，暂不使用多线程）
    * @param pCallback 要调用的操作
    * @param pStep 目标进度值，小于0表示为最大进度值，默认-1
     * @param pSeconds 预计耗时，单位秒，负数表示1秒
     * @param pInterval 进度更新间隔，单位毫秒，负数表示100毫秒
    * @param pTips 进度提示，空表示不修改进度提示，默认为空
    */
    void AutoStepWithImage(std::function<void()> pCallback, int pStep=-1, int pSeconds=-1, int pInterval=-1, QString pTips="");
    /**
     * @brief 自动进度处理结束
     */
    void AutoStepFinish();
    /**
     * @brief 结束进度条
     */
    void Finish();
    /**
     * @brief 是否退出
     */
    bool IsQuit();
signals:
    void signal_autoStepFinish(void* pRet);
    void signal_noname();
private:
    ProgressWidget* mProgressWidget;
    bool mFinish;
    bool mAutoStepFinish;
    int mMax;
    int mMin;
    int mStep;
    void* mAutoStepRet;
};

#endif // PROGRESSUTIL_H
