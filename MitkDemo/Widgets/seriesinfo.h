#ifndef SERIESINFO_H
#define SERIESINFO_H

#include <QString>

/**
 * @brief 序列信息
 */
struct SeriesTuple
{
    /**
     * @brief 序列id
     */
    int         SeriesId;
    /**
     * @brief 病案id
     */
    int         CaseId;
    /**
     * @brief 访问号
     */
    QString     AccessNo;
    /**
     * @brief 序列号
     */
    QString     SeriesNo;
    /**
     * @brief 图像数量
     */
    int         ImageCount;
    /**
     * @brief 窗宽
     */
    int         WindowWidth;
    /**
     * @brief 窗位
     */
    int         WindowCenter;
};

#endif // SERIESINFO_H
