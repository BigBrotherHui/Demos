#pragma execution_character_set("utf-8")
#ifndef WINDOWUTIL_H
#define WINDOWUTIL_H

#include <QFileInfoList>
#include <QObject>
#include <QTranslator>

namespace SystemUtil {
/**
 * @brief install event fitler
 * @param parent
 */
void InstallEventFitler(QObject *parent, QObjectList objectList = {});
/**
 * @brief 加载样式
 */
void LoadStyle();
/**
 * @brief 加载本地化
 */
//void LoadTranslator();
/**
 * @brief 检查日志目录是否存在，不存在则新建
 * @param pDir 目录名
 * @return 目录是否存在
 */
bool CheckLogDir(QString pDir);
/**
 * @brief 检查目录是否存在，不存在则新建
 * @param pDir 目录名
 * @param pCreate 不存在时是否新建
 * @return 目录是否存在
 */
bool CheckDir(QString pDir, bool pCreate);
/**
 * @brief 检查目录是否存在，不存在则新建
 * @param pDir 目录名
 * @param pDirDesc 目录描述
 * @param pCreate 不存在时是否新建
 * @return 目录是否存在
 */
bool CheckDir(QString pDir, QString pDirDesc = "", bool pCreate = true);
/**
 * @brief 创建目录
 * @param pDir 目录名
 * @param pDirDesc 目录描述
 * @return 目录创建是否成功
 */
bool CreateDir(QString pDir, QString pDirDesc = "");
/**
 * @brief 删除目录(包括目录下的所有文件及子目录)
 * @param pDir 目录名
 * @param pDirDesc 目录描述
 * @return 目录删除是否成功
 */
bool DeleteDir(QString pDir, QString pDirDesc = "");
/**
 * @brief 获取全路径
 * @param pDir 原路径
 * @return 全路径
 */
QString GetFullDir(QString pDir);
/**
 * @brief 复制文件
 * @param pSourceFile 源文件（全路径）
 * @param pTargetFile 目标文件（全路径）
 * @param pOveride 是否强制覆盖，默认是
 * @return
 */
bool FileCopy(QString pSourceFile, QString pTargetFile, bool pOveride = true);
/**
 * @brief
 * 复制文件夹下所有内容到目标文件夹，包含子目录，目标文件夹不存在时自动创建
 * @param pSourceDir 源文件夹
 * @param pTargetDir 目标文件夹
 * @param pOveride 是否强制覆盖，默认是
 * @return 复制文件数，-1表示失败
 */
int CopyDir(QString pSourceDir, QString pTargetDir, bool pOveride = true);
/**
 * @brief
 * 复制文件夹下指定后缀名的文件到目标文件夹，包含子目录，目标文件夹不存在时自动创建
 * @param pSourceDir 源文件夹
 * @param pTargetDir 目标文件夹
 * @param pSuffix 后缀名
 * @param pOveride 是否强制覆盖，默认是
 * @return 复制文件数，-1表示失败
 */
int CopyDir(QString pSourceDir, QString pTargetDir, QString pSuffix, bool pOveride = true);
/**
 * @brief
 * 复制文件夹下指定文件到目标文件夹，包含子目录，目标文件夹不存在时自动创建
 * @param pSourceDir 源文件夹
 * @param pTargetDir 目标文件夹
 * @param pFilters 过滤器
 * @param pOveride 是否强制覆盖，默认是
 * @return 复制文件数，-1表示失败
 */
int CopyDir(QString pSourceDir, QString pTargetDir, QStringList pFilters, bool pOveride = true);
/**
 * @brief 目录是否为空
 * @param pDir 目录名
 * @return false:目录不为空, true:目录为空或目录不存在
 */
bool IsDirEmpty(QString pDir);
/**
 * @brief 清空目录
 * @param pDir 目录名
 * @param pEmptySub 是否同时删除子目录，默认是
 * @param pSuffix 要删除的文件后缀
 * @return false:清空失败, true:清空成功
 */
bool EmptyDir(QString pDir,
              QString pSuffix, bool pEmptySub = true);
/**
 * @brief 清空目录
 * @param pDir 目录名
 * @param pEmptySub 是否同时删除子目录，默认是
 * @param pNameFilters 名称过滤
 * @param pExcept 要排除的文件
 * @return false:清空失败, true:清空成功
 */
bool EmptyDir(QString pDir, bool pEmptySub = true,
              QStringList pNameFilters = {});
/**
 * @brief 获取文件所在目录
 * @param pFilePath 文件路径
 * @return 文件所在目录
 */
QString GetDir(QString pFilePath);
/**
 * @brief 分割文件所在目录
 * @param pFilePath 文件全路径（分割完成后，不含目录）
 * @return 文件所在目录
 */
QString SplitDir(QString &pFilePath);
/**
 * @brief 获取本机IP
 * @return 本机IP
 */
QStringList GetLocalIp();
/**
 * @brief Round 保留小数位
 * @param pSource 原数值
 * @param pPrecision 精度
 * @return 保留小数位后的数值
 */
double Round(double pSource, int pPrecision);
/**
 * @brief 键值转键名
 * @param pKeyValue 键值
 * @return 键名
 */
QString KeyValueToKeyName(int pKeyValue);
/**
 * @brief 键名转键值
 * @param pKeyName 键名
 * @return 键值
 */
int KeyNameToKeyValue(QString pKeyName);
/**
 * @brief 等待指定时间，不锁线程
 * @param pMS 等待毫秒数
 */
void Wait(int pMS);
/**
 * @brief 以文本形式读取文件内容
 * @param pFile 要读取的文件
 * @return 文件内容
 */
QString ReadFile(QString pFile);
/**
 * @brief 按后缀获取目录下的文件列表
 * @param pDir 目录
 * @param pSuffix 后缀
 * @return 文件列表
 */
QFileInfoList GetFileList(QString pDir, QString pSuffix);
/**
 * @brief 按名称过滤器获取目录下的文件列表
 * @param pDir 目录
 * @param pSuffix 后缀
 * @return 文件列表
 */
QFileInfoList GetFileList(QString pDir, QStringList pNameFilters = {});
};  // namespace SystemUtil

#endif  // WINDOWUTIL_H
