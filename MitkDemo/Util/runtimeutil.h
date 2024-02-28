#ifndef RUNTIMEUTIL_H
#define RUNTIMEUTIL_H

#include <QColor>
#include <QString>

namespace RuntimeUtil
{
/**
 * @brief 获取授权类型
 * @return 授权类型
 */
int GetLicenseType();
/**
 * @brief 获取授权码
 * @return 授权码
 */
QString GetLicense();
/**
 * @brief 设置授权码
 * @param license 授权码
 */
void SetLicense(QString license);
/**
 * @brief 获取序列号
 * @return 序列号
 */
QString GetSerialNo();
/**
 * @brief 设置序列号
 * @param serialNo 序列号
 */
void SetSerialNo(QString serialNo);
/**
 * @brief 获取是否debug模式
 * @return true:开启debug模式
 */
bool GetDebug();
/**
 * @brief 获取是否第一次运行
 * @return true:第一次运行
 */
bool GetFirstRun();
/**
 * @brief 设置非第一次运行
 */
void DisableFirstRun();
/**
 * @brief 获取是否记住登录信息
 * @return true:记住登录信息
 */
bool GetRememberMe();
/**
 * @brief 设置是否记录信息信息
 * @param rememberMe true:记录信息信息
 */
void SetRememberMe(bool rememberMe);
/**
 * @brief 获取上次登录名
 * @return 登录名
 */
QString GetLoginName();
/**
 * @brief 获取上次密码
 * @return 密码
 */
QString GetPassword();
/**
 * @brief 设置上次登录信息
 * @param loginName 登录名
 * @param password 密码
 */
void SetLastLogin(QString loginName, QString password);
/**
 * @brief 设置系统运行状态
 */
void SetRunning(bool pRunning);
/**
 * @brief 获取系统运行状态
 * @return 系统运行状态
 */
bool GetRunning();

};

#endif // RUNTIMEUTIL_H
