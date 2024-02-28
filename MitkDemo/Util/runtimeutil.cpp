#include "runtimeutil.h"
#include "iniutil.h"

namespace RuntimeUtil
{
/**
 * @brief 文件
 */
static const QString RUNTIME_FILE = "./runtime.ini";
/**
 * @brief 否
 */
static const QString BOOL_FALSE = "N";
/**
 * @brief 是
 */
static const QString BOOL_TRUE = "Y";
/**
 * @brief 授权组
 */
static const QString LICENSE_GROUP = "License";
/**
 * @brief 授权类型
 */
static const QString LICENSE_TYPE_KEY = "type";
/**
 * @brief 授权码
 */
static const QString LICENSE_KEY = "license";
/**
 * @brief 序列号
 */
static const QString SERIAL_NO_KEY = "serialNo";
/**
 * @brief 系统组
 */
static const QString SYSTEM_GROUP = "System";
/**
 * @brief DEBUG模式
 */
static const QString DEBUG_KEY = "debug";
/**
 * @brief 是否第一次运行
 */
static const QString FIRST_RUN_KEY = "firstRun";
/**
 * @brief 登录组
 */
static const QString LOGIN_GROUP = "Login";
/**
 * @brief 是否记住登录信息
 */
static const QString REMEMBER_ME_KEY = "rememberMe";
/**
 * @brief 是否正在运行
 */
static const QString RUNNING_KEY = "running";
/**
 * @brief 登录名
 */
static const QString LOGIN_NAME_KEY = "loginName";
/**
 * @brief 密码
 */
static const QString PASSWORD_KEY = "password";

IniUtil *mIniUtil = new IniUtil(RUNTIME_FILE);

QString getStringConfig(QString group, QString key, QString defaultValue="", bool logging=true){
    return mIniUtil->GetStringValue(group, key, defaultValue, logging);
}

int getIntConfig(QString group, QString key, int defaultValue=0, bool logging=true){
    return mIniUtil->GetIntValue(group, key, defaultValue, logging);
}

double getDoubleConfig(QString group, QString key, double defaultValue=0, bool logging=true){
    return mIniUtil->GetDoubleValue(group, key, defaultValue, logging);
}

void SetConfig(QString group, QString key, QVariant value){
    mIniUtil->SetValue(group, key, value);
}

int GetLicenseType(){
    return getIntConfig(LICENSE_GROUP, LICENSE_TYPE_KEY);
}

QString GetLicense(){
    return getStringConfig(LICENSE_GROUP, LICENSE_KEY);
}

void SetLicense(QString license){
    SetConfig(LICENSE_GROUP, LICENSE_KEY, license);
}

QString GetSerialNo()
{
    return getStringConfig(LICENSE_GROUP, SERIAL_NO_KEY);
}

void SetSerialNo(QString serialNo){
    QString serialNo_ = getStringConfig(LICENSE_GROUP, SERIAL_NO_KEY);
    if(serialNo==serialNo_){
        return;
    }
    SetConfig(LICENSE_GROUP, SERIAL_NO_KEY, serialNo);
}

bool GetDebug(){
    QString debugStr = getStringConfig(SYSTEM_GROUP, DEBUG_KEY);
    return "1"==debugStr || "YES"==debugStr.toUpper() || "Y"==debugStr.toUpper();
}

bool GetFirstRun()
{
    QString firstRunStr = getStringConfig(SYSTEM_GROUP, FIRST_RUN_KEY);
    return BOOL_FALSE!=firstRunStr.toUpper();
}

void DisableFirstRun()
{
    SetConfig(SYSTEM_GROUP, FIRST_RUN_KEY, BOOL_FALSE);
}

bool GetRememberMe()
{
    QString rememberMe = getStringConfig(LOGIN_GROUP, REMEMBER_ME_KEY, BOOL_FALSE);
    return BOOL_TRUE==rememberMe.toUpper();
}

void SetRememberMe(bool rememberMe)
{
    SetConfig(LOGIN_GROUP, REMEMBER_ME_KEY, rememberMe?BOOL_TRUE:BOOL_FALSE);
}

QString GetLoginName()
{
    return getStringConfig(LOGIN_GROUP, LOGIN_NAME_KEY);
}

QString GetPassword()
{
    return getStringConfig(LOGIN_GROUP, PASSWORD_KEY);
}

void SetLastLogin(QString loginName, QString password)
{
    SetConfig(LOGIN_GROUP, LOGIN_NAME_KEY, loginName);
    SetConfig(LOGIN_GROUP, PASSWORD_KEY, password);
}

void SetRunning(bool pRunning)
{
    SetConfig(SYSTEM_GROUP, RUNNING_KEY, pRunning?BOOL_TRUE:BOOL_FALSE);
}

bool GetRunning()
{
    QString running = getStringConfig(SYSTEM_GROUP, RUNNING_KEY, BOOL_FALSE);
    return BOOL_TRUE==running.toUpper();

}

}
