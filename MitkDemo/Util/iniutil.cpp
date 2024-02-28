#pragma execution_character_set("utf-8")
#include "iniutil.h"
#include <QDebug>

#ifndef tr
#define tr QObject::tr
#endif

/**
 * @brief 否
 */
static const QString BOOL_FALSE = "N";
/**
 * @brief 是
 */
static const QString BOOL_TRUE = "Y";

IniUtil::IniUtil(QString pIniFile, QString pIniTips)
{
    mIniFile = pIniFile;
    mIniTips = pIniTips;
    Reload();
}

void IniUtil::Reload()
{
    if(mSettings!=nullptr)
    {
        delete mSettings;
    }
    mSettings = new QSettings(mIniFile, QSettings::IniFormat);
}

QVariant IniUtil::GetValue(QString group, QString key, QVariant defaultValue, bool logging)
{
    if(mSettings==nullptr)
    {
        if(logging)
        {
            qDebug()<<(tr("获取%1文件信息失败").arg(mIniTips));
        }
        return defaultValue;
    }
    //默认分组(General)被内置，使用默认分组时不能设置分组，否则会被转换为%Genera
    if(group!=DEFAULT_GROUP)
    {
        mSettings->beginGroup(group);
    }
    auto value = mSettings->value(key, defaultValue);
    if(group!=DEFAULT_GROUP)
    {
        mSettings->endGroup();
    }
    if(logging)
    {
        qDebug()<<tr("获取%1信息:%2-%3=%4").arg(mIniTips).arg(group).arg(key).arg(value.toString());
    }
    return value;
}

QString IniUtil::GetStringValue(QString group, QString key, QString defaultValue, bool logging)
{
    auto value = GetValue(group, key, defaultValue, logging);
    if(value==NULL)
    {
        return NULL;
    }
    return value.toString();
}

double IniUtil::GetDoubleValue(QString group, QString key, double defaultValue, bool logging)
{
    auto value = GetValue(group, key, defaultValue, logging);
    if(value==NULL)
    {
        return NULL;
    }
    return value.toDouble();
}

bool IniUtil::GetBoolValue(QString group, QString key, bool defaultValue, bool logging)
{
    auto value = GetStringValue(group, key, defaultValue?BOOL_TRUE:BOOL_FALSE, logging);
    if(value==NULL)
    {
        return NULL;
    }
    return BOOL_TRUE == value.toUpper();
}

int IniUtil::GetIntValue(QString group, QString key, int defaultValue, bool logging)
{
    auto value = GetValue(group, key, defaultValue, logging);
    if(value==NULL)
    {
        return NULL;
    }
    return value.toInt();
}

void IniUtil::SetValue(QString group, QString key, QVariant value, bool logging)
{
    if(mSettings==nullptr){
        if(logging)
        {
            qDebug()<<tr("获取%1文件信息失败").arg(mIniTips);
        }
        return;
    }
    if(value.userType()==QVariant::Bool)
    {
        value = value.toBool()?BOOL_TRUE:BOOL_FALSE;
    }
    //默认分组(General)被内置，使用默认分组时不能设置分组，否则会被转换为%General
    if(group!=DEFAULT_GROUP)
    {
        mSettings->beginGroup(group);
    }
    mSettings->setValue(key, value);
    if(group!=DEFAULT_GROUP)
    {
        mSettings->endGroup();
    }
    Reload();
    if(logging)
    {
        qDebug()<<tr("写%1信息:%2-%3=%4").arg(mIniTips).arg(group).arg(key).arg(value.toString());
    }
}
