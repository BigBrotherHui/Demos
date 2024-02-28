#pragma execution_character_set("utf-8")
#ifndef INIUTIL_H
#define INIUTIL_H
#include <QString>
#include <QSettings>

class IniUtil
{
public:
    IniUtil(QString pIniFile, QString pIniTips="配置");
    void Reload();
    QVariant GetValue(QString group, QString key, QVariant defaultValue, bool logging=true);
    QString GetStringValue(QString group, QString key, QString defaultValue, bool logging=true);
    double GetDoubleValue(QString group, QString key, double defaultValue, bool logging=true);
    bool GetBoolValue(QString group, QString key, bool defaultValue=false, bool logging = true);
    int GetIntValue(QString group, QString key, int defaultValue, bool logging=true);
    void SetValue(QString group, QString key, QVariant value, bool logging=true);

private:
    QString DEFAULT_GROUP = "General";
    QString mIniFile,mIniTips;
    QSettings* mSettings=nullptr;
};

#endif // INIUTIL_H
