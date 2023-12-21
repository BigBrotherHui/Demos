#pragma once
#include "Plugins_global.h"
#include <QObject>
#include "PluginInterface.h"
#include <QDebug>
class PLUGINEXAMPLE_EXPORT PluginExample : public QObject,public PluginUIInterface
{
    Q_OBJECT
    Q_INTERFACES(PluginUIInterface)
    Q_PLUGIN_METADATA(IID "example")
    static int typeId;
public:
    Q_INVOKABLE PluginExample();
    ~PluginExample();
    QString get_name() const 
    {
        return "PluginExample";
    }
    QString show_text() const
    {
        return "this is PluginExample";
    }
    virtual void recMsgfromManager(PluginMetaData metaData) Q_DECL_OVERRIDE
    {
        qDebug()<<"插件PluginExample接收到消息："<< metaData.msg;
    }
    QWidget * createWidget() override;
    int m_test;

signals:
    void sendMsg2Manager(PluginMetaData) Q_DECL_OVERRIDE;
private:
};
