#ifndef PLUGIN_H
#define PLUGIN_H

#include "Plugins_global.h"
#include <QObject>
#include "PluginInterface.h"
#include <QDebug>
class PLUGINEXAMPLE_EXPORT Plugin : public QObject,public PluginUIInterface
{
    Q_OBJECT
    Q_INTERFACES(PluginUIInterface)
    Q_PLUGIN_METADATA(IID "example")

public:
    Plugin();
    ~Plugin();
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

signals:
    void sendMsg2Manager(PluginMetaData) Q_DECL_OVERRIDE;
private:
};

#endif // PLUGIN_H
