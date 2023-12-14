#ifndef PLUGIN_H
#define PLUGIN_H

#include "Plugins_global.h"
#include <QObject>
#include "PluginInterface.h"
#include <QDebug>
class PlUGINSHAREDDAtASTORAGE_EXPORT PluginSharedDataStorage : public QObject,public PluginUIInterface
{
    Q_OBJECT
    Q_INTERFACES(PluginUIInterface)
    Q_PLUGIN_METADATA(IID "shareddatastorage")
    static int typeId;
public:
    Q_INVOKABLE PluginSharedDataStorage();
    ~PluginSharedDataStorage();
    QString get_name() const
    {
        return "PluginSharedDataStorage";
    }
    QString show_text() const
    {
        return "this is PluginSharedDataStorage";
    }
    virtual void recMsgfromManager(PluginMetaData metaData) Q_DECL_OVERRIDE
    {
        qDebug()<<"插件PluginSharedDataStorage接收到消息："<< metaData.msg;
    }
    QWidget * createWidget() override;

signals:
    void sendMsg2Manager(PluginMetaData) Q_DECL_OVERRIDE;
private:
};

#endif // PLUGIN_H
