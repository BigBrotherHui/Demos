#pragma once

#include "Plugins_global.h"
#include <QObject>
#include "PluginInterface.h"
#include <QDebug>
class PluginManager;

class PluginCgalDemo : public QObject, public PluginUIInterface
{
    Q_OBJECT
    Q_INTERFACES(PluginUIInterface)
    Q_PLUGIN_METADATA(IID "cgaldemo")

public:
    PluginCgalDemo();
    ~PluginCgalDemo();
    QString get_name() const
    {
        return "PluginCgalDemo";
    }
    QString show_text() const
    {
        return "this is PluginCgalDemo";
    }
    virtual void recMsgfromManager(PluginMetaData metaData) Q_DECL_OVERRIDE
    {
        qDebug() << "插件PluginCgalDemo接收到消息：" << metaData.msg;
    }
    QWidget* createWidget() override;

signals:
    void sendMsg2Manager(PluginMetaData) Q_DECL_OVERRIDE;
private slots:

private:
    QWidget* m_widget{nullptr};
};

