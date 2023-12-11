#ifndef PluginPolyDataToImage_H
#define PluginPolyDataToImage_H

#include "Plugins_global.h"
#include <QObject>
#include "PluginInterface.h"
#include <QDebug>
class PluginManager;

class PluginPolyDataToImage : public QObject, public PluginUIInterface
{
    Q_OBJECT
    Q_INTERFACES(PluginUIInterface)
    Q_PLUGIN_METADATA(IID "pti")

public:
    PluginPolyDataToImage();
    ~PluginPolyDataToImage();
    QString get_name() const
    {
        return "PluginPolyDataToImage";
    }
    QString show_text() const
    {
        return "this is PluginPolyDataToImage";
    }
    virtual void recMsgfromManager(PluginMetaData metaData) Q_DECL_OVERRIDE
    {
        qDebug() << "插件PluginPolyDataToImage接收到消息：" << metaData.msg;
    }
    QWidget* createWidget() override;

signals:
    void sendMsg2Manager(PluginMetaData) Q_DECL_OVERRIDE;
private slots:

private:
};

#endif // WIDGET_H
