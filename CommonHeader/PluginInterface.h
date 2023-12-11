#ifndef PLUGININTERFACE_H
#define PLUGININTERFACE_H

#include <QtPlugin>
#include <QJsonObject>
#include <qwidget.h>
enum MsgType {
    MSG_SHOWWIDGET
};
Q_DECLARE_METATYPE(MsgType);//确保类型可以通过信号槽传递

struct PluginMetaData
{
    QString from;//消息来源
    QString dest;//消息目的地
    QString msg;
    MsgType msgtype;
    QObject *object = nullptr;
    QJsonObject info = QJsonObject();
};
Q_DECLARE_METATYPE(PluginMetaData);//确保类型可以通过信号槽传递

class WidgetBase : public QWidget {
public:
    WidgetBase(QWidget* parent = nullptr) :QWidget(parent) {}
    virtual ~WidgetBase(){}
protected:
    virtual void processEvent(const PluginMetaData&)=0;
};

class PluginInterface
{
public:
    virtual ~PluginInterface() {}
    virtual QString get_name() const = 0;
    virtual QString show_text() const = 0;
    virtual void recMsgfromManager(PluginMetaData) = 0;//接收到来自创建管理器的消息
    virtual void sendMsg2Manager(PluginMetaData)   = 0;//给插件管理器发消息
};
Q_DECLARE_INTERFACE(PluginInterface, "org.galaxyworld.plugins.PluginInterface/1.0")

class PluginUIInterface : public PluginInterface
{
public:
    virtual ~PluginUIInterface() {
        if (m_widget) {
            delete m_widget;
            m_widget = nullptr;
        }
    }
    virtual QWidget* createWidget() {
        return nullptr;
    }
protected:
    WidgetBase* m_widget{ nullptr };
};
Q_DECLARE_INTERFACE(PluginUIInterface, "org.galaxyworld.plugins.PluginUIInterface/1.0")
#endif // PLUGININTERFACE_H
