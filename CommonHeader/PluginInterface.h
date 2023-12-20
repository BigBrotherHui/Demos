#ifndef PLUGININTERFACE_H
#define PLUGININTERFACE_H

#include <QtPlugin>
#include <QJsonObject>
#include <qwidget.h>

//Q_DECLARE_METATYPE(MsgType);//确保类型可以通过信号槽传递

struct PluginMetaData
{
    QString from;//消息来源
    QString dest;//消息目的地
    QString msg;
    QObject *object = nullptr;
    QJsonObject info = QJsonObject();
    int priority{ 0 };
    bool blocking{ 1 };//阻塞
};
Q_DECLARE_METATYPE(PluginMetaData);//确保类型可以通过信号槽传递
class PluginInterface;
class WidgetBase : public QWidget {
public:
    WidgetBase(QWidget* parent = nullptr) :QWidget(parent){}
    virtual ~WidgetBase(){}
protected:
    PluginInterface* m_interface{ nullptr };
};

class PluginInterface
{
public:
    virtual ~PluginInterface() {}
    virtual QString get_name() const = 0;
    virtual QString show_text() const = 0;
    virtual void recMsgfromManager(PluginMetaData m) = 0//接收到来自创建管理器的消息
    {
        qDebug() << m.msg;
    }
    virtual void sendMsg2Manager(PluginMetaData)   = 0;//给插件管理器发消息
protected:
    enum PluginState {
        Sleep=0,//休眠状态不接收任何信息
        Active
    };
    PluginState m_state=Active;
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
    QWidget* m_widget{ nullptr };
};
Q_DECLARE_INTERFACE(PluginUIInterface, "org.galaxyworld.plugins.PluginUIInterface/1.0")
#endif // PLUGININTERFACE_H
